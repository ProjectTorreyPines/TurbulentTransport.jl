#!/usr/bin/env python3
# coding: utf-8
"""Combine every .onnx in the working directory into one TorchScript committee.

Regular models: mean and variance of the members, then output affine
(ym, ysigma) so the public interface is in TGLF units.

Residual models (directory name ``parent_suffix``, parent xnames+ynames
match local xnames): load ``../parent/committee.pt``, concatenate its
mean onto the public inputs, run the local members as error factors,
return ``err * parent_mean`` and ``err**2 * parent_var``.

Verification (default) is against ONNX Runtime on the member files, not
Julia/Flux. The ONNX files are the export snapshot; this script only
needs to prove ConvertModel + committee wrap + trace match that snapshot.
Pass --no-verify to skip.

Requires xm.txt (and xsigma.txt; ym.txt/ysigma.txt for regular models).
"""

from __future__ import annotations

import argparse
import os
import re
import sys

import numpy as np
import onnx
import torch
import torch.nn as nn
from onnx2pytorch import ConvertModel


class FixedClip(nn.Module):
    """onnx2pytorch Clip can be called with extra args; keep min/max only."""

    def __init__(self, old):
        super().__init__()
        self.old = old

    def forward(self, *args):
        if len(args) > 2:
            args = args[:2]
        return self.old(*args)


def patch_clips(model):
    for name, module in list(model.named_children()):
        if module.__class__.__name__ == "Clip":
            setattr(model, name, FixedClip(module))
        else:
            patch_clips(module)


def _die(message):
    print(f"error: {message}", file=sys.stderr)
    raise SystemExit(1)


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--no-verify",
        action="store_true",
        help="Skip ONNX Runtime checks of members and of committee.pt",
    )
    p.add_argument("--atol", type=float, default=1.0e-4)
    p.add_argument("--rtol", type=float, default=1.0e-4)
    p.add_argument(
        "--n-random",
        type=int,
        default=8,
        help="Extra Gaussian samples around xm used in the checks",
    )
    return p.parse_args(argv)


def list_onnx(directory="."):
    files = sorted(f for f in os.listdir(directory) if f.endswith(".onnx"))
    if not files:
        _die(f"no .onnx files in {os.path.abspath(directory)}")
    return [os.path.join(directory, f) for f in files]


def load_txt(path):
    return np.loadtxt(path).astype(np.float32)


def to_numpy(y):
    if isinstance(y, (tuple, list)):
        y = y[0]
    if isinstance(y, torch.Tensor):
        return y.detach().cpu().numpy()
    return np.asarray(y)


def module_numpy(mod, x):
    xt = torch.as_tensor(np.ascontiguousarray(x, dtype=np.float32))
    with torch.no_grad():
        return to_numpy(mod(xt))


def convert_member(path):
    onnx_model = onnx.load(path)
    model = ConvertModel(onnx_model)
    patch_clips(model)
    model.eval()
    return model


class Committee(nn.Module):
    def __init__(self, models, xm, xs, ym, ys):
        super().__init__()
        self.models = nn.ModuleList(models)
        self.register_buffer("xm", torch.as_tensor(xm))
        self.register_buffer("xs", torch.as_tensor(xs))
        self.register_buffer("ym", torch.as_tensor(ym))
        self.register_buffer("ys", torch.as_tensor(ys))

    def forward(self, inp):
        scaled = (inp - self.xm) / self.xs
        stacked = torch.stack([m(scaled) for m in self.models])
        mean = torch.mean(stacked, dim=0)
        var = torch.var(stacked, dim=0)
        return self.ys * mean + self.ym, self.ys**2 * var


class ResidualCommittee(nn.Module):
    def __init__(self, models, parent, xm, xs):
        super().__init__()
        self.models = nn.ModuleList(models)
        self.parent = parent
        self.register_buffer("xm", torch.as_tensor(xm))
        self.register_buffer("xs", torch.as_tensor(xs))

    def forward(self, inp):
        parent_mean, parent_var = self.parent(inp)
        local = torch.cat([inp, parent_mean], dim=1)
        scaled = (local - self.xm) / self.xs
        stacked = torch.stack([m(scaled) for m in self.models])
        err_mean = torch.mean(stacked, dim=0)
        return err_mean * parent_mean, err_mean**2 * parent_var


def detect_residual(dir_name):
    """Return (is_residual, parent_dir, parent_committee or None, parent_ynames)."""
    if "_" not in dir_name or re.search(r"-\d+$", dir_name):
        return False, None, None, None
    parent_name = dir_name.rsplit("_", 1)[0]
    parent_dir = os.path.join("..", parent_name)
    if not os.path.isdir(parent_dir):
        print("✗ Parent directory not found → treating as regular model")
        return False, None, None, None
    try:
        parent_xnames = open(os.path.join(parent_dir, "xnames.txt")).read().split()
        parent_ynames = open(os.path.join(parent_dir, "ynames.txt")).read().split()
        current_xnames = open("xnames.txt").read().split()
    except OSError as exc:
        print(f"✗ Could not read xnames files: {exc} → treating as regular model")
        return False, None, None, None
    if parent_xnames + parent_ynames != current_xnames:
        print("✗ xnames check failed → treating as regular model")
        return False, None, None, None
    committee_file = os.path.join(parent_dir, "committee.pt")
    if not os.path.exists(committee_file):
        print(
            f"⚠ Residual model detected (parent: {parent_name}), "
            "but committee.pt not found."
        )
        print("   Please run this script in the parent directory first.")
        return False, None, None, None
    parent_committee = torch.jit.load(committee_file)
    parent_committee.eval()
    print(f"✓ Residual model detected. Parent: {parent_name}")
    return True, parent_dir, parent_committee, parent_ynames


def make_physical_inputs(xm, xs, n_random, seed=0):
    xm = np.asarray(xm, dtype=np.float32).reshape(-1)
    xs = np.asarray(xs, dtype=np.float32).reshape(-1)
    rows = [xm.copy()]
    if n_random > 0:
        rng = np.random.default_rng(seed)
        noise = rng.standard_normal((n_random, xm.size)).astype(np.float32)
        rows.append(xm + 0.5 * xs * noise)
    return np.vstack(rows)


def require_ort():
    try:
        import onnxruntime as ort
    except ImportError:
        _die(
            "onnxruntime is required to verify against the member ONNX files "
            "(pass --no-verify to skip)"
        )
    return ort


def run_onnx_members(ort, paths, scaled):
    scaled = np.ascontiguousarray(scaled, dtype=np.float32)
    outs = []
    for path in paths:
        sess = ort.InferenceSession(path, providers=["CPUExecutionProvider"])
        name = sess.get_inputs()[0].name
        outs.append(sess.run(None, {name: scaled})[0])
    return np.stack(outs, axis=0)


def report_diff(label, got, ref, atol, rtol):
    got = np.asarray(got, dtype=np.float64)
    ref = np.asarray(ref, dtype=np.float64)
    if got.shape != ref.shape:
        _die(f"{label}: shape {got.shape} != {ref.shape}")
    abs_err = np.max(np.abs(got - ref))
    rel_err = np.max(np.abs(got - ref) / (np.abs(ref) + 1.0e-8))
    ok = np.allclose(got, ref, atol=atol, rtol=rtol, equal_nan=False)
    status = "ok" if ok else "FAIL"
    print(f"  {status} {label}: max abs {abs_err:.3e}  max rel {rel_err:.3e}")
    if not ok:
        _die(
            f"{label} disagrees with ONNX Runtime "
            f"(atol={atol}, rtol={rtol})"
        )


def numpy_regular_committee(members, ym, ys):
    mean = members.mean(axis=0)
    var = members.var(axis=0, ddof=1)
    return ys * mean + ym, (ys**2) * var


def numpy_residual_committee(err_members, parent_mean, parent_var):
    err_mean = err_members.mean(axis=0)
    return err_mean * parent_mean, (err_mean**2) * parent_var


def verify_members(ort, paths, models, scaled, atol, rtol):
    print(f"Checking {len(paths)} members vs ONNX Runtime, batch={scaled.shape[0]}")
    stacked_ort = run_onnx_members(ort, paths, scaled)
    stacked_pt = []
    for path, model, y_ort in zip(paths, models, stacked_ort):
        y_pt = module_numpy(model, scaled)
        report_diff(os.path.basename(path), y_pt, y_ort, atol, rtol)
        stacked_pt.append(y_pt)
    return stacked_ort, np.stack(stacked_pt, axis=0)


def verify_committee_pair(label, mean_got, var_got, mean_ref, var_ref, atol, rtol):
    report_diff(f"{label} mean", mean_got, mean_ref, atol, rtol)
    report_diff(f"{label} var", var_got, var_ref, atol, rtol)


def main(argv=None):
    args = parse_args(argv)
    dir_name = os.path.basename(os.getcwd())
    print(f"Processing: {dir_name}")

    is_residual, parent_dir, parent_committee, parent_ynames = detect_residual(
        dir_name
    )
    if not is_residual:
        print("✓ Treating as regular model")

    xm = load_txt("xm.txt")
    xs = load_txt("xsigma.txt")
    if not is_residual:
        ym = load_txt("ym.txt")
        ys = load_txt("ysigma.txt")

    if is_residual:
        n_parent_out = len(parent_ynames)
        public_xm = xm[: len(xm) - n_parent_out]
        public_xs = xs[: len(xs) - n_parent_out]
        example_input = torch.as_tensor(
            public_xm.reshape(1, -1), dtype=torch.float32
        )
    else:
        public_xm, public_xs = xm, xs
        example_input = torch.as_tensor(xm.reshape(1, -1), dtype=torch.float32)

    onnx_paths = list_onnx(".")
    print(f"Converting {len(onnx_paths)} ONNX members")
    pytorch_models = []
    for path in onnx_paths:
        print(f"  {os.path.basename(path)}")
        pytorch_models.append(convert_member(path))

    if is_residual:
        committee = ResidualCommittee(
            pytorch_models, parent_committee, xm, xs
        )
    else:
        committee = Committee(pytorch_models, xm, xs, ym, ys)
    committee.eval()

    if not args.no_verify:
        ort = require_ort()
        x_phys = make_physical_inputs(public_xm, public_xs, args.n_random)
        if is_residual:
            parent_paths = list_onnx(parent_dir)
            parent_xm = load_txt(os.path.join(parent_dir, "xm.txt"))
            parent_xs = load_txt(os.path.join(parent_dir, "xsigma.txt"))
            parent_ym = load_txt(os.path.join(parent_dir, "ym.txt"))
            parent_ys = load_txt(os.path.join(parent_dir, "ysigma.txt"))
            parent_scaled = (x_phys - parent_xm) / parent_xs
            parent_ort = run_onnx_members(ort, parent_paths, parent_scaled)
            parent_mean_ref, parent_var_ref = numpy_regular_committee(
                parent_ort, parent_ym, parent_ys
            )
            local = np.concatenate([x_phys, parent_mean_ref], axis=1)
            scaled = (local - xm) / xs
            stacked_ort, _ = verify_members(
                ort, onnx_paths, pytorch_models, scaled, args.atol, args.rtol
            )
            mean_ref, var_ref = numpy_residual_committee(
                stacked_ort, parent_mean_ref, parent_var_ref
            )
        else:
            scaled = (x_phys - xm) / xs
            stacked_ort, _ = verify_members(
                ort, onnx_paths, pytorch_models, scaled, args.atol, args.rtol
            )
            mean_ref, var_ref = numpy_regular_committee(stacked_ort, ym, ys)

        with torch.no_grad():
            eager_mean, eager_var = committee(torch.as_tensor(x_phys))
        eager_mean = to_numpy(eager_mean)
        eager_var = to_numpy(eager_var)
        print("Checking eager committee vs ONNX Runtime composition")
        verify_committee_pair(
            "eager", eager_mean, eager_var, mean_ref, var_ref, args.atol, args.rtol
        )

    print("Tracing")
    traced = torch.jit.trace(committee, example_input)
    traced.eval()

    if not args.no_verify:
        with torch.no_grad():
            tr_mean, tr_var = traced(torch.as_tensor(x_phys))
        tr_mean, tr_var = to_numpy(tr_mean), to_numpy(tr_var)
        print("Checking traced committee vs ONNX Runtime composition")
        verify_committee_pair(
            "traced", tr_mean, tr_var, mean_ref, var_ref, args.atol, args.rtol
        )
        print("Checking traced committee vs eager")
        verify_committee_pair(
            "trace-vs-eager",
            tr_mean,
            tr_var,
            eager_mean,
            eager_var,
            args.atol,
            args.rtol,
        )

    out = "committee.pt"
    traced.save(out)
    print(f"Committee TorchScript saved to {out}")
    if is_residual:
        print("This is a residual model")
    else:
        print("This is a regular model")


if __name__ == "__main__":
    main()
