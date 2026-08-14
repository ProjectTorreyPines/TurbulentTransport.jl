[![CI](https://github.com/ProjectTorreyPines/TurbulentTransport.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/ProjectTorreyPines/TurbulentTransport.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/github/ProjectTorreyPines/TurbulentTransport.jl/graph/badge.svg?token=6CdJLykYa9)](https://codecov.io/github/ProjectTorreyPines/TurbulentTransport.jl)

# TurbulentTransport.jl

TurbulentTransport.jl runs tokamak turbulent-transport models from a common
interface: the first-principles codes **TGLF**, **TJLF** (the Julia
reimplementation), and **QLGYRO**, and a family of fast **neural-network
surrogates** (TGLF-NN, GKNN, QLNN, FINN, ModeID-NN) trained to reproduce them in
milliseconds. It builds the model inputs directly from IMAS data structures, so
it plugs straight into FUSE and related ProjectTorreyPines tools.

If you are new here, the fastest way to get oriented is the
[example notebooks](#examples-and-tutorials) below — start with
[`examples/run_TGLFNN.ipynb`](https://github.com/ProjectTorreyPines/TurbulentTransport.jl/blob/master/examples/run_TGLFNN.ipynb).
For the full model list and API reference, see the
[online documentation](https://projecttorreypines.github.io/TurbulentTransport.jl/dev).

## Quick start

Add the package (it lives in the ProjectTorreyPines registry):

```julia
using Pkg
Pkg.add("TurbulentTransport")   # or: Pkg.develop("TurbulentTransport") for a dev checkout
```

Run a neural-network TGLF surrogate from an IMAS ODS:

```julia
using TurbulentTransport, IMAS

dd = IMAS.json2imas("path/to/ods.json")
rho = [0.3, 0.5, 0.7]

# Build TGLF inputs at the requested radii (sat3, electromagnetic)
input_tglfs = TurbulentTransport.InputTGLF(dd, rho, :sat3, true, false).tglfs

# Evaluate a TGLF-NN model (returns one GACODE.FluxSolution per radius)
fluxes = run_tglfnn(input_tglfs; model_filename="sat3_em_d3d_azf-1", warn_nn_train_bounds=true)
@show fluxes[1].ENERGY_FLUX_e fluxes[1].ENERGY_FLUX_i
```

Not sure which model to use? Start with [`model_selector`](#picking-a-model-with-model_selector)
to discover the best-fitting model for your case before committing to a
`model_filename`.

> **Model weights:** weights are stored under `models/` as Git LFS files. A plain
> `Pkg` install does not run `git lfs pull`, so the package transparently
> downloads (and SHA-verifies) each weight file the first time it is needed — no
> manual step required as long as you have network access.

## Examples and tutorials

In-repo notebooks (the place to start):

- [`examples/run_TGLFNN.ipynb`](https://github.com/ProjectTorreyPines/TurbulentTransport.jl/blob/master/examples/run_TGLFNN.ipynb)
  — **start here**: write/load an `input.tglf`, run TGLF / TJLF / TGLF-NN, use
  `model_selector`, and a ForwardDiff sensitivity demo.
- [`examples/uncertainty_input_output.ipynb`](https://github.com/ProjectTorreyPines/TurbulentTransport.jl/blob/master/examples/uncertainty_input_output.ipynb)
  — propagate input/output uncertainty using NN ensembles and `Measurements.jl`.
- [`examples/spot_check.ipynb`](https://github.com/ProjectTorreyPines/TurbulentTransport.jl/blob/master/examples/spot_check.ipynb)
  — sanity-check a model's flux outputs against expected values.

Utility notebook:

- [`utilities/convert_nn.ipynb`](https://github.com/ProjectTorreyPines/TurbulentTransport.jl/blob/master/utilities/convert_nn.ipynb)
  — convert trained NN weights into the `.bson` format this package loads.

FUSE worked examples (the [FuseExamples](https://github.com/ProjectTorreyPines/FuseExamples)
repository), which exercise these models end-to-end:

- [`fluxmatcher.ipynb`](https://github.com/ProjectTorreyPines/FuseExamples/blob/master/fluxmatcher.ipynb)
  — the canonical flux-matcher tutorial driving TGLF-NN / QLNN.
- [`study_TGLFdb.ipynb`](https://github.com/ProjectTorreyPines/FuseExamples/blob/master/study_TGLFdb.ipynb)
  and [`study_database_generator.ipynb`](https://github.com/ProjectTorreyPines/FuseExamples/blob/master/study_database_generator.ipynb)
  — TGLF database generation and analysis.
- [`tutorial.ipynb`](https://github.com/ProjectTorreyPines/FuseExamples/blob/master/tutorial.ipynb)
  — general FUSE entry point for context.

## Model families explained

All NN surrogates take the same TGLF inputs as the first-principles codes, so you
can swap between them without rebuilding your problem.

- **TGLF-NN** (`run_tglfnn`, `run_tglfnn_onnx`) — a drop-in surrogate for TGLF
  quasi-linear fluxes (electron/ion energy, particle, momentum) evaluated from an
  `InputTGLF`. Roughly milliseconds per call versus seconds for TGLF. Ensemble models
  return uncertainty with `uncertain=true` (via `Measurements.jl`), and inference is
  ForwardDiff-compatible so it works inside the FUSE flux matcher. `run_tglfnn_onnx`
  is the ONNXRuntime variant.

- **GKNN** (`run_tglfnn(...; fidelity=:GKNN)`) — a multi-fidelity *correction* layered
  on a TGLF-NN base model that nudges the fluxes toward gyrokinetic (CGYRO) truth.
  Includes core / near-edge / edge-blended DIII-D variants. Reach for it when you want
  higher fidelity than TGLF alone.

- **QLNN** (`run_qlnn`, `loadqlnnbundle`) — predicts per-`ky` **quasi-linear weights and
  eigenvalues** (energy / particle / momentum / eigenvalue regressors, plus an optional
  stability classifier) and feeds them through TJLF's saturation rule, so `SAT_RULE`,
  `ALPHA_ZF`, and `UNITS` from the input still apply. This is more physics-faithful than
  flux-level TGLF-NN and can additionally reconstruct 2D fluctuation spectra
  (`qlnn_fluctuation_spectra`, with a GPU path via the `CUDA` extension). QLNN bundles
  live in `models/QLNN*` directories.

- **FINN** (`run_finn`) — a flux-matcher *inversion* network. Instead of iterating a
  transport solve, it predicts the converged flux-matched gradients (a/L_Te, a/L_Ti,
  a/L_ne, and ExB shear) directly from geometry and sources in a single forward pass.

- **ModeID** (`run_modeid_nn`, `run_modeid_qlnn`) — classifies the **dominant turbulence
  mode** (ITG / TEM / KBM / ETG / MTM) at each radius. `run_modeid_nn` is a direct softmax
  classifier from TGLF inputs; `run_modeid_qlnn` derives the mode from QLNN QL-weight
  ratios and eigenvalue frequency, mirroring TJLF's `identify_modes`.

The first-principles drivers `run_tglf` (Fortran TGLF), `run_tjlf` (Julia TJLF), and
`run_qlgyro` (QLGYRO) share the same input types and `FluxSolution` outputs, which makes
them convenient ground truth for the surrogates.

### Model naming and discovery

Trained TGLF-NN / GKNN models follow the convention

```
<sat_rule>_<es|em>_<device>_azf<±1>[_variant]
```

for example `sat3_em_d3d_azf-1` (SAT3, electromagnetic, DIII-D, `ALPHA_ZF=-1`). Devices
include `d3d`, `mastu`, `nstx`, `ukstep`, `iter`, `fpp`/`stfpp`, multi-machine combinations
(`d3d+mastu+nstx`), edge variants (`d3dedge`, `d3dnearedge`), and `withnegD` (negative
triangularity). Suffixes such as `gknn*`, `tglfnn24` denote correction/companion models.

List what is installed with `available_models()` (TGLF-NN / GKNN model files) and
`available_qlnn_bundles()` (QLNN bundle directories).

### Picking a model with `model_selector`

With many models available, `model_selector` is the recommended way to choose one. Given
an ODS path (or an `IMAS.dd`) and a `rho_grid`, it runs the available TGLF-NN models, then
ranks them per radius by confidence (ensemble uncertainty) and — when `ground_truth=true`
(default) — by accuracy against TJLF / Fortran TGLF, returning the top `max_models`:

```julia
results = model_selector("path/to/ods.json"; rho_grid=range(0.1, 0.9, 9))

# best models at the first radius
results.rankings[1].top_models
```

Useful keywords: `sat_rule`, `electromagnetic`, `filter_sat_rule`, `max_models`,
`ground_truth`. See
[`examples/run_TGLFNN.ipynb`](https://github.com/ProjectTorreyPines/TurbulentTransport.jl/blob/master/examples/run_TGLFNN.ipynb)
for a full walkthrough.

## Utilities and CLI

Helpers in `src/utils.jl`:

- `load` / `save` — read and write `input.tglf`, `input.cgyro`, and `input.qlgyro` files.
- `scan` — build a vector of `InputTGLF` over a parameter scan.
- `compare_two_input_tglfs` — diff two inputs field-by-field.
- `parse_out_tglf_gbflux` — parse a TGLF `out.tglf.gbflux` file.

GACODE-style command-line drivers in `bin/` run a simulation directory through the
Julia implementations:

```bash
tjlf   -e .     # run TJLF on input.tglf/input.tjlf in the current directory
qlgyro -e .     # run QLGYRO
```

## Use in FUSE, TJLF, and TJLFEP

- **FUSE** drives these models through its transport actors: `ActorTGLF` calls
  `run_tglfnn` / `run_tglfnn_onnx` / `run_qlnn`, the flux matcher uses the
  ForwardDiff-compatible `run_tglfnn` / `run_qlnn` paths, and `ActorFINN`, `ActorModeID`,
  `ActorQLGYRO`, and `ActorTGLFEP` map to `run_finn`, `run_modeid_nn` / `run_modeid_qlnn`,
  `run_qlgyro`, and the TGLF-EP scan driver.

- **TJLF** provides the saturation-rule physics and the `InputTGLF` / `InputTJLF` types
  that this package builds and reuses, including for the QLNN saturation-rule integration.

- **TJLFEP** (the energetic-particle stability code) builds its TGLF-EP inputs from
  `TurbulentTransport.InputTGLFEP` and can be launched on SLURM via `run_tjlfep`.

## Citation
If this software contributes to an academic publication, please cite it as follows:

TGLF-NN:
> Tom Neiser, Orso Meneghini, Sterling Smith, Joseph McClenaghan, David Orozco, Joseph Hall, Gary Staebler, Emily Belli, Jeff Candy, _Database generation for validation of TGLF and retraining of neural network accelerated TGLF-NN_, APS Division of Plasma Physics Meeting Abstracts (2022)

GKNN:
> Tom Neiser, Orso Meneghini, Sterling Smith, Joseph McClenaghan, Tim Slendebroek, David Orozco, Brian Sammuli, Gary Staebler, Joseph Hall, Emily Belli, Jeff Candy, _Multi-fidelity neural network representation of gyrokinetic turbulence_, APS Division of Plasma Physics Meeting Abstracts (2023)
