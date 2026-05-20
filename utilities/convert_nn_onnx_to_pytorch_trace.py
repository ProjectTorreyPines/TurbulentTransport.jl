#!/usr/bin/env python
# coding: utf-8

# Takes all .onnx files in working directory and converts them to a single
# .pt file of comittee PyTorch trace. Assumes 'xm.txt' is present and can
# be used as example model input (to generate trace).
# 
# Working model directories are treated as either a 'regular' model where
# only the local files are required or a 'residual' model, where a 'regular'
# parent model committee trace is assumed to already exist (run this script
# in the parent directory first in that case).

import onnx
import torch
import torch.nn as nn
import numpy as np
from onnx2pytorch import ConvertModel
import os
import re

dir_name = os.path.basename(os.getcwd())

print(f"Processing: {dir_name}")

# ====================== RESIDUAL DETECTION ======================
is_residual = False
parent_committee = None
if '_' in dir_name and not re.search(r'-\d+$', dir_name): # possible parent suffix without "-number"
    potential_parent_name = dir_name.rsplit('_', 1)[0]
    potential_parent_path = os.path.join('..', potential_parent_name)

    if os.path.exists(potential_parent_path):
        try:
            parent_xnames = open(os.path.join(potential_parent_path, 'xnames.txt')).read().strip().split()
            parent_ynames = open(os.path.join(potential_parent_path, 'ynames.txt')).read().strip().split()
            current_xnames = open('xnames.txt').read().strip().split()
            # Input and output names of parent model match inputs of local directory.
            if parent_xnames + parent_ynames == current_xnames:
                is_residual = True
                parent_path = potential_parent_path

                committee_file = os.path.join(parent_path, 'committee.pt')
                if os.path.exists(committee_file):
                    parent_committee = torch.jit.load(committee_file)
                    print(f"✓ Residual model detected. Parent: {potential_parent_name}")
                else:
                    print(f"⚠ Residual model detected (parent: {potential_parent_name}), but committee.pt not found.")
                    print("   Please run this script in the parent directory first.")
                    is_residual = False  # treat as regular model
            else:
                print("✗ xnames check failed → treating as regular model")
        except Exception as e:
            print(f"✗ Could not read xnames files: {e} → treating as regular model")
    else:
        print("✗ Parent directory not found → treating as regular model")

if not is_residual:
    print("✓ Treating as regular model")

# === LOAD SCALINGS ===
# Load input mean, sd values
xm       = np.loadtxt('xm.txt').astype(np.float32)
model_xm = torch.from_numpy(xm)
model_xs = torch.from_numpy(np.loadtxt('xsigma.txt').astype(np.float32))

# Load scalings required to get result in standard TGLF units.
# N.B. should not be used in residual models.
if not is_residual:
    model_ym = torch.from_numpy(np.loadtxt('ym.txt').astype(np.float32))
    model_ys = torch.from_numpy(np.loadtxt('ysigma.txt').astype(np.float32))

# === PUBLIC INTERFACE ===
# Use xm as example input to construct the public interface.
# Residual models should have the same number of inputs as the parent model.
if is_residual and parent_committee is not None:
    # Number of parent outputs = length of parent_ynames
    num_parent_outputs = len(parent_ynames)
    # Take only the original inputs (strip the last num_parent_outputs values)
    public_xm = xm[0 : len(xm) - num_parent_outputs]
    example_input = torch.as_tensor(public_xm.reshape(1, len(public_xm)), dtype=torch.float32)
else:
    example_input = torch.as_tensor(xm.reshape(1, len(xm)), dtype=torch.float32)

# === LOAD ONNX MODELS AND CONVERT TO PT ===
# Get list of onnx model files.
onnx_model_files = [f for f in os.listdir('.') if f.endswith('.onnx')]

pytorch_models = []  # List to store the converted PyTorch models
for model_file_in in onnx_model_files:
    print(f" {model_file_in}")
    onnx_model = onnx.load(model_file_in)
    pytorch_models.append(ConvertModel(onnx_model)) # Convert and add to the list

# === CONSTRUCT MODEL ===

# For residual models, call forward method of parent model and concatenate its
# output with original inputs.
if is_residual:
    class ResidualCommittee(nn.Module):
        def __init__(self, models, parent):
            super().__init__()
            self.models = nn.ModuleList(models)
            self.parent = parent

        def forward(self, input):
            parent_mean, parent_var = self.parent(input)

            local_input    = torch.cat([input, parent_mean], dim=1)
            scaled_input   = (local_input - model_xm) / model_xs

            outputs_list = []
            i = 0
            while i < len(self.models):
                m = self.models[i]
                output = m(scaled_input)
                outputs_list.append(output)
                i = i + 1

            all_outputs = torch.stack(outputs_list)

            error_factor_mean = torch.mean(all_outputs, dim=0)
            error_factor_var  = torch.var(all_outputs, dim=0)

            return error_factor_mean * parent_mean, error_factor_mean ** 2 * parent_var
    committee = ResidualCommittee(pytorch_models, parent_committee)

else:
    class Committee(nn.Module):
        def __init__(self, models):
            super().__init__()
            self.models = nn.ModuleList(models)

        def forward(self, input):
            scaled_input = (input - model_xm) / model_xs

            outputs_list = []
            i = 0
            while i < len(self.models):
                m = self.models[i]
                output = m(scaled_input)
                outputs_list.append(output)
                i = i + 1

            all_outputs = torch.stack(outputs_list)

            mean = torch.mean(all_outputs, dim=0)
            var = torch.var(all_outputs, dim=0)

            return model_ys * mean + model_ym, model_ys ** 2 * var

    committee = Committee(pytorch_models)

# === Trace and save ===

# Trace the committee
traced_committee = torch.jit.trace(committee, example_input)

# Save the TorchScript committee
output_filename = "committee.pt"
traced_committee.save(output_filename)
print(f"Committee TorchScript saved to {output_filename}")

if is_residual and parent_committee is not None:
    print("This is a residual model")
else:
    print("This is a regular model")
