# ONNX to PyTorch and using models in Fortran

The script `convert_nn_onnx_to_pytorch_trace.py` takes multiple ONNX model files and combines them into a single PyTorch traced model, saving the model to `committee.pt`. The combined model is scaled according to `xm.txt`, `xsigma.txt`, `ym.txt` and `ysigma.txt`, so that the inputs and outputs correspond to those used in TGLF. The model file `committee.pt` can then be used in PyTorch, or in Fortran via [FTorch](https://github.com/Cambridge-ICCS/FTorch).

The traced model takes one input array with the features in `xnames.txt` (model-dependent length) and returns two arrays (mean and variance) of length 4, corresponding to the names in `ynames.txt`. The first dimension is batch: a single point is shape `(1, n_in)`, a profile is `(n, n_in)`. See [Batching](#batching).

## Running the Script

The `convert_nn_onnx_to_pytorch_trace.py` script combines all `.onnx` files in the local directory into `committee.pt`. Given the prerequisites below, you can run the script from a model directory (e.g. from `models/sat3_em_nstx_azf-1/`):
```bash
python ../../utilities/convert_nn_onnx_to_pytorch_trace.py
```

By default this also checks the conversion against ONNX Runtime (see [Verification](#verification)). To skip that check (and the optional `onnxruntime` dependency):
```bash
python ../../utilities/convert_nn_onnx_to_pytorch_trace.py --no-verify
```

For residual models (e.g. `models/sat3_em_d3d+mastu+nstx_azf-1_gknn31/`) it is necessary to first run the script in the parent model directory since the parent `committee.pt` is required in conversion.

### Prerequisites

Ensure the following dependencies are installed in your Python environment:

- Python 3.8+
- `onnx` (to load the member graphs)
- `onnx2pytorch` (for converting ONNX to PyTorch)
- `torch` and `torchvision` (PyTorch for model manipulation and tracing)
- `numpy` (for handling the input scaling file)

Optional, used only for the post-trace check:

- `onnxruntime` (ONNX Runtime). Not required to write `committee.pt`. If it is not installed, pass `--no-verify`.

Note that one can also use the python installation of `torch` for FTorch in Fortran. For use on CPU, it is suggested to first install
compatible `torch` and `torchvision` with
```bash
pip install torch torchvision --index-url https://download.pytorch.org/whl/cpu

```
then install `onnx2pytorch` afterwards. Install `onnxruntime` the same way if you want the verification.

### Verification

The `.onnx` files are the export snapshot. The script does not go back to Julia/Flux; it only needs to show that `ConvertModel`, the committee wrap, and `torch.jit.trace` match that snapshot.

Unless `--no-verify` is passed, it:

1. Runs each member through ONNX Runtime and through the converted PyTorch module, on `xm` plus a few Gaussian samples around it (`--n-random`, default 8).
2. Composes those ONNX outputs into committee mean and variance (and, for residual models, the parent members the same way).
3. Checks the eager committee against that composition.
4. Traces, then checks the traced module against the same ONNX composition and against eager.

Tolerances default to `1e-4` (`--atol`, `--rtol`). A mismatch exits before `committee.pt` is written.

`onnxruntime` is an optional dependency: it is imported only for these checks. Skip them with `--no-verify` if it is not installed or if you only want to produce the trace.

### Notes

- The script assumes all `.onnx` files have compatible input/output structures.
- The script requires `xm.txt`, `xsigma.txt`, `ym.txt` and `ysigma.txt` to be in the model directory.
- If there are issues with `onnx2pytorch`, it might be replaced with [`onnx2torch`](https://github.com/ENOT-AutoDL/onnx2torch)
 and its method `convert()` in place of `ConvertModel()`.

## Batching

The Fortran example below is a **single-batch** call: tensors of shape `(1, n_in)` and `(1, 4)`, one physical point per `run_model`. That matches a scalar evaluate.

The traces are not limited to `n = 1`. They take a leading batch dimension `(n, n_in)` and return mean and variance each `(n, 4)`. A single row pays a large dispatcher cost (about 1 ms for a regular 20-member committee). Independent points — a radial profile, a grid — are much cheaper per point if they go in one call: 16 rows of the same net are only about 1.5 ms total, so looping 16 scalar evaluates is roughly ten times more expensive. Residual traces behave the same way at about 2.4× those times.

If the host already holds a vector of states, pass them as one `(n, n_in)` array rather than looping the single-row example.

## Using the PyTorch traced model in Fortran

The traced PyTorch model (`committee.pt`) can be used in Fortran with [FTorch](https://github.com/torch/FTorch), which provides Fortran bindings for LibTorch (the C++ backend of PyTorch). The listing below is the single-batch case (`n = 1`); see [Batching](#batching) if you have more than one point.

### Prerequisites

- **LibTorch**: Download the appropriate version (CPU or GPU) from the [PyTorch website](https://pytorch.org/get-started/locally/) and ensure it is accessible in your environment. CPU versions of the `LibTorch` and `Pip` packages have been tested. The `LibTorch` version requires no Python to install or run. It is suggested to look at the `FTorch` instructions below first.
- **FTorch**: Install the FTorch library following the instructions in the [FTorch repository](https://github.com/torch/FTorch). This also provides a compiler specific module (`ftorch.mod`).
- **Fortran Compiler**: Use a modern Fortran compiler (e.g., `gfortran` or `ifort`) compatible with FTorch.
- **CMake**: Check for version required. Required to build FTorch. Not essential, but helpful for building final Fortran code.

### Example Fortran Code

Below is an example Fortran *module* that deals with loading and running the model via FTorch and an example *program* that utilizes it.

The example program is specific to `sat3_em_nstx_azf-1`, because it inputs 31 values according to the `xnames.txt` file found in that model directory. The module is not model specific, but assumes 4 outputs, as is presently the case for all models.
The input size should be specified according to the model.

The code assumes `committee.pt` is present in local the run directory, but can be adapted as required (might be given as input, for example).

The module does not return the model variance, but can be easily adapted to do so.

```fortran
module tglf_torchscript

  use, intrinsic :: iso_fortran_env, only : wp => real32
  use ftorch

  implicit none

  ! Models
  type(torch_model), save :: model

  ! Tensors to use in torch
  ! n_in / n_out are FTorch tensor counts, not batch size:
  ! one input array per call, two outputs (mean and variance).
  integer, parameter :: n_in = 1, n_out = 2
  type(torch_tensor), dimension(n_in),  save :: input_tensor
  type(torch_tensor), dimension(n_out), save :: output_tensor

  ! Arrays to use in Fortran (input_array, output_mean, output_var)
  ! First index is batch. This example uses n = 1, so shape (1,:).
  ! The trace accepts (n, n_in); a profile of points in one call is
  ! much cheaper per point than looping this scalar path.
  integer, parameter :: in_dims = 2     ! dimensions
  integer :: layout(in_dims) = [1,2] ! array dimension ordering
  real (wp), allocatable, dimension(:,:), target :: input_array
  real (wp), allocatable, dimension(:,:), target :: output_mean, output_var

  private
  public :: setup_model, run_model

  contains

    ! call setup once before running model
    subroutine setup_model(input_size)

      ! Provide number of model inputs and outputs
      ! e.g. 31 inputs and 4 outputs for sat3_em_nstx_azf-1
      integer, intent(in) :: input_size
      integer, parameter :: output_size = 4

      ! Load model, assuming "committee.pt" is in the local directory;
      ! modify accordingly.
      call torch_model_load(model, "committee.pt", torch_kCPU)

      ! Allocate input and output arrays
      ! First dim = 1 is a single point. For a profile, use (n, input_size)
      ! / (n, output_size) and rebuild the tensors once for that n.
      allocate(input_array(1,input_size)) ! single batch
      allocate(output_mean(1,output_size), output_var(1,output_size))

      ! Set up pointers to inputs and outputs
      call torch_tensor_from_array(input_tensor(1),  input_array, layout, torch_kCPU)
      call torch_tensor_from_array(output_tensor(1), output_mean, layout, torch_kCPU)
      call torch_tensor_from_array(output_tensor(2), output_var,  layout, torch_kCPU)

    end subroutine setup_model

    subroutine run_model(input,output) ! add variance if required

      ! These inputs and outputs need not be the same type of real as
      ! associated with the torch calls below.
      real, dimension(:), intent(in)  :: input
      real, dimension(:), intent(out) :: output

      ! Copy one point into row 1. Looping this is the slow path;
      ! pass a rank-2 (n, n_in) array instead if you have many points.
      input_array(1,:) = input

      ! Call traced model
      call torch_model_forward(model, input_tensor, output_tensor)

      ! Map returned outputs (add variance if required)
      output(1:4) = output_mean(1,1:4)

    end subroutine run_model

end module tglf_torchscript
```

```
program call_tglf_torchscript

  ! optional: possibly re-use structures from TGLF?
  !use tglf_interface

  use tglf_torchscript, only : setup_model, run_model
  implicit none

  ! sat3_em_nstx_azf-1 specific
  integer, parameter :: in_size = 31
  real, allocatable, dimension(:) :: model_input  ! Need not be single
  real, allocatable, dimension(:) :: model_output ! precision reals.

  allocate(model_input(in_size))  ! single point
  allocate(model_output(4))       ! 4 Fluxes returned

  call setup_model(in_size)

  ! Construct inputs according to model "xnames.txt".
  ! (here the contents of "xm.txt" have been used)
  model_input(:) =[ 0.6406908,   & ! AS_2, 
                    0.057125505, & ! AS_3
                   -2.0141907,   & ! BETAE_log10
                   -2.071454,    & ! DEBYE_log10
                    0.22197191,  & ! DELTA_LOC
                   -0.33162355,  & ! DRMAJDX_LOC
                   -0.031959284, & ! DZMAJDX_LOC
                    2.2358067,   & ! KAPPA_LOC
                   -0.017725442, & ! P_PRIME_LOC
                    3.7026527,   & ! Q_LOC
                   74.82468,     & ! Q_PRIME_LOC
                    0.82765865,  & ! RLNS_1
                    1.0271144,   & ! RLNS_2
                    0.39525518,  & ! RLNS_3
                    1.957151,    & ! RLTS_1
                    2.1438692,   & ! RLTS_2
                    2.1776273,   & ! RLTS_3
                    1.619759,    & ! RMAJ_LOC
                    0.5672825,   & ! RMIN_LOC
                    0.28694844,  & ! S_DELTA_LOC
                    0.0745877,   & ! S_KAPPA_LOC
                    0.018578079, & ! S_ZETA_LOC
                    1.0358093,   & ! TAUS_2
                    1.0382642,   & ! TAUS_3
                   -0.08510396,  & ! VEXB_SHEAR
                    0.37218216,  & ! VPAR_1
                    0.7683608,   & ! VPAR_SHEAR_1
                   -0.23693003,  & ! XNUE_log10
                    2.7133608,   & ! ZEFF
                   -0.014444544, & ! ZETA_LOC
                   -0.032985397 ]  ! ZMAJ_LOC

  ! run the forward() method and return means in model_output
  call run_model(model_input,model_output)

  ! Outputs, as in "ynames.txt" (TGLF units)
  write(*,*) 'OUT_G_elec: ', model_output(1)
  write(*,*) 'OUT_P_ions: ', model_output(2)
  write(*,*) 'OUT_Q_elec: ', model_output(3)
  write(*,*) 'OUT_Q_ions: ', model_output(4)

end program call_tglf_torchscript
```
### Notes

- **TGLF Input Mapping**: Details of the expected inputs in `xnames.txt` and their mappings to TGLF can be found in the [TGLF documentation](https://gafusion.github.io/doc/tglf/tglf_table.html).

For questions or contributions, please open an issue or pull request on the [TurbulentTransport.jl GitHub repository](https://github.com/ProjectTorreyPines/TurbulentTransport.jl).
