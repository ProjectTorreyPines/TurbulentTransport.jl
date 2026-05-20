# ONNX to PyTorch and using models in Fortran

The script `convert_nn_onnx_to_pytorch_trace.py` takes multiple ONNX model files and combines them into a single PyTorch traced model, saving the model to `committee.pt`. The combined model is scaled according to `xm.txt`, `xsigma.txt`, `ym.txt` and `ysigma.txt`, so that the inputs and outputs correspond to those used in TGLF. The model file `committee.pt` can then be used in PyTorch, or in Fortran via [FTorch](https://github.com/Cambridge-ICCS/FTorch).

The traced model requires 1 input array that contains the inputs as in `xnames.txt` (model dependent length) and outputs 2 arrays (prediction and variance) of length 4, corresponding to the names in `ynames.txt`.

## Running the Script

The `convert_nn_onnx_to_pytorch_trace.py` script combines all `.onnx` files in the local directory into `committee.pt`. Given the prerequisites below, you can run the script from a model directory (e.g. from `models/sat3_em_nstx_azf-1/`):
```bash
python ../../utilities/convert_nn_onnx_to_pytorch_trace.py
```

For residual models (e.g. `models/sat3_em_d3d+mastu+nstx_azf-1_gknn31/`) it is necessary to first run the script in the parent model directory since the parent `committee.pt` is required in conversion.

### Prerequisites

Ensure the following dependencies are installed in your Python environment:

- Python 3.8+
- `onnx2pytorch` (for converting ONNX to PyTorch)
- `torch` and `torchvision` (PyTorch for model manipulation and tracing)
- `numpy` (for handling the input scaling file)

Note that one can also use the python installation of `torch` for FTorch in Fortran. For use on CPU, it is suggested to first install
compatible `torch` and `torchvision` with
```bash
pip install torch torchvision --index-url https://download.pytorch.org/whl/cpu

```
then install `onnx2pytorch` afterwards.

### Notes

- The script assumes all `.onnx` files have compatible input/output structures.
- The script requires `xm.txt`, `xsigma.txt`, `ym.txt` and `ysigma.txt` to be in the model directory.
- If there are issues with `onnx2pytorch`, it might be replaced with [`onnx2torch`](https://github.com/ENOT-AutoDL/onnx2torch)
 and its method `convert()` in place of `ConvertModel()`.

## Using the PyTorch traced model in Fortran

The traced PyTorch model (`committee.pt`) can be used in Fortran with [FTorch](https://github.com/torch/FTorch), which provides Fortran bindings for LibTorch (the C++ backend of PyTorch).

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
  ! models require 1 array of inputs per evaluation (n_in = 1)
  ! models output both mean and variance (n_out = 2)
  integer, parameter :: n_in = 1, n_out = 2
  type(torch_tensor), dimension(n_in),  save :: input_tensor
  type(torch_tensor), dimension(n_out), save :: output_tensor

  ! Arrays to use in Fortran (input_array, output_mean, output_var)
  ! Note that we need to use 2D arrays that are of shape (1,:) to conform with
  ! the traced model.
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
      allocate(input_array(1,input_size))
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

      ! Copy inputs
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

  allocate(model_input(in_size))
  allocate(model_output(4))      ! 4 Fluxes returned

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
