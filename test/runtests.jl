using Test
using TurbulentTransport
using TurbulentTransport: InputTGLF, InputTJLF, flux_array, flux_solution, loadmodel, available_models
using TurbulentTransport: TGLFNNmodel, TGLFNNensemble
import Flux
using TurbulentTransport.Measurements
using TurbulentTransport.GACODE

# Relative tolerance for regression tests (allows cross-platform floating-point variance)
const REGRESSION_RTOL = 1e-12

# Include test fixtures
include("fixtures.jl")

# Check if specific test files are requested via ARGS
if !isempty(ARGS)
    for testfile in ARGS
        @info "Running test file: $testfile"
        include(testfile)
    end
else
    include("test_load.jl")
    include("test_presets.jl")
    include("test_models.jl")

    include("test_readonly_cache.jl")
    include("test_flux_array.jl")
    include("test_run_tglfnn.jl")
    include("test_finn.jl")
    include("test_modeID.jl")
    include("test_utils.jl")
    include("test_input_tglfs.jl")
    include("test_pooled_chain.jl")
    include("test_dual_pool_compatibility.jl")
    include("test_model_selector.jl")
    include("test_qlgyro.jl")
    include("test_qlnn.jl")
    include("test_qlnn_fluctuation_spectra.jl")
    include("test_tglf.jl")
    include("test_tjlf_modes.jl")
    include("test_slurm_utils.jl")
    include("test_tglf_ep.jl")
    include("test_dd_inputs.jl")
end
