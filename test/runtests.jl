using Test
using TurbulentTransport
using TurbulentTransport: InputTGLF, flux_array, flux_solution, loadmodel, available_models
using TurbulentTransport: TGLFNNmodel, TGLFNNensemble
import Flux
using TurbulentTransport.Measurements

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
    include("test_flux_array.jl")
    include("test_run_tglfnn.jl")
    include("test_utils.jl")
    include("test_input_tglfs.jl")
    include("test_pooled_chain.jl")
    include("test_model_selector.jl")
end
