# Test fixtures and helper functions for TurbulentTransport tests

# Path to sample input files
const TEST_DATA_DIR = joinpath(@__DIR__, "data")
const SAMPLE_INPUT_PATH = joinpath(TEST_DATA_DIR, "sample_input.tglf")

# Known good model filenames for testing
const TEST_MODEL_SINGLE = "sat3_em_d3d_azf-1"
const TEST_MODEL_ENSEMBLE = "sat3_em_d3d_azf-1"  # This is an ensemble model

# Expected field values after loading sample_input.tglf
const EXPECTED_LOAD_VALUES = (
    NS = 3,
    SAT_RULE = 3,
    BETAE = 0.00362972,
    Q_LOC = 2.00545,
    KAPPA_LOC = 1.40438,
    DELTA_LOC = 0.0681444,
    RMAJ_LOC = 2.86212,
    RMIN_LOC = 0.573129,
    MASS_1 = 0.000272445,
    MASS_2 = 1.0,
    MASS_3 = 6.0,
    ZS_1 = -1.0,
    ZS_2 = 1.0,
    ZS_3 = 6.0,
    AS_1 = 1.0,
    AS_2 = 0.784867,
    AS_3 = 0.0302081,
    USE_BPER = true,
    USE_BPAR = true,
    UNITS = "GYRO",
)

# Load sample InputTGLF
function load_sample_input()
    TurbulentTransport.load(InputTGLF(), SAMPLE_INPUT_PATH)
end

# Helper function to generate valid test input for a model
# The model's xbounds are in the TRANSFORMED space (after log10 for _log10 fields)
# but flux_array expects ORIGINAL values (before log10)
function generate_valid_input(model)
    n_inputs = length(model.xnames)
    x = zeros(Float64, n_inputs)

    for (i, name) in enumerate(model.xnames)
        # Get midpoint in the transformed (training) space
        mid_transformed = (model.xbounds[i, 1] + model.xbounds[i, 2]) / 2

        if contains(name, "_log10")
            # xbounds are in log10 space, convert back to original
            # e.g., bounds [-4, -1] in log10 space → midpoint -2.5 → original value 10^(-2.5)
            x[i] = 10.0^mid_transformed
        else
            x[i] = mid_transformed
        end
    end

    return x
end

# Helper to generate matrix input
function generate_valid_input_matrix(model, n_samples::Int)
    n_inputs = length(model.xnames)
    x = zeros(Float64, n_inputs, n_samples)

    base_input = generate_valid_input(model)
    for j in 1:n_samples
        x[:, j] = base_input
    end

    return x
end
