# Test fixtures and helper functions for TurbulentTransport tests

# Path to sample input files
const TEST_DATA_DIR = joinpath(@__DIR__, "data")
const SAMPLE_INPUT_PATH = joinpath(TEST_DATA_DIR, "sample_input.tglf")

# Known good model filenames for testing
const TEST_MODEL_SINGLE = "sat3_em_d3d_azf-1"
const TEST_MODEL_ENSEMBLE = "sat3_em_d3d_azf-1"  # This is an ensemble model
const TEST_MODEL_GKNN = "sat3_em_d3d_azf-1_gknne24"  # GKNN correction model

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

# Expected baseline values for regression testing
# Captured from v1.0.14 to detect any numerical changes

# flux_array expected outputs for sat0_em_d3d single model (first model in ensemble)
const EXPECTED_FLUX_ARRAY_SINGLE = [0.021994637644215054, 0.14526707591951704, 0.4718316916684522, 0.4132556522695987]

# flux_array expected outputs for sat0_em_d3d single model with fidelity=:GKNN
# Note: This is a TGLF-NN model, so fidelity=:GKNN just skips denormalization (raw NN output)
const EXPECTED_FLUX_ARRAY_SINGLE_RAW = [-0.23004951851133254, -0.21326135534711577, -0.31156058557178723, -0.32071940967901724]

# GKNN correction model (sat3_em_d3d_azf-1_gknne24) expected outputs
# These models have ynames of length 2, and fidelity=:GKNN outputs div(ynames, 2) = 1 value
const EXPECTED_GKNN_MODEL_SINGLE = [0.93651175002966]
const EXPECTED_GKNN_MODEL_ENSEMBLE = [0.9496942106346887]

# flux_array expected outputs for sat3_em_d3d_azf-1 ensemble
const EXPECTED_FLUX_ARRAY_ENSEMBLE = [-0.019443083518266357, 0.04080973584865324, 0.23715745308655625, 0.07256906359102389]

# run_tglfnn expected outputs with sample_input.tglf
const EXPECTED_RUN_TGLFNN_SAT3 = (
    ENERGY_FLUX_e = 3.612103806515151,
    ENERGY_FLUX_i = 6.160152661030724,
    PARTICLE_FLUX_e = 0.6037634688605953,
    STRESS_TOR_i = 2.5092916159367995,
)

const EXPECTED_RUN_TGLFNN_SAT2 = (
    ENERGY_FLUX_e = 3.3986959391828035,
    ENERGY_FLUX_i = 6.599432537271315,
    PARTICLE_FLUX_e = 0.6049814789721364,
    STRESS_TOR_i = 3.4500222580700415,
)

const EXPECTED_RUN_TGLFNN_SAT3_GKNN = (
    ENERGY_FLUX_e = 2.3269702143661934,
    ENERGY_FLUX_i = 3.8352993604857386,
    PARTICLE_FLUX_e = 0.5105086253045001,
    STRESS_TOR_i = 2.5997926311395063,
)
