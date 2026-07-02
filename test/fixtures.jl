# Test fixtures and helper functions for TurbulentTransport tests

# Path to sample input files
const TEST_DATA_DIR = joinpath(@__DIR__, "data")
const SAMPLE_INPUT_PATH = joinpath(TEST_DATA_DIR, "sample_input.tglf")

# A small, self-contained IMAS `dd` with both equilibrium and core_profiles,
# used to exercise the `dd`-based input constructors (InputTGLF/InputCGYRO/
# InputTGLFEP) without pulling in FUSE. Copied verbatim from IMASdd's
# `sample/omas_sample.json`.
const SAMPLE_DD_PATH = joinpath(TEST_DATA_DIR, "sample_dd.json")

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

# Load the sample IMAS `dd` and select a valid global time. IMAS is reached via
# the TurbulentTransport namespace so no extra test-project dependency is needed.
function load_sample_dd()
    dd = TurbulentTransport.IMAS.IMASdd.json2imas(SAMPLE_DD_PATH; show_warnings=false)
    if !isempty(dd.equilibrium.time)
        dd.global_time = dd.equilibrium.time[end]
    end
    return dd
end

# Run `f()` with a throwaway `sacct` shim on PATH so the SLURM-polling helpers
# can be exercised without a real scheduler. The emitted job State is `state`;
# pass "__FAIL__" to make the shim exit non-zero (exercising the catch path).
function with_fake_sacct(f, state::AbstractString)
    mktempdir() do bin
        shim = joinpath(bin, "sacct")
        open(shim, "w") do io
            println(io, "#!/usr/bin/env bash")
            println(io, "if [ \"\$FAKE_SACCT_STATE\" = \"__FAIL__\" ]; then exit 1; fi")
            println(io, "printf '%s\\n' \"\$FAKE_SACCT_STATE\"")
        end
        chmod(shim, 0o755)
        withenv("PATH" => bin * ":" * get(ENV, "PATH", ""), "FAKE_SACCT_STATE" => state) do
            f()
        end
    end
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

# ============================================
# Regression Test Expected Values
# Captured on 2025-12-06 for single/vector equivalence
# ============================================

"""
Create test input variations for regression testing.
Returns 3 inputs: base, modified Q_LOC/BETAE, modified RLTS_2/VEXB_SHEAR.
"""
function create_regression_inputs()
    input1 = load_sample_input()

    # Variation 2: Modified Q_LOC and BETAE
    input2 = deepcopy(input1)
    input2.Q_LOC = input1.Q_LOC * 1.5
    input2.BETAE = input1.BETAE * 0.8

    # Variation 3: Modified RLTS_2 and VEXB_SHEAR
    input3 = deepcopy(input1)
    input3.RLTS_2 = input1.RLTS_2 * 1.2
    input3.VEXB_SHEAR = input1.VEXB_SHEAR * 0.5

    return [input1, input2, input3]
end

# Expected values for regression tests across models and fidelity modes
const REGRESSION_EXPECTED_VALUES = Dict(
    # sat3_em_d3d_azf-1, fidelity=:TGLFNN
    ("sat3_em_d3d_azf-1", :TGLFNN, 1) => (
        ENERGY_FLUX_e = 3.612103806515152,
        ENERGY_FLUX_i = 6.160152661030723,
        PARTICLE_FLUX_e = 0.6037634688605953,
        STRESS_TOR_i = 2.5092916159367995,
    ),
    ("sat3_em_d3d_azf-1", :TGLFNN, 2) => (
        ENERGY_FLUX_e = 1.8837089002569793,
        ENERGY_FLUX_i = 2.3364969232015023,
        PARTICLE_FLUX_e = 0.1257954283852079,
        STRESS_TOR_i = 1.3059147156837425,
    ),
    ("sat3_em_d3d_azf-1", :TGLFNN, 3) => (
        ENERGY_FLUX_e = 4.632091591270553,
        ENERGY_FLUX_i = 9.48671993800187,
        PARTICLE_FLUX_e = 0.957060317337908,
        STRESS_TOR_i = 4.255046963250391,
    ),

    # sat3_em_d3d_azf-1, fidelity=:GKNN (gknne/i/g/p24 branch)
    ("sat3_em_d3d_azf-1", :GKNN, 1) => (
        ENERGY_FLUX_e = 2.326970214366194,
        ENERGY_FLUX_i = 3.8352993604857395,
        PARTICLE_FLUX_e = 0.5105086253044998,
        STRESS_TOR_i = 2.5997926311395054,
    ),
    ("sat3_em_d3d_azf-1", :GKNN, 2) => (
        ENERGY_FLUX_e = 1.61453271832532,
        ENERGY_FLUX_i = 1.30421968907876,
        PARTICLE_FLUX_e = 0.339137038337307,
        STRESS_TOR_i = 1.3007845401706544,
    ),
    ("sat3_em_d3d_azf-1", :GKNN, 3) => (
        ENERGY_FLUX_e = 3.1545384803862357,
        ENERGY_FLUX_i = 6.948598145901189,
        PARTICLE_FLUX_e = 0.9030322917041694,
        STRESS_TOR_i = 4.248745514181584,
    ),

    # sat3_em_d3d+mastu+nstx_azf-1, fidelity=:GKNN (gknn31 branch)
    ("sat3_em_d3d+mastu+nstx_azf-1", :GKNN, 1) => (
        ENERGY_FLUX_e = 2.2558135138639814,
        ENERGY_FLUX_i = 3.8954310883831766,
        PARTICLE_FLUX_e = 0.5192591528467271,
        STRESS_TOR_i = 2.7706766501101456,
    ),
    ("sat3_em_d3d+mastu+nstx_azf-1", :GKNN, 2) => (
        ENERGY_FLUX_e = 2.0174103876797473,
        ENERGY_FLUX_i = 2.137456777691702,
        PARTICLE_FLUX_e = 0.1745448150550685,
        STRESS_TOR_i = 1.5831589871630858,
    ),
    ("sat3_em_d3d+mastu+nstx_azf-1", :GKNN, 3) => (
        ENERGY_FLUX_e = 3.79254031378578,
        ENERGY_FLUX_i = 8.437753960524402,
        PARTICLE_FLUX_e = 1.2255810079609373,
        STRESS_TOR_i = 5.632476086598156,
    ),

    # sat3_em_d3d_azf-1_gkdb, fidelity=:GKNN (gknn31 + cgyro branch)
    ("sat3_em_d3d_azf-1_gkdb", :GKNN, 1) => (
        ENERGY_FLUX_e = 0.2393618170708511,
        ENERGY_FLUX_i = 4.645704500677719,
        PARTICLE_FLUX_e = -0.2002497215001747,
        STRESS_TOR_i = -1.460232726488118,
    ),
    ("sat3_em_d3d_azf-1_gkdb", :GKNN, 2) => (
        ENERGY_FLUX_e = 0.07893419857398694,
        ENERGY_FLUX_i = 1.559657642656188,
        PARTICLE_FLUX_e = -0.0020799893550763034,
        STRESS_TOR_i = -0.3851482268883565,
    ),
    ("sat3_em_d3d_azf-1_gkdb", :GKNN, 3) => (
        ENERGY_FLUX_e = 0.33133716581952266,
        ENERGY_FLUX_i = 7.81927305867567,
        PARTICLE_FLUX_e = -0.20606777890050687,
        STRESS_TOR_i = -2.505297891574754,
    ),

    # sat3_em_d3d_azf-1_withnegD, fidelity=:GKNN (gknn31 for core, gknn37 for nearedge/edge)
    # Core region: standard regression inputs have RMIN_LOC ~ 0.573 (< 0.881)
    ("sat3_em_d3d_azf-1_withnegD", :GKNN, 1) => (
        ENERGY_FLUX_e = 2.6684677958208214,
        ENERGY_FLUX_i = 3.6847767062720567,
        PARTICLE_FLUX_e = 0.4532594773065663,
        STRESS_TOR_i = 2.6189945343365966,
    ),
    ("sat3_em_d3d_azf-1_withnegD", :GKNN, 2) => (
        ENERGY_FLUX_e = 1.2031340124396976,
        ENERGY_FLUX_i = 1.2022329022959128,
        PARTICLE_FLUX_e = 0.005548728196802677,
        STRESS_TOR_i = 0.934178367714521,
    ),
    ("sat3_em_d3d_azf-1_withnegD", :GKNN, 3) => (
        ENERGY_FLUX_e = 3.170809458786961,
        ENERGY_FLUX_i = 5.73307225277411,
        PARTICLE_FLUX_e = 0.807380415310261,
        STRESS_TOR_i = 4.127336366002528,
    ),

    # sat3_em_d3d+mastu_azf-1, fidelity=:GKNN (gknn36 branch)
    ("sat3_em_d3d+mastu_azf-1", :GKNN, 1) => (
        ENERGY_FLUX_e = 2.263637009987901,
        ENERGY_FLUX_i = 3.9894442429327945,
        PARTICLE_FLUX_e = 0.6549655928795123,
        STRESS_TOR_i = 2.717778773990871,
    ),
    ("sat3_em_d3d+mastu_azf-1", :GKNN, 2) => (
        ENERGY_FLUX_e = 1.5714492824954072,
        ENERGY_FLUX_i = 2.2383117710101845,
        PARTICLE_FLUX_e = 0.14861916072222559,
        STRESS_TOR_i = 1.9524814752010464,
    ),
    ("sat3_em_d3d+mastu_azf-1", :GKNN, 3) => (
        ENERGY_FLUX_e = 4.299515290536031,
        ENERGY_FLUX_i = 8.967724726518043,
        PARTICLE_FLUX_e = 1.2741434422510884,
        STRESS_TOR_i = 5.672982151946214,
    ),
)

# Expected GKNN outputs for sat3_em_d3d_azf-1_withnegD near-edge and edge regions.
# These use the same 3 regression inputs but with RMIN_LOC overridden.
const EXPECTED_WITHNEGD_GKNN_NEAREDGE = [
    (ENERGY_FLUX_e = 7.0707026293322075,  ENERGY_FLUX_i = 9.168772659151367,   PARTICLE_FLUX_e = 0.5154563007863154, STRESS_TOR_i = 10.268832008465322),
    (ENERGY_FLUX_e = 12.076422778679625,  ENERGY_FLUX_i = 21.937627911188027,  PARTICLE_FLUX_e = 2.3654111750957756, STRESS_TOR_i = 21.95335587598116),
    (ENERGY_FLUX_e = 7.031283327361469,   ENERGY_FLUX_i = 12.345653679696769,  PARTICLE_FLUX_e = 0.8550890719132245, STRESS_TOR_i = 11.803256087757445),
]

const EXPECTED_WITHNEGD_GKNN_EDGE = [
    (ENERGY_FLUX_e = 46.09076728069341,  ENERGY_FLUX_i = 69.0036589060271,    PARTICLE_FLUX_e = 12.028206971680389, STRESS_TOR_i = 77.66283910913276),
    (ENERGY_FLUX_e = 90.81229427158553,  ENERGY_FLUX_i = 102.34896197878268,  PARTICLE_FLUX_e = 23.975390432947474, STRESS_TOR_i = 76.4223323296172),
    (ENERGY_FLUX_e = 44.76493473777022,  ENERGY_FLUX_i = 75.97364835822846,   PARTICLE_FLUX_e = 12.749389144163326, STRESS_TOR_i = 86.82843914976294),
]

# ============================================
# FINN Test Constants
# ============================================

const TEST_FINN_MODEL = "finn_sat3_d3d_withnegD"

# ============================================
# ModeID Test Constants
# ============================================

const TEST_MODEID_MODEL = "modeid_qlgyro_sat3_azf-1"
const MODEID_N_INPUTS = 34
const MODEID_N_CLASSES = 5
const MODEID_YNAMES = ["ETG", "ITG", "KBM", "MTM", "TEM"]  # alphabetical order

# Midpoint of training bounds input → expected outputs
# Captured 2026-04-01. ynames order: RLNS_1, RLNS_2, RLTS_1, RLTS_2, VEXB_SHEAR
const EXPECTED_FINN_MIDPOINT = (
    RLTS_1     =  2.7906991467656845,
    RLTS_2     =  1.70271305966974,
    RLNS_1     =  1.5691446730152112,
    VEXB_SHEAR = -0.0036021624654632954,
)

# Column-1 of matrix prediction (same input, first of three columns)
const EXPECTED_FINN_MATRIX_COL1 = (
    RLTS_1     =  2.7906991467656836,
    RLTS_2     =  1.70271305966974,
    RLNS_1     =  1.5691446730152108,
    VEXB_SHEAR = -0.003602162465463306,
)

# Model configurations for regression testing
const REGRESSION_MODEL_CONFIGS = [
    ("sat3_em_d3d_azf-1", :TGLFNN, "TGLFNN baseline"),
    ("sat3_em_d3d_azf-1", :GKNN, "GKNN gknne/i/g/p24"),
    ("sat3_em_d3d+mastu+nstx_azf-1", :GKNN, "GKNN gknn31"),
    ("sat3_em_d3d_azf-1_gkdb", :GKNN, "GKNN gknn31+cgyro"),
    ("sat3_em_d3d+mastu_azf-1", :GKNN, "GKNN gknn36"),
    ("sat3_em_d3d_azf-1_withnegD", :GKNN, "GKNN gknn31/gknn37 radial (core)"),
]
