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
    TurbulentTransport.load(InputTGLF{Float64}(), SAMPLE_INPUT_PATH)
end

# The sample input is SAT_RULE=3 in GYRO units. TJLF 2 defines SAT2/3 in CGYRO units
# only and rejects the combination as soon as an InputTJLF is built from it
# (`checkInput` runs inside `update_input_tjlf!`, before any presets), so tests that
# need an InputTJLF must switch the units first — exactly what `apply_presets!` and
# the run_tglf/run_tjlf paths do at runtime.
function load_sample_input_cgyro()
    input_tglf = load_sample_input()
    input_tglf.UNITS = "CGYRO"
    return input_tglf
end

load_sample_input_tjlf() = InputTJLF{Float64}(load_sample_input_cgyro())

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

# flux_array expected outputs for sat0_em_d3d single model (first model in ensemble = lowest BSON key, sorted in dict2ens)
const EXPECTED_FLUX_ARRAY_SINGLE = [0.020119015925778916, 0.23107953401865972, 0.7410156911426711, 0.6620782527214439]

# GKNN correction model (sat3_em_d3d_azf-1_gknne24) expected outputs
# These models have ynames of length 2, and fidelity=:GKNN outputs div(ynames, 2) = 1 value
const EXPECTED_GKNN_MODEL_SINGLE = [1.0372520644956456]
const EXPECTED_GKNN_MODEL_ENSEMBLE = [0.9496942106346887]

# flux_array expected outputs for sat3_em_d3d_azf-1 ensemble
const EXPECTED_FLUX_ARRAY_ENSEMBLE = [-0.019443083518266357, 0.04080973584865324, 0.23715745308655625, 0.07256906359102389]

# run_tglfnn expected outputs with sample_input.tglf
const EXPECTED_RUN_TGLFNN_SAT3 = (
    ENERGY_FLUX_e = 2.5022810307621777,
    ENERGY_FLUX_i = 4.0710317823476965,
    PARTICLE_FLUX_e = 0.31299530330023073,
    STRESS_TOR_i = 2.6504124282780994,
)
# Same pin captured before v1.4 with SIGN_BT=+1 (what the legacy model was fed pre-1.4)
const EXPECTED_RUN_TGLFNN_SAT3_SIGNBT_PLUS = (
    ENERGY_FLUX_e = 3.612103806515151,
    ENERGY_FLUX_i = 6.160152661030724,
    PARTICLE_FLUX_e = 0.6037634688605953,
    STRESS_TOR_i = 2.5092916159367995,
)

const EXPECTED_RUN_TGLFNN_SAT2 = (
    ENERGY_FLUX_e = 3.181768974789443,
    ENERGY_FLUX_i = 6.419348689658621,
    PARTICLE_FLUX_e = 0.5128159266571789,
    STRESS_TOR_i = 4.740717871630627,
)

const EXPECTED_RUN_TGLFNN_SAT3_GKNN = (
    ENERGY_FLUX_e = 1.7108496154701096,
    ENERGY_FLUX_i = 2.928802971436727,
    PARTICLE_FLUX_e = 0.3873157614655667,
    STRESS_TOR_i = 2.816091554622234,
)

# ============================================
# Regression Test Expected Values
# Captured on 2025-12-06 for single/vector equivalence; recaptured 2026-09-10 for the
# VEXB_SHEAR sign-convention semantics of v1.4 (sample_input.tglf has SIGN_BT=-1, so
# every :legacy model now sees -VEXB_SHEAR; :tgyro GKNN nets see +VEXB_SHEAR)
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
        ENERGY_FLUX_e = 2.5022810307621777,
        ENERGY_FLUX_i = 4.0710317823476965,
        PARTICLE_FLUX_e = 0.31299530330023073,
        STRESS_TOR_i = 2.6504124282780994,
    ),
    ("sat3_em_d3d_azf-1", :TGLFNN, 2) => (
        ENERGY_FLUX_e = 1.1464155090787584,
        ENERGY_FLUX_i = 1.1823629172426628,
        PARTICLE_FLUX_e = 0.03001626942527181,
        STRESS_TOR_i = 0.8699248670066521,
    ),
    ("sat3_em_d3d_azf-1", :TGLFNN, 3) => (
        ENERGY_FLUX_e = 3.914423003478955,
        ENERGY_FLUX_i = 7.970802164755081,
        PARTICLE_FLUX_e = 0.7549948560122873,
        STRESS_TOR_i = 4.3881647898852645,
    ),

    # sat3_em_d3d_azf-1, fidelity=:GKNN (gknne/i/g/p24 branch)
    ("sat3_em_d3d_azf-1", :GKNN, 1) => (
        ENERGY_FLUX_e = 1.7108496154701096,
        ENERGY_FLUX_i = 2.928802971436727,
        PARTICLE_FLUX_e = 0.3873157614655667,
        STRESS_TOR_i = 2.816091554622234,
    ),
    ("sat3_em_d3d_azf-1", :GKNN, 2) => (
        ENERGY_FLUX_e = 1.1305235765710857,
        ENERGY_FLUX_i = 0.705612949163249,
        PARTICLE_FLUX_e = 0.13556179355846043,
        STRESS_TOR_i = 0.9570636383545278,
    ),
    ("sat3_em_d3d_azf-1", :GKNN, 3) => (
        ENERGY_FLUX_e = 2.836397865746182,
        ENERGY_FLUX_i = 6.235097819256284,
        PARTICLE_FLUX_e = 0.8049915364393726,
        STRESS_TOR_i = 4.502481613797547,
    ),

    # sat3_em_d3d+mastu+nstx_azf-1, fidelity=:GKNN (gknn31 branch)
    ("sat3_em_d3d+mastu+nstx_azf-1", :GKNN, 1) => (
        ENERGY_FLUX_e = 1.6057286290104413,
        ENERGY_FLUX_i = 2.5658085806948026,
        PARTICLE_FLUX_e = 0.3039934201281886,
        STRESS_TOR_i = 2.1057627687179292,
    ),
    ("sat3_em_d3d+mastu+nstx_azf-1", :GKNN, 2) => (
        ENERGY_FLUX_e = 1.8880492898066112,
        ENERGY_FLUX_i = 1.0761109235060802,
        PARTICLE_FLUX_e = -0.04401874299202968,
        STRESS_TOR_i = 1.016548970697034,
    ),
    ("sat3_em_d3d+mastu+nstx_azf-1", :GKNN, 3) => (
        ENERGY_FLUX_e = 3.234835705381657,
        ENERGY_FLUX_i = 7.280002679298558,
        PARTICLE_FLUX_e = 1.01798485482718,
        STRESS_TOR_i = 5.187377472078515,
    ),

    # sat3_em_d3d_azf-1_gkdb, fidelity=:GKNN (gknn31 + cgyro branch)
    ("sat3_em_d3d_azf-1_gkdb", :GKNN, 1) => (
        ENERGY_FLUX_e = 0.18060817026130743,
        ENERGY_FLUX_i = 3.8590266689462274,
        PARTICLE_FLUX_e = -0.13249063806797146,
        STRESS_TOR_i = -2.0000705489486386,
    ),
    ("sat3_em_d3d_azf-1_gkdb", :GKNN, 2) => (
        ENERGY_FLUX_e = 0.05336854092681487,
        ENERGY_FLUX_i = 0.9973791034564599,
        PARTICLE_FLUX_e = 0.010773943426737716,
        STRESS_TOR_i = -0.3259587358358478,
    ),
    ("sat3_em_d3d_azf-1_gkdb", :GKNN, 3) => (
        ENERGY_FLUX_e = 0.3065487347792778,
        ENERGY_FLUX_i = 7.495995814992572,
        PARTICLE_FLUX_e = -0.17521143467070988,
        STRESS_TOR_i = -2.989582771207009,
    ),

    # sat3_em_d3d_azf-1_withnegD, fidelity=:GKNN (gknn31 for core, gknn37 for nearedge/edge)
    # Core region: standard regression inputs have RMIN_LOC ~ 0.573 (< 0.881)
    ("sat3_em_d3d_azf-1_withnegD", :GKNN, 1) => (
        ENERGY_FLUX_e = 2.0441457451405443,
        ENERGY_FLUX_i = 3.05197556845397,
        PARTICLE_FLUX_e = 0.3127717369389831,
        STRESS_TOR_i = 3.455981657357479,
    ),
    ("sat3_em_d3d_azf-1_withnegD", :GKNN, 2) => (
        ENERGY_FLUX_e = 0.8415245603819719,
        ENERGY_FLUX_i = 0.7694785525872214,
        PARTICLE_FLUX_e = -0.028621644422678033,
        STRESS_TOR_i = 0.8073217886815465,
    ),
    ("sat3_em_d3d_azf-1_withnegD", :GKNN, 3) => (
        ENERGY_FLUX_e = 2.922706922544521,
        ENERGY_FLUX_i = 5.468578856510341,
        PARTICLE_FLUX_e = 0.7315935431487315,
        STRESS_TOR_i = 4.780893765420086,
    ),

    # sat3_em_d3d+mastu_azf-1, fidelity=:GKNN (gknn36 branch)
    ("sat3_em_d3d+mastu_azf-1", :GKNN, 1) => (
        ENERGY_FLUX_e = 1.6401607290583815,
        ENERGY_FLUX_i = 3.0080087254003236,
        PARTICLE_FLUX_e = 0.37580845478818553,
        STRESS_TOR_i = 2.8161530024783903,
    ),
    ("sat3_em_d3d+mastu_azf-1", :GKNN, 2) => (
        ENERGY_FLUX_e = 1.930473635015596,
        ENERGY_FLUX_i = 1.4232614410771292,
        PARTICLE_FLUX_e = -0.03208218139146896,
        STRESS_TOR_i = 1.4459849657116663,
    ),
    ("sat3_em_d3d+mastu_azf-1", :GKNN, 3) => (
        ENERGY_FLUX_e = 4.524206748434979,
        ENERGY_FLUX_i = 8.383671756005539,
        PARTICLE_FLUX_e = 1.1125355300509006,
        STRESS_TOR_i = 6.154043191277987,
    ),
)

# Expected GKNN outputs for sat3_em_d3d_azf-1_withnegD near-edge and edge regions.
# These use the same 3 regression inputs but with RMIN_LOC overridden.
const EXPECTED_WITHNEGD_GKNN_NEAREDGE = [
    (ENERGY_FLUX_e = 6.35252037682603,  ENERGY_FLUX_i = 7.601011260699546,  PARTICLE_FLUX_e = 0.406510689827285, STRESS_TOR_i = 9.109841774620573),
    (ENERGY_FLUX_e = 12.447401373926079,  ENERGY_FLUX_i = 18.35740309754034,  PARTICLE_FLUX_e = 1.8302181632536318, STRESS_TOR_i = 19.266176194018893),
    (ENERGY_FLUX_e = 6.481949776701626,  ENERGY_FLUX_i = 11.370015500689936,  PARTICLE_FLUX_e = 0.7791756115185693, STRESS_TOR_i = 11.12035966667817),
]

const EXPECTED_WITHNEGD_GKNN_EDGE = [
    (ENERGY_FLUX_e = 52.24785774657087,  ENERGY_FLUX_i = 65.27704564923825,  PARTICLE_FLUX_e = 10.837458696559732, STRESS_TOR_i = 88.60742510477147),
    (ENERGY_FLUX_e = 96.09448698049593,  ENERGY_FLUX_i = 98.33961784161689,  PARTICLE_FLUX_e = 22.374125669172734, STRESS_TOR_i = 95.32580187505751),
    (ENERGY_FLUX_e = 48.49986372917798,  ENERGY_FLUX_i = 74.10251438679151,  PARTICLE_FLUX_e = 12.162910846136697, STRESS_TOR_i = 92.99622802948817),
]

# ============================================
# FINN Test Constants
# ============================================

const TEST_FINN_MODEL = "finn_sat3_d3d_withnegD"

# ============================================
# ModeID Test Constants
# ============================================

const TEST_MODEID_MODEL = "modeid_qlgyro_sat3_azf-1"
const TEST_MODEID_TGLF_MODEL = "modeid_tglf_sat3_azf-1_d3d"   # TGLF-labelled DIII-D model, legacy VEXB sign
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
