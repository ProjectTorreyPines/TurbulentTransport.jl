# Test fixtures and helper functions for TurbulentTransport tests

# Complete input.tglf content from example notebook
const TEST_INPUT_TGLF_CONTENT = """
ADIABATIC_ELEC = .false.
ALPHA_E = 1.0
ALPHA_MACH = 0.0
ALPHA_P = 1.0
ALPHA_QUENCH = 0
ALPHA_ZF = -1.0
AS_1 = 1.0
AS_2 = 0.784867
AS_3 = 0.0302081
BETAE = 0.00362972
BETA_LOC = 0.0
DAMP_PSI = 0.0
DAMP_SIG = 0.0
DEBYE = 0.0217677
DEBYE_FACTOR = 1.0
DELTA_LOC = 0.0681444
DRMAJDX_LOC = -0.189065
DRMINDX_LOC = 1.0
DZMAJDX_LOC = 0.00278328
ETG_FACTOR = 1.25
FILTER = 2.0
FIND_WIDTH = .true.
GCHAT = 1.0
GHAT = 1.0
GRADB_FACTOR = 0.0
IBRANCH = -1
IFLUX = .true.
KAPPA_LOC = 1.40438
KX0_LOC = 0.0
KY = 0.3
KYGRID_MODEL = 4
LINSKER_FACTOR = 0.0
MASS_1 = 0.000272445
MASS_2 = 1.0
MASS_3 = 6.0
NBASIS_MAX = 6
NBASIS_MIN = 2
NEW_EIKONAL = .true.
NKY = 12
NMODES = 2
NS = 3
NWIDTH = 21
NXGRID = 16
PARK = 1.0
P_PRIME_LOC = -0.00355359
Q_LOC = 2.00545
Q_PRIME_LOC = 14.7947
RLNP_CUTOFF = 18.0
RLNS_1 = 0.513787
RLNS_2 = 0.758616
RLNS_3 = -0.872531
RLTS_1 = 2.03987
RLTS_2 = 2.20153
RLTS_3 = 2.20153
RMAJ_LOC = 2.86212
RMIN_LOC = 0.573129
SAT_RULE = 3
SIGN_BT = -1
SIGN_IT = 1
S_DELTA_LOC = 0.116297
S_KAPPA_LOC = 0.125574
S_ZETA_LOC = -0.0258657
TAUS_1 = 1.0
TAUS_2 = 1.39296
TAUS_3 = 1.39296
THETA_TRAPPED = 0.7
UNITS = GYRO
USE_AVE_ION_GRID = .false.
USE_BISECTION = .true.
USE_BPAR = .true.
USE_BPER = .true.
USE_INBOARD_DETRAPPED = .false.
USE_MHD_RULE = .false.
VEXB_SHEAR = 0.080234
VPAR_1 = 0.419061
VPAR_2 = 0.419061
VPAR_3 = 0.419061
VPAR_MODEL = 0
VPAR_SHEAR_1 = 0.803536
VPAR_SHEAR_2 = 0.803536
VPAR_SHEAR_3 = 0.803536
VPAR_SHEAR_MODEL = 1
WDIA_TRAPPED = 1.0
WD_ZERO = 0.1
WIDTH = 1.65
WIDTH_MIN = 0.3
XNUE = 0.0948099
XNU_FACTOR = 1.0
XNU_MODEL = 3
ZEFF = 1.90624
ZETA_LOC = -0.0148888
ZMAJ_LOC = -0.0576768
ZS_1 = -1.0
ZS_2 = 1.0
ZS_3 = 6.0
"""

# Known good model filenames for testing
const TEST_MODEL_SINGLE = "sat3_em_d3d_azf-1"
const TEST_MODEL_ENSEMBLE = "sat3_em_d3d_azf-1"  # This is an ensemble model

# Expected field values after loading TEST_INPUT_TGLF_CONTENT
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

# Helper function to create a temporary input.tglf file
function create_temp_input_tglf(content::String=TEST_INPUT_TGLF_CONTENT)
    tmpdir = mktempdir()
    filepath = joinpath(tmpdir, "input.tglf")
    open(filepath, "w") do f
        write(f, content)
    end
    return filepath, tmpdir
end

# Helper function to cleanup temporary directory
function cleanup_temp_dir(tmpdir::String)
    rm(tmpdir; force=true, recursive=true)
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
