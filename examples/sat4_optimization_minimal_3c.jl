# ============================================================================
# SAT4 Coefficient Optimization - 3-COEFFICIENT VERSION
# VERSION: 2025-12-02
#
# Optimizes 3 SAT4 coefficients (C_NORM, C_EXP, C_COEFF) with C_ETG fixed.
# ============================================================================

# ============================================================================
# COMMAND LINE ARGUMENT PARSING
# ============================================================================

using ArgParse

function parse_commandline()
    s = ArgParseSettings(description = "SAT4 coefficient optimization (C_NORM, C_EXP, C_COEFF only, C_ETG fixed)")

    @add_arg_table! s begin
        "--optimizer", "-o"
            help = "Optimizer type: samin, particleswarm, or nelder_mead"
            arg_type = String
            default = "samin"
        "--grid", "-g"
            help = "Grid search size (0=disabled, 3=coarse)"
            arg_type = Int
            default = 0
        "--shots", "-s"
            help = "Maximum number of shots to use (default: 20)"
            arg_type = Int
            default = 20
        "--test"
            help = "Enable test mode with 10 shots"
            action = :store_true
    end

    return parse_args(s)
end

# Parse arguments
args = parse_commandline()

# ============================================================================
# CONFIGURATION
# ============================================================================

const TEST_MODE = args["test"]
const N_TEST_SHOTS = 10
const N_WORKERS = TEST_MODE ? 10 : 120

# Grid search size from command line
const GRID_SIZE = args["grid"]

# Always run refinement
const RUN_REFINEMENT = true

# Maximum shots to use (from command line)
const MAX_PRODUCTION_SHOTS = args["shots"]

# ============================================================================
# NOTE: This version optimizes only 3 parameters
# C_ETG is ALWAYS fixed at its default value (1.25)
# ============================================================================

# ============================================================================
# REFINEMENT WITH CONFIGURABLE OPTIMIZER
# ============================================================================
# Optimizer options:
#   :samin        - Simulated Annealing (robust, handles noise, respects bounds)
#                   Converges when: f_calls_limit reached OR f_tol met OR iterations exhausted
#   :particleswarm - Particle Swarm Optimization (global search)
#                   Converges when: f_calls_limit reached OR f_tol met OR iterations exhausted
#   :nelder_mead  - Nelder-Mead with Fminbox (original, but can violate bounds)
#                   Converges when: f_calls_limit reached OR f_tol/x_tol met OR iterations exhausted
#
# CONVERGENCE CRITERIA (set in refine_with_optimizer function):
#   f_calls_limit - Maximum cost function evaluations (NOT individual TGLF calls!)
#                   With N shots: actual TGLF calls = f_calls_limit × N
#   f_tol         - Relative function tolerance: stops if |f_new - f_old|/|f_old| < f_tol
#   x_tol         - Parameter tolerance: stops if ||x_new - x_old|| < x_tol
#
# Parse optimizer type from command line
const OPTIMIZER_TYPE = Symbol(args["optimizer"])

# Fortran TGLF defaults
const C_NORM_DEFAULT = 1.82770384
const C_EXP_DEFAULT = 1.39786897
const C_COEFF_DEFAULT = 0.36017009
const C_ETG_DEFAULT = 1.25

bound_lo = 0.05
bound_hi = 20.0

const C_NORM_BOUNDS = (C_NORM_DEFAULT * bound_lo, C_NORM_DEFAULT * bound_hi)
const C_EXP_BOUNDS = (C_EXP_DEFAULT * bound_lo, C_EXP_DEFAULT * bound_hi)
const C_COEFF_BOUNDS = (C_COEFF_DEFAULT * bound_lo, C_COEFF_DEFAULT * bound_hi)
# C_ETG is NOT optimized in this version - always uses default

# Create unique output directory for this configuration
const OUTPUT_DIR = "sat4_results_3c_$(args["optimizer"])_grid$(args["grid"])_shots$(args["shots"])"
const DATA_FILE = "cgyrodb_master_with_caseid.json"
const VALID_INDICES_FILE = "valid_shot_indices.txt"

# ============================================================================
# WORKER SETUP
# ============================================================================

using Distributed

if nprocs() == 1
    println("Adding $N_WORKERS worker processes...")
    addprocs(N_WORKERS,
             exeflags="--project=/global/homes/t/trifolio/.julia/dev/TJLF",
             topology=:all_to_all)
    println("Workers added: $(nworkers()) workers\n")
end

# Activate TJLF dev environment on main process
using Pkg
Pkg.activate("/global/homes/t/trifolio/.julia/dev/TJLF")

# Load MINIMAL packages on main process
println("Loading packages on main process (TJLF only)...")
using TJLF, JSON, Optim, Statistics, DelimitedFiles, Printf
println("Main process ready!\n")

# Load on workers
if nworkers() > 0
    println("Loading packages on $(nworkers()) workers...")
    # Don't activate on workers - they inherit from exeflags
    @everywhere begin
        using TJLF, JSON, Optim, Statistics, DelimitedFiles, Printf
    end
    println("Workers ready!\n")

    # Share constants
    @eval @everywhere const C_NORM_DEFAULT = $C_NORM_DEFAULT
    @eval @everywhere const C_EXP_DEFAULT = $C_EXP_DEFAULT
    @eval @everywhere const C_COEFF_DEFAULT = $C_COEFF_DEFAULT
    @eval @everywhere const C_ETG_DEFAULT = $C_ETG_DEFAULT
end

# ============================================================================
# DATA LOADING
# ============================================================================

@everywhere function load_and_filter_data(json_file::String, valid_indices_file::String, max_shots::Union{Int,Nothing}=nothing)
    # Parse JSON - returns Dict with shot IDs as keys
    data_dict = JSON.parsefile(json_file)

    if !isa(data_dict, Dict)
        error("Expected JSON to contain a dictionary of shots")
    end

    # Convert dict to array - MUST match filter_valid_shots.jl indexing!
    # filter_valid_shots.jl uses collect(values(data_dict)) which is unordered
    all_shots = collect(values(data_dict))

    # Load valid shot indices (1-based indices from filter_valid_shots.jl)
    valid_indices = readdlm(valid_indices_file, Int)[:, 1]
    println("Loaded $(length(valid_indices)) pre-filtered valid shot indices from $valid_indices_file")

    # Filter to valid shots
    valid_shots = [all_shots[i] for i in valid_indices if i <= length(all_shots)]
    println("Total shots in dataset: $(length(all_shots))")
    println("Valid shots after ky spectrum filtering: $(length(valid_shots))")

    if max_shots !== nothing
        println("Limiting to first $max_shots shots for testing")
        return valid_shots[1:min(max_shots, length(valid_shots))]
    else
        return valid_shots
    end
end

# ============================================================================
# TGLF INPUT GENERATION
# ============================================================================

# Template input.tglf string (minimal CGYRO defaults)
@everywhere input_tglf_template = """
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
KYGRID_MODEL = 0
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
UNITS = 'GYRO'
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
C_B = 0.315
SIG_B = 0.34
BOUNCE_COEFF = 3.0
ZS_1 = -1.0
ZS_2 = 1.0
ZS_3 = 6.0
C_NORM = 1.82770384
C_EXP = 1.39786897
C_COEFF = 0.36017009
C_ETG = 1.25
""";

@everywhere function create_tglf_input(shot_data::Dict, c_norm, c_exp, c_coeff, c_etg)
    # Parse template
    input_lines = split(input_tglf_template, '\n')

    # Collect params to update from shot["tglf"]
    params_to_update = copy(shot_data["tglf"])

    # Override with optimization-specific settings
    params_to_update["SAT_RULE"] = 4
    params_to_update["UNITS"] = "CGYRO"
    params_to_update["FIND_WIDTH"] = false
    params_to_update["NKY"] = shot_data["input_cgyro"]["N_TOROIDAL"]
    params_to_update["KY"] = shot_data["input_cgyro"]["N_TOROIDAL"] * shot_data["input_cgyro"]["KY"]

    # Set SAT4 coefficients
    params_to_update["C_NORM"] = c_norm
    params_to_update["C_EXP"] = c_exp
    params_to_update["C_COEFF"] = c_coeff
    params_to_update["C_ETG"] = c_etg

    # Update template lines with params_to_update
    for i in 1:length(input_lines)
        m = match(r"^([A-Za-z0-9_]+)\s*=\s*(.*)", input_lines[i])
        if m !== nothing
            pname = m.captures[1]
            template_val_str = strip(m.captures[2])
            if haskey(params_to_update, pname)
                val = params_to_update[pname]
                # Handle boolean conversion
                if val === true
                    val = "true"
                elseif val === false
                    val = "false"
                end
                input_lines[i] = "$(pname) = $(val)"
            end
        end
    end

    # Write to temporary file
    filepath = tempname() * ".tglf"
    open(filepath, "w") do f
        write(f, join(input_lines, "\n"))
    end

    # Load as InputTJLF using readInput
    input_tjlf = TJLF.readInput(filepath)

    # Clean up temp file
    rm(filepath, force=true)

    return input_tjlf
end

# ============================================================================
# TGLF RUNNER
# ============================================================================

@everywhere function run_tglf_get_fluxes(shot_data::Dict, c_params::Vector{Float64}; _version::Int=1)
    c_norm, c_exp, c_coeff = c_params
    c_etg = C_ETG_DEFAULT  # Always use default for C_ETG

    try
        # Create input
        input_tjlf = create_tglf_input(shot_data, c_norm, c_exp, c_coeff, c_etg)

        # Run TJLF pipeline
        satParams = TJLF.get_sat_params(input_tjlf)
        input_tjlf.KY_SPECTRUM .= TJLF.get_ky_spectrum(input_tjlf, satParams.grad_r0)
        outputHermite = TJLF.gauss_hermite(input_tjlf)
        QL_weights, eigenvalue = TJLF.tjlf_TM(input_tjlf, satParams, outputHermite)

        # Get zonal mixing for SAT4 # DONT USE
        # sat_params_with_mixing = TJLF.get_zonal_mixing(input_tjlf, satParams, eigenvalue, QL_weights)

        # Sum over KY spectrum to get integrated fluxes
        QL_flux_out, flux_spectrum = TJLF.sum_ky_spectrum(input_tjlf, satParams, eigenvalue[:, :, 1], QL_weights)

        # Extract fluxes
        # QL_flux_out dimensions: [field (1=particle, 2=energy), species, radial_location]
        q_electrons = sum(QL_flux_out[:, 1, 2])  # Electron energy flux
        q_ions = sum(QL_flux_out[:, 2:end, 2])   # Ion energy flux (sum over all ion species)
        g_electrons = sum(QL_flux_out[:, 1, 1])  # Electron particle flux

        return (q_ions, q_electrons, g_electrons)

    catch e
        @warn "TGLF run failed" exception=e
        return (NaN, NaN, NaN)
    end
end

# ============================================================================
# LOSS COMPUTATION
# ============================================================================

@everywhere function compute_shot_loss(shot::Dict, c_params::Vector{Float64})
    # Get TGLF predictions
    pred_q_ions, pred_q_electrons, pred_g_electrons = run_tglf_get_fluxes(shot, c_params)

    # Check for NaN
    if isnan(pred_q_ions) || isnan(pred_q_electrons) || isnan(pred_g_electrons)
        return 1e6  # Large penalty for failed runs
    end

    # Get CGYRO targets
    target_q_ions = shot["flux"]["q"][1]  # Ion heat flux
    target_q_electrons = shot["flux"]["q"][2]  # Electron heat flux
    target_g_electrons = shot["flux"]["g"][2]  # Electron particle flux

    # Compute asinh-transformed errors
    err_q_ions = asinh(pred_q_ions) - asinh(target_q_ions)
    err_q_electrons = asinh(pred_q_electrons) - asinh(target_q_electrons)
    err_g_electrons = asinh(pred_g_electrons) - asinh(target_g_electrons)

    # Total loss for this shot (sum of squared errors)
    loss = err_q_ions^2 + err_q_electrons^2 + err_g_electrons^2

    return loss
end

# Diagnostic version that returns detailed breakdown
@everywhere function compute_shot_loss_detailed(shot::Dict, c_params::Vector{Float64})
    # Get TGLF predictions
    pred_q_ions, pred_q_electrons, pred_g_electrons = run_tglf_get_fluxes(shot, c_params)

    # Check for NaN
    if isnan(pred_q_ions) || isnan(pred_q_electrons) || isnan(pred_g_electrons)
        return (loss=1e6, failed=true, q_ions_loss=NaN, q_elec_loss=NaN, g_elec_loss=NaN,
                pred_q_ions=NaN, pred_q_elec=NaN, pred_g_elec=NaN,
                target_q_ions=NaN, target_q_elec=NaN, target_g_elec=NaN)
    end

    # Get CGYRO targets
    target_q_ions = shot["flux"]["q"][1]
    target_q_electrons = shot["flux"]["q"][2]
    target_g_electrons = shot["flux"]["g"][2]

    # Compute asinh-transformed errors
    err_q_ions = asinh(pred_q_ions) - asinh(target_q_ions)
    err_q_electrons = asinh(pred_q_electrons) - asinh(target_q_electrons)
    err_g_electrons = asinh(pred_g_electrons) - asinh(target_g_electrons)

    # Individual channel losses
    q_ions_loss = err_q_ions^2
    q_elec_loss = err_q_electrons^2
    g_elec_loss = err_g_electrons^2

    # Total loss
    loss = q_ions_loss + q_elec_loss + g_elec_loss

    return (loss=loss, failed=false,
            q_ions_loss=q_ions_loss, q_elec_loss=q_elec_loss, g_elec_loss=g_elec_loss,
            pred_q_ions=pred_q_ions, pred_q_elec=pred_q_electrons, pred_g_elec=pred_g_electrons,
            target_q_ions=target_q_ions, target_q_elec=target_q_electrons, target_g_elec=target_g_electrons)
end

function evaluate_cost_function(c_params::Vector{Float64}, shot_list::Vector)
    # Parallel reduction over shots
    total_loss = @distributed (+) for shot in shot_list
        compute_shot_loss(shot, c_params)
    end
    return total_loss
end

function evaluate_diagnostics(c_params::Vector{Float64}, shot_list::Vector, label::String="")
    """
    Evaluate parameters with detailed per-shot diagnostics.
    Returns per-shot loss breakdown and prints statistics.
    """
    println("\n" * "="^70)
    println("DIAGNOSTIC EVALUATION: $label")
    println("="^70)
    println("Parameters: C_NORM=$(c_params[1]), C_EXP=$(c_params[2]), C_COEFF=$(c_params[3]) (C_ETG fixed at $C_ETG_DEFAULT)")

    # Collect detailed results from all shots
    results = pmap(shot -> compute_shot_loss_detailed(shot, c_params), shot_list)

    # Separate successful and failed shots
    successful = [r for r in results if !r.failed]
    n_failed = count(r -> r.failed, results)

    if n_failed > 0
        println("\nWARNING: $n_failed shots failed!")
    end

    if isempty(successful)
        println("ERROR: All shots failed!")
        return nothing
    end

    # Extract channel-wise losses
    q_ions_losses = [r.q_ions_loss for r in successful]
    q_elec_losses = [r.q_elec_loss for r in successful]
    g_elec_losses = [r.g_elec_loss for r in successful]
    total_losses = [r.loss for r in successful]

    # Compute statistics
    println("\n" * "-"^70)
    println("LOSS STATISTICS ($(length(successful)) shots)")
    println("-"^70)
    println(@sprintf("%-20s %12s %12s %12s %12s", "Channel", "Mean", "Median", "Min", "Max"))
    println("-"^70)
    println(@sprintf("%-20s %12.4f %12.4f %12.4f %12.4f", "q_ions",
            mean(q_ions_losses), median(q_ions_losses), minimum(q_ions_losses), maximum(q_ions_losses)))
    println(@sprintf("%-20s %12.4f %12.4f %12.4f %12.4f", "q_electrons",
            mean(q_elec_losses), median(q_elec_losses), minimum(q_elec_losses), maximum(q_elec_losses)))
    println(@sprintf("%-20s %12.4f %12.4f %12.4f %12.4f", "g_electrons",
            mean(g_elec_losses), median(g_elec_losses), minimum(g_elec_losses), maximum(g_elec_losses)))
    println("-"^70)
    println(@sprintf("%-20s %12.4f %12.4f %12.4f %12.4f", "TOTAL",
            mean(total_losses), median(total_losses), minimum(total_losses), maximum(total_losses)))
    println("-"^70)

    # Calculate total loss
    total_loss = sum(total_losses)
    println("\nTotal loss: $total_loss")

    # Channel contribution percentages
    total_q_ions = sum(q_ions_losses)
    total_q_elec = sum(q_elec_losses)
    total_g_elec = sum(g_elec_losses)

    println("\nChannel contributions:")
    println(@sprintf("  q_ions:      %6.2f%% (%.4f)", 100*total_q_ions/total_loss, total_q_ions))
    println(@sprintf("  q_electrons: %6.2f%% (%.4f)", 100*total_q_elec/total_loss, total_q_elec))
    println(@sprintf("  g_electrons: %6.2f%% (%.4f)", 100*total_g_elec/total_loss, total_g_elec))

    # Find worst shots
    println("\n" * "-"^70)
    println("WORST 10 SHOTS")
    println("-"^70)
    sorted_indices = sortperm(total_losses, rev=true)
    println(@sprintf("%-6s %12s %12s %12s %12s", "Rank", "Total Loss", "q_ions", "q_elec", "g_elec"))
    println("-"^70)
    for i in 1:min(10, length(sorted_indices))
        idx = sorted_indices[i]
        r = successful[idx]
        println(@sprintf("%-6d %12.4f %12.4f %12.4f %12.4f", i, r.loss, r.q_ions_loss, r.q_elec_loss, r.g_elec_loss))
    end

    # Find best shots
    println("\n" * "-"^70)
    println("BEST 10 SHOTS")
    println("-"^70)
    println(@sprintf("%-6s %12s %12s %12s %12s", "Rank", "Total Loss", "q_ions", "q_elec", "g_elec"))
    println("-"^70)
    for i in 1:min(10, length(sorted_indices))
        idx = sorted_indices[end-i+1]
        r = successful[idx]
        println(@sprintf("%-6d %12.4f %12.4f %12.4f %12.4f", i, r.loss, r.q_ions_loss, r.q_elec_loss, r.g_elec_loss))
    end

    println("="^70)

    return (results=results, total_loss=total_loss)
end

# ============================================================================
# GRID SEARCH
# ============================================================================

function grid_search(shot_list::Vector, grid_size::Int)
    println("\n" * "="^70)
    println("GRID SEARCH")
    println("="^70)
    println("Optimizing 3 parameters: C_NORM, C_EXP, C_COEFF")
    println("Grid size: $grid_size points per dimension")
    println("Total evaluations: $(grid_size^3)")

    # Create grid points for 3 parameters
    c_norm_grid = range(C_NORM_BOUNDS[1], C_NORM_BOUNDS[2], length=grid_size)
    c_exp_grid = range(C_EXP_BOUNDS[1], C_EXP_BOUNDS[2], length=grid_size)
    c_coeff_grid = range(C_COEFF_BOUNDS[1], C_COEFF_BOUNDS[2], length=grid_size)

    best_loss = Inf
    best_params = [C_NORM_DEFAULT, C_EXP_DEFAULT, C_COEFF_DEFAULT]

    eval_count = 0
    total_evals = grid_size^3
    start_time = time()

    # Triple nested loop over 3 parameters
    for c_norm in c_norm_grid
        for c_exp in c_exp_grid
            for c_coeff in c_coeff_grid
                eval_count += 1
                params = [c_norm, c_exp, c_coeff]
                loss = evaluate_cost_function(params, shot_list)

                if loss < best_loss
                    best_loss = loss
                    best_params = copy(params)
                    println("New best at eval $eval_count/$total_evals: loss = $loss")
                    println("  C_NORM=$c_norm, C_EXP=$c_exp, C_COEFF=$c_coeff")
                end

                # Progress update every 100 evaluations
                if eval_count % 100 == 0
                    elapsed = time() - start_time
                    rate = eval_count / elapsed
                    remaining = (total_evals - eval_count) / rate
                    println("Progress: $eval_count/$total_evals ($(round(100*eval_count/total_evals, digits=1))%), " *
                            "$(round(rate, digits=2)) eval/s, ~$(round(remaining/60, digits=1)) min remaining")
                end
            end
        end
    end

    elapsed_time = time() - start_time
    println("\nGrid search complete!")
    println("Best loss: $best_loss")
    println("Best parameters: C_NORM=$(best_params[1]), C_EXP=$(best_params[2]), C_COEFF=$(best_params[3])")
    println("Time: $(round(elapsed_time/60, digits=2)) minutes")

    return best_params, best_loss
end


function refine_with_optimizer(initial_params::Vector{Float64}, shot_list::Vector)
    println("\n" * "="^70)
    println("OPTIMIZER REFINEMENT")
    println("="^70)
    println("Optimizer: $OPTIMIZER_TYPE")
    println("Optimizing 3 parameters: C_NORM, C_EXP, C_COEFF")
    println("Starting from: C_NORM=$(initial_params[1]), C_EXP=$(initial_params[2]), C_COEFF=$(initial_params[3])")

    # Define objective function
    function objective(params)
        return evaluate_cost_function(params, shot_list)
    end

    # Set up bounds for 3 parameters
    lower = [C_NORM_BOUNDS[1], C_EXP_BOUNDS[1], C_COEFF_BOUNDS[1]]
    upper = [C_NORM_BOUNDS[2], C_EXP_BOUNDS[2], C_COEFF_BOUNDS[2]]

    # Nudge starting point off boundaries to avoid Fminbox issues
    # (Fminbox can fail or leak when starting exactly on boundary)
    epsilon = 0.05  # margin from boundary
    nudged_params = copy(initial_params)
    param_names = ["C_NORM", "C_EXP", "C_COEFF"]
    for i in 1:3
        if nudged_params[i] >= upper[i] * (1 - 0.001)  # Within 0.1% of upper bound
            nudged_params[i] = upper[i] * (1 - epsilon)
            println("  Nudging $(param_names[i]) off upper boundary: $(initial_params[i]) -> $(nudged_params[i])")
        elseif nudged_params[i] <= lower[i] * (1 + 0.001)  # Within 0.1% of lower bound
            nudged_params[i] = lower[i] * (1 + epsilon)
            println("  Nudging $(param_names[i]) off lower boundary: $(initial_params[i]) -> $(nudged_params[i])")
        end
    end

    # Run optimization based on selected optimizer
    start_time = time()

    if OPTIMIZER_TYPE == :samin
        # Simulated Annealing - excellent for noisy, expensive functions
        # Respects bounds natively, robust to local minima
        # With 20 shots: f_calls_limit=1000 ≈ 2.7 hours wall clock, 20k total TGLF calls
        result = optimize(objective, lower, upper, nudged_params, SAMIN(),
                         Optim.Options(
                             iterations = 10^5,          # Max SA iterations (usually converges earlier)
                             f_calls_limit = 1000,       # Max cost function evals (CRITICAL for expensive funcs)
                             f_reltol = 1e-3,               # Stop if relative improvement < 0.1%
                             show_trace = true,
                             show_every = 100
                         ))

    elseif OPTIMIZER_TYPE == :particleswarm
        # Particle Swarm - good global optimizer for box-constrained problems
        # With 20 shots: f_calls_limit=500 ≈ 1.4 hours wall clock, 10k total TGLF calls
        result = optimize(objective, lower, upper, nudged_params, ParticleSwarm(n_particles=20),
                         Optim.Options(
                             iterations = 1000,
                             f_calls_limit = 500,        # Reasonable limit for particle swarm
                             f_reltol = 1e-3,               # Stop if relative improvement < 0.1%
                             show_trace = true,
                             show_every = 10
                         ))

    elseif OPTIMIZER_TYPE == :nelder_mead
        # Original Nelder-Mead with Fminbox (can violate bounds as seen in logs)
        # With 20 shots: f_calls_limit=500 ≈ 1.4 hours wall clock, 10k total TGLF calls
        result = optimize(objective, lower, upper, nudged_params, Fminbox(NelderMead()),
                         Optim.Options(
                             iterations = 1000,
                             f_calls_limit = 500,
                             f_reltol = 1e-4,               # Stop if relative improvement < 0.01%
                             x_tol = 1e-4,               # Stop if parameter changes tiny
                             show_trace = true,
                             show_every = 10
                         ))

    else
        error("Unknown optimizer type: $OPTIMIZER_TYPE")
    end

    elapsed_time = time() - start_time

    refined_params = Optim.minimizer(result)
    refined_loss = Optim.minimum(result)

    # Check if bounds were violated
    bounds_violated = any(refined_params .< lower) || any(refined_params .> upper)
    if bounds_violated
        println("\nWARNING: Refined parameters violate bounds!")
        for i in 1:3
            if refined_params[i] < lower[i]
                println("  $(param_names[i]): $(refined_params[i]) < lower bound $(lower[i])")
            elseif refined_params[i] > upper[i]
                println("  $(param_names[i]): $(refined_params[i]) > upper bound $(upper[i])")
            end
        end
    end

    println("\nRefinement complete!")
    println("Refined loss: $refined_loss")
    println("Refined parameters: C_NORM=$(refined_params[1]), C_EXP=$(refined_params[2]), C_COEFF=$(refined_params[3])")
    println("Time: $(round(elapsed_time/60, digits=2)) minutes")
    println("Iterations: $(Optim.iterations(result))")
    println("Converged: $(Optim.converged(result))")

    return refined_params, refined_loss
end

# ============================================================================
# BASELINE EVALUATION
# ============================================================================

function evaluate_baseline(shot_list::Vector)
    println("\n" * "="^70)
    println("BASELINE EVALUATION (Fortran defaults)")
    println("="^70)

    baseline_params = [C_NORM_DEFAULT, C_EXP_DEFAULT, C_COEFF_DEFAULT]
    baseline_loss = evaluate_cost_function(baseline_params, shot_list)

    println("Baseline loss: $baseline_loss")
    println("Parameters: C_NORM=$C_NORM_DEFAULT, C_EXP=$C_EXP_DEFAULT, C_COEFF=$C_COEFF_DEFAULT (C_ETG=$C_ETG_DEFAULT fixed)")

    return baseline_loss
end

# ============================================================================
# RESULTS SAVING
# ============================================================================

function save_results(baseline_loss, grid_params, grid_loss, refined_params, refined_loss, shot_list)
    mkpath(OUTPUT_DIR)

    # Save summary
    open(joinpath(OUTPUT_DIR, "optimization_summary.txt"), "w") do f
        println(f, "SAT4 Coefficient Optimization Results (3-parameter version)")
        println(f, "="^70)
        println(f, "\nNumber of shots: $(length(shot_list))")
        println(f, "\nBaseline (Fortran defaults):")
        println(f, "  C_NORM = $C_NORM_DEFAULT")
        println(f, "  C_EXP = $C_EXP_DEFAULT")
        println(f, "  C_COEFF = $C_COEFF_DEFAULT")
        println(f, "  C_ETG = $C_ETG_DEFAULT (FIXED)")
        println(f, "  Loss = $baseline_loss")
        println(f, "\nGrid search best:")
        println(f, "  C_NORM = $(grid_params[1])")
        println(f, "  C_EXP = $(grid_params[2])")
        println(f, "  C_COEFF = $(grid_params[3])")
        println(f, "  C_ETG = $C_ETG_DEFAULT (FIXED)")
        println(f, "  Loss = $grid_loss")
        println(f, "\nRefined ($OPTIMIZER_TYPE):")
        println(f, "  C_NORM = $(refined_params[1])")
        println(f, "  C_EXP = $(refined_params[2])")
        println(f, "  C_COEFF = $(refined_params[3])")
        println(f, "  C_ETG = $C_ETG_DEFAULT (FIXED)")
        println(f, "  Loss = $refined_loss")
        println(f, "\nImprovement:")
        println(f, "  Baseline → Grid: $(round(100*(baseline_loss - grid_loss)/baseline_loss, digits=2))%")
        println(f, "  Grid → Refined: $(round(100*(grid_loss - refined_loss)/grid_loss, digits=2))%")
        println(f, "  Baseline → Refined: $(round(100*(baseline_loss - refined_loss)/baseline_loss, digits=2))%")
    end

    # Save parameters as CSV
    open(joinpath(OUTPUT_DIR, "optimized_parameters.csv"), "w") do f
        println(f, "parameter,baseline,grid_search,refined")
        println(f, "C_NORM,$C_NORM_DEFAULT,$(grid_params[1]),$(refined_params[1])")
        println(f, "C_EXP,$C_EXP_DEFAULT,$(grid_params[2]),$(refined_params[2])")
        println(f, "C_COEFF,$C_COEFF_DEFAULT,$(grid_params[3]),$(refined_params[3])")
        println(f, "C_ETG,$C_ETG_DEFAULT,$C_ETG_DEFAULT,$C_ETG_DEFAULT")
    end

    println("\nResults saved to $OUTPUT_DIR/")
end

# ============================================================================
# MAIN EXECUTION
# ============================================================================

function main()
    println("\n" * "="^70)
    println("SAT4 COEFFICIENT OPTIMIZATION (3-PARAMETER VERSION)")
    println("="^70)
    println("Optimizing: C_NORM, C_EXP, C_COEFF (C_ETG fixed at $C_ETG_DEFAULT)")
    println("Mode: $(TEST_MODE ? "TEST ($N_TEST_SHOTS shots)" : "PRODUCTION")")
    println("Workers: $(nworkers())")
    if GRID_SIZE > 0
        println("Grid search: $GRID_SIZE per dimension ($(GRID_SIZE^3) total evaluations)")
    else
        println("Grid search: DISABLED (starting from defaults)")
    end
    println("Optimizer refinement: $(RUN_REFINEMENT ? "ENABLED ($OPTIMIZER_TYPE)" : "DISABLED")")
    if !TEST_MODE && MAX_PRODUCTION_SHOTS !== nothing
        println("Shot limit: $MAX_PRODUCTION_SHOTS")
    end

    # Load data
    println("\nLoading data from $DATA_FILE using valid indices from $VALID_INDICES_FILE...")
    max_shots = TEST_MODE ? N_TEST_SHOTS : MAX_PRODUCTION_SHOTS
    shot_list = load_and_filter_data(DATA_FILE, VALID_INDICES_FILE, max_shots)
    println("Using $(length(shot_list)) shots for optimization")

    # Check if we have any valid shots
    if length(shot_list) == 0
        error("No valid shots found! Cannot proceed with optimization.")
    end

    # Evaluate baseline
    baseline_loss = evaluate_baseline(shot_list)

    # Detailed diagnostic for baseline
    baseline_params = [C_NORM_DEFAULT, C_EXP_DEFAULT, C_COEFF_DEFAULT]
    evaluate_diagnostics(baseline_params, shot_list, "BASELINE (Fortran defaults)")

    # Grid search (skip if GRID_SIZE == 0)
    if GRID_SIZE > 0
        grid_params, grid_loss = grid_search(shot_list, GRID_SIZE)
        # Detailed diagnostic for best grid result
        evaluate_diagnostics(grid_params, shot_list, "GRID SEARCH BEST")
    else
        println("\nSkipping grid search (GRID_SIZE=0), will start optimizer from defaults")
        grid_params = [C_NORM_DEFAULT, C_EXP_DEFAULT, C_COEFF_DEFAULT]
        grid_loss = baseline_loss
    end

    # Refine with optimizer (optional)
    if RUN_REFINEMENT
        refined_params, refined_loss = refine_with_optimizer(grid_params, shot_list)
        # Detailed diagnostic for refined result
        evaluate_diagnostics(refined_params, shot_list, "OPTIMIZER REFINED")
    else
        println("\nSkipping optimizer refinement (RUN_REFINEMENT=false)")
        refined_params = grid_params
        refined_loss = grid_loss
    end

    # Save results
    save_results(baseline_loss, grid_params, grid_loss, refined_params, refined_loss, shot_list)

    # Final summary
    println("\n" * "="^70)
    println("FINAL RESULTS")
    println("="^70)
    println("Baseline loss:      $baseline_loss")
    println("Grid search loss:   $grid_loss")
    println("Refined loss:       $refined_loss")
    println("\nImprovement: $(round(100*(baseline_loss - refined_loss)/baseline_loss, digits=2))%")
    println("\nOptimized parameters:")
    println("  C_NORM  = $(refined_params[1]) (default: $C_NORM_DEFAULT)")
    println("  C_EXP   = $(refined_params[2]) (default: $C_EXP_DEFAULT)")
    println("  C_COEFF = $(refined_params[3]) (default: $C_COEFF_DEFAULT)")
    println("  C_ETG   = $C_ETG_DEFAULT (FIXED)")
end

# Run main
main()
