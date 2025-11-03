using Distributed

# ============================================================================
# IMPORTANT: IF YOU EDIT @everywhere FUNCTIONS, YOU MUST RESTART JULIA!
# ============================================================================
# Workers load function definitions when this script first runs. If you edit
# any @everywhere function and re-run this file WITHOUT restarting Julia,
# the workers will still use the OLD code, causing mysterious errors.
#
# To apply changes to @everywhere functions:
#   1. Exit Julia completely
#   2. Restart Julia
#   3. Re-run this script
# ============================================================================

# Add worker processes FIRST before loading packages
# You can modify the number based on available cores
if nprocs() == 1
    n_workers = 70  # For full nested parallelism: 4 shots × ~20 ky points
    # Important: Pass the project environment to workers via exeflags
    # Use all_to_all topology for nested parallelism (required for parallel_mode="both")
    # This allows workers to communicate with each other, not just the master
    println("Adding $n_workers worker processes...")
    addprocs(n_workers,
             exeflags="--project=/global/homes/t/trifolio/.julia/dev/TurbulentTransport",
             topology=:all_to_all)  # Changed from :master_worker for nested parallelism
    println("Added $n_workers worker processes\n")
end

# Now load packages (must be at top level, after addprocs)
using Pkg
Pkg.activate("/global/homes/t/trifolio/.julia/dev/TurbulentTransport")

println("Loading packages on main process (this may take a minute)...")
using TurbulentTransport
using JSON
using Plots
using LaTeXStrings
using TJLF
using Interpolations
using Optim
using BlackBoxOptim
using Statistics
using DelimitedFiles
println("Main process packages loaded successfully!\n")

# Load packages on all workers
# Using @everywhere loads on all workers in parallel, but since we already
# precompiled on the main process, this should be fast and not cause conflicts
if nworkers() > 0
    println("Loading packages on $(nworkers()) workers...")
    @everywhere begin
        using Pkg
        Pkg.activate("/global/homes/t/trifolio/.julia/dev/TurbulentTransport")
    end

    @everywhere begin
        using TurbulentTransport
        using JSON
        using Plots
        using LaTeXStrings
        using TJLF
        using Interpolations
        using Optim
        using BlackBoxOptim
        using Statistics
        using DelimitedFiles
    end
    println("All workers ready!\n")
end

# Load data on all workers
@everywhere data = JSON.parsefile("qlnn_training_subset_merged_result_dict_dmn36p5_dict.json")

# Template input.tglf content (defined on all workers)
@everywhere input_tglf_lines = """
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
ZS_1 = -1.0
ZS_2 = 1.0
ZS_3 = 6.0
""";

@everywhere function setup_input_file(shot::Int, filepath::String)
    """Setup input.tglf file for a specific shot"""

    input_lines = split(input_tglf_lines, '\n')

    # Extract parameter names from template
    param_names = String[]
    for line in input_lines
        m = match(r"^([A-Za-z0-9_]+)\s*=", line)
        if m !== nothing
            push!(param_names, m.captures[1])
        end
    end

    # Build update dictionary
    params_to_update = Dict{String,Any}()
    for pname in param_names
        if haskey(data, pname) && isa(data[pname], AbstractVector) && length(data[pname]) >= shot
            params_to_update[pname] = data[pname][shot]
        end
    end

    # Always set these parameters
    params_to_update["SAT_RULE"] = 3
    params_to_update["NBASIS_MAX"] = 6
    params_to_update["NKY"] = 9
    params_to_update["UNITS"] = "CGYRO"
    params_to_update["FIND_WIDTH"] = false

    # Update lines in the template
    for i in 1:length(input_lines)
        m = match(r"^([A-Za-z0-9_]+)\s*=\s*(.*)", input_lines[i])
        if m !== nothing
            pname = m.captures[1]
            template_val_str = strip(m.captures[2])
            if haskey(params_to_update, pname)
                val = params_to_update[pname]
                # Convert Julia Bool to Fortran-style .true./.false.
                if val === true
                    val = ".true."
                elseif val === false
                    val = ".false."
                else
                    # Try to match the type of the template value
                    if occursin(r"^\-?\d+\.\d*([eE][\+\-]?\d+)?$", template_val_str)
                        val = float(val)
                    elseif occursin(r"^\-?\d+$", template_val_str)
                        val = Int(round(val))
                    end
                end
                input_lines[i] = "$(pname) = $(val)"
            end
        end
    end

    new_input_tglf_lines = join(input_lines, "\n")

    open(filepath, "w") do f
        write(f, new_input_tglf_lines)
    end

    return filepath
end

function to_scalar(val)
    """Safely convert a value to scalar Float64, handling nested arrays"""
    result = val
    while result isa AbstractArray && length(result) >= 1
        result = result[1]
    end
    return Float64(result)
end

@everywhere function to_scalar(val)
    """Safely convert a value to scalar Float64, handling nested arrays"""
    result = val
    while result isa AbstractArray && length(result) >= 1
        result = result[1]
    end
    return Float64(result)
end

@everywhere function run_tjlf_with_width_spectrum(input_tjlf, width_spectrum::Vector{Float64})
    """Run TJLF with a custom WIDTH_SPECTRUM vector

    Returns:
        Tuple of (gamma, frequency) arrays for all ky points
    """

    # Create a copy to avoid modifying the original
    input = deepcopy(input_tjlf)

    # Set the WIDTH_SPECTRUM - this should be a vector of length equal to ky grid
    input.WIDTH_SPECTRUM = width_spectrum

    # Make sure we're not using the scalar WIDTH when WIDTH_SPECTRUM is provided
    # (The TJLF code should prioritize WIDTH_SPECTRUM over WIDTH)

    satParams = get_sat_params(input)
    input.KY_SPECTRUM .= get_ky_spectrum(input, satParams.grad_r0)
    outputHermite = gauss_hermite(input)
    fluxes, eigenvalue = tjlf_TM(input, satParams, outputHermite)

    # Return both gamma (growth rate) and frequency
    gamma = eigenvalue[1, :, 1]
    frequency = eigenvalue[1, :, 2]
    return (gamma, frequency)
end

@everywhere function optimize_single_ky(
    ky_index::Int,
    input_tjlf,
    ky_grid::Vector{Float64},
    cgyro_gamma_target::Float64,
    baseline_width::Float64,
    lambda::Float64=0.01;  # Regularization parameter
    _version::Int=3  # Version marker to force recompilation - INCREMENT THIS WHEN DEBUGGING
)
    """
    VERSION: 2024-10-22-v3 - Added version parameter to force worker recompilation

    Optimize WIDTH for a single ky point using GLOBAL optimization (Differential Evolution)

    This version uses BlackBoxOptim.jl with differential evolution to avoid getting
    stuck in local minima. The WIDTH vs gamma relationship in TJLF is highly non-monotonic
    with multiple local minima, so global optimization is essential.

    Args:
        ky_index: Index of the ky point to optimize
        input_tjlf: TJLF input structure
        ky_grid: Full ky grid
        cgyro_gamma_target: Target gamma value from CGYRO for this ky
        baseline_width: Baseline WIDTH value for initial guess and bounds
        lambda: Regularization strength (default: 0.01)
                Higher values → stay closer to baseline
                Lower values → better fit but more variation

    Returns:
        Optimal WIDTH value for this ky point, optimal error
    """

    println("*** FUNCTION ENTRY: optimize_single_ky v3 called for ky_index=$ky_index ***")
    flush(stdout)  # Force output immediately

    n_ky_total = length(ky_grid)

    # Track best solution ourselves to avoid BlackBoxOptim's vector-to-scalar issues
    best_width = Ref(baseline_width)
    best_error = Ref(Inf)

    # Objective function for this single ky point with regularization
    function objective_1d(width_value)
        # Create WIDTH_SPECTRUM with baseline values everywhere except this ky
        width_spectrum = fill(baseline_width, n_ky_total)
        width_spectrum[ky_index] = width_value

        # Run TJLF
        try
            gamma_tjlf, freq_tjlf = run_tjlf_with_width_spectrum(input_tjlf, width_spectrum)
            gamma_this_ky = gamma_tjlf[ky_index]

            # CRITICAL FIX: Penalize cases where TJLF returns gamma=0
            # TJLF clamps negative growth rates to 0, but this creates a flat region
            # in the objective function that confuses the optimizer
            if gamma_this_ky == 0.0 && cgyro_gamma_target > 0.0
                # Return a very large error that's always worse than baseline
                # This prevents the optimizer from choosing WIDTH values that over-stabilize
                return 1e8
            end

            # Calculate squared error for this ky point
            mse = (gamma_this_ky - cgyro_gamma_target)^2

            # Add regularization term: penalize deviation from baseline
            # Using normalized squared difference to make lambda scale-independent
            regularization = lambda * ((width_value - baseline_width) / baseline_width)^2

            total_error = mse + regularization

            # Track best solution ourselves
            if total_error < best_error[]
                best_error[] = total_error
                best_width[] = width_value
            end

            return total_error
        catch e
            return 1e10  # Return large error if TJLF fails
        end
    end

    # 1D optimization with bounds
    lower = 0.1 * baseline_width
    upper = 10.0 * baseline_width

    # WORKAROUND: Use grid search + local refinement instead of BlackBoxOptim
    # BlackBoxOptim 0.6.3 has internal vector-to-scalar conversion issues with Julia 1.11.4
    # Grid search provides global exploration, then refine locally
    println("DEBUG ky=$ky_index: Starting grid search optimization")
    flush(stdout)

    # Phase 1: Coarse grid search (global exploration)
    n_grid = 20
    grid_points = range(lower, upper, length=n_grid)
    for width_val in grid_points
        objective_1d(width_val)  # This updates best_width[] and best_error[]
    end

    println("DEBUG ky=$ky_index: Grid search complete, best so far: width=$(best_width[]), error=$(best_error[])")
    flush(stdout)

    # Phase 2: Fine grid search around best point
    if best_error[] < Inf
        # Search +/- 20% around best point
        refine_lower = max(lower, best_width[] * 0.8)
        refine_upper = min(upper, best_width[] * 1.2)
        refine_points = range(refine_lower, refine_upper, length=n_grid)
        for width_val in refine_points
            objective_1d(width_val)
        end
    end

    println("DEBUG ky=$ky_index: Refinement complete")
    flush(stdout)

    # Use our tracked best solution
    width_out = Float64(best_width[])
    error_out = Float64(best_error[])

    println("DEBUG ky=$ky_index: Returning width=$width_out, error=$error_out")
    flush(stdout)

    # Return as tuple
    return (width_out, error_out)
end

@everywhere function optimize_width_spectrum_per_ky(
    shot::Int;
    ky_filter_mode::String="less_than",
    ky_threshold_low::Union{Float64,Nothing}=nothing,
    ky_threshold_high::Union{Float64,Nothing}=nothing,
    lambda::Float64=0.01,  # Regularization strength
    use_parallel::Bool=true  # Whether to use parallel processing for ky points
)
    """
    Optimize WIDTH_SPECTRUM by optimizing each ky point independently in parallel

    This is much faster than joint optimization when ky modes are independent.

    Args:
        shot: Shot index
        ky_filter_mode: How to filter ky values
        ky_threshold_low: Lower ky threshold
        ky_threshold_high: Upper ky threshold
        lambda: Regularization strength (default: 0.01)
                0.0 = no regularization (original behavior)
                0.01 = light regularization (recommended)
                0.1 = moderate regularization
                1.0 = strong regularization (stays very close to baseline)
    """

    println("=" ^ 60)
    println("PER-KY PARALLEL OPTIMIZATION FOR SHOT $shot")
    println("Shot ID: $(data["id"][shot])")
    println("ky filter mode: $ky_filter_mode")
    println("Regularization lambda: $lambda")
    println("=" ^ 60)

    # Setup input file
    tmpdir = mktempdir()
    filepath = joinpath(tmpdir, "input.tglf")
    setup_input_file(shot, filepath)

    # Load and setup TJLF
    input_tglf = TurbulentTransport.load(InputTGLF(), filepath)
    TurbulentTransport.apply_presets!(input_tglf)
    TurbulentTransport.run_tjlf(input_tglf)

    input_tjlf = readInput(filepath)
    satParams = get_sat_params(input_tjlf)
    input_tjlf.KY_SPECTRUM .= get_ky_spectrum(input_tjlf, satParams.grad_r0)

    # Get CGYRO data
    ky_cgyro = Float64.(data["cgyro_growthrate_spectra"][shot][1])
    cgyro_gamma = Float64.(data["cgyro_growthrate_spectra"][shot][2])
    cgyro_freq = Float64.(data["cgyro_growthrate_spectra"][shot][3])

    # Check for insufficient CGYRO data
    n_cgyro_points = length(ky_cgyro)
    if n_cgyro_points < 2
        println("WARNING: Shot $shot has insufficient CGYRO data (only $n_cgyro_points points)")
        println("         Skipping optimization for this shot - need at least 2 points for interpolation")
        println("=" ^ 60)
        # Return baseline WIDTH spectrum (no optimization)
        baseline_width = input_tjlf.WIDTH
        width_spectrum_optimal = fill(baseline_width, length(input_tjlf.KY_SPECTRUM))
        result = (converged=false, iterations=0, minimum=Inf)
        return width_spectrum_optimal, Float64[], Float64[], result, input_tjlf, Bool[], Float64[], Float64[]
    end

    println("CGYRO data: $n_cgyro_points points covering ky ∈ [$(round(minimum(ky_cgyro), digits=3)), $(round(maximum(ky_cgyro), digits=3))]")

    # Setup optimization domain
    ky_grid = input_tjlf.KY_SPECTRUM

    # Start with CGYRO data range as base filter
    mask = (ky_grid .>= minimum(ky_cgyro)) .& (ky_grid .<= maximum(ky_cgyro))

    # Apply additional ky filtering based on mode
    if ky_filter_mode == "less_than"
        threshold = isnothing(ky_threshold_high) ? 0.8 : ky_threshold_high
        mask = mask .& (ky_grid .< threshold)
        println("ky filter: ky < $threshold")
    elseif ky_filter_mode == "greater_than"
        if isnothing(ky_threshold_low)
            error("ky_threshold_low must be specified for 'greater_than' mode")
        end
        mask = mask .& (ky_grid .> ky_threshold_low)
        println("ky filter: ky > $ky_threshold_low")
    elseif ky_filter_mode == "both"
        low = isnothing(ky_threshold_low) ? 0.0 : ky_threshold_low
        high = isnothing(ky_threshold_high) ? 0.8 : ky_threshold_high
        mask = mask .& (ky_grid .> low) .& (ky_grid .< high)
        println("ky filter: $low < ky < $high")
    elseif ky_filter_mode == "full"
        println("ky filter: full CGYRO range (no additional filtering)")
    else
        error("Unknown ky_filter_mode: $ky_filter_mode")
    end

    ky_grid_in = ky_grid[mask]
    n_ky_opt = length(ky_grid_in)
    opt_indices = findall(mask)

    # Check if we have any ky points to optimize after filtering
    if n_ky_opt == 0
        println("WARNING: No ky points to optimize after filtering for Shot $shot")
        println("         CGYRO range: [$(round(minimum(ky_cgyro), digits=3)), $(round(maximum(ky_cgyro), digits=3))]")
        println("         TJLF ky_grid range: [$(round(minimum(ky_grid), digits=3)), $(round(maximum(ky_grid), digits=3))]")
        println("         Skipping optimization for this shot")
        println("=" ^ 60)
        # Return baseline WIDTH spectrum (no optimization)
        baseline_width = input_tjlf.WIDTH
        width_spectrum_optimal = fill(baseline_width, length(ky_grid))
        result = (converged=false, iterations=0, minimum=Inf)
        return width_spectrum_optimal, Float64[], Float64[], result, input_tjlf, mask, Float64[], Float64[]
    end

    println("Total ky points: $(length(ky_grid))")
    println("Optimization ky points: $n_ky_opt")
    println("Optimization ky range: $(round(minimum(ky_grid_in), digits=3)) to $(round(maximum(ky_grid_in), digits=3))")

    # Interpolate CGYRO gamma and frequency onto TJLF grid
    itp_gamma = interpolate((ky_cgyro,), cgyro_gamma, Gridded(Linear()))
    cgyro_gamma_interp = itp_gamma.(ky_grid_in)

    itp_freq = interpolate((ky_cgyro,), cgyro_freq, Gridded(Linear()))
    cgyro_freq_interp = itp_freq.(ky_grid_in)

    # Get baseline WIDTH value
    baseline_width = input_tjlf.WIDTH
    println("Baseline WIDTH: $baseline_width")
    println("\nOptimizing each ky point independently in parallel...")

    # Prepare arguments for parallel optimization
    opt_args = [(idx, input_tjlf, ky_grid, cgyro_gamma_interp[i], baseline_width, lambda)
                for (i, idx) in enumerate(opt_indices)]

    # Optimize all ky points (parallel or sequential based on use_parallel)
    start_time = time()
    if use_parallel
        println("Using parallel processing for ky points...")
        results = pmap(args -> optimize_single_ky(args...), opt_args;
                       retry_delays=ExponentialBackOff(n=3, first_delay=1.0, max_delay=10.0),
                       on_error=ex -> (println("Warning: Worker failed with error: $ex"); nothing))

        # Check for failures
        if any(isnothing, results)
            error("Some ky optimizations failed. Check worker connections and reduce n_workers.")
        end
    else
        println("Using sequential processing for ky points...")
        results = map(args -> optimize_single_ky(args...), opt_args)
    end
    elapsed_time = time() - start_time

    # Extract optimal widths and errors
    # CRITICAL FIX: Manually unwrap arrays in case workers don't have to_scalar
    # This is a failsafe for distributed computing issues
    function safe_extract(val)
        v = val
        while v isa AbstractArray && length(v) >= 1
            v = v[1]
        end
        return Float64(v)
    end

    optimal_widths = [safe_extract(r[1]) for r in results]
    optimal_errors = [safe_extract(r[2]) for r in results]

    # Create full WIDTH_SPECTRUM
    width_spectrum_optimal = fill(baseline_width, length(ky_grid))
    width_spectrum_optimal[opt_indices] = optimal_widths

    # Calculate total error
    total_error = sum(optimal_errors)

    # Calculate baseline error for comparison
    baseline_width_spectrum = fill(baseline_width, length(ky_grid))
    gamma_baseline, freq_baseline = run_tjlf_with_width_spectrum(input_tjlf, baseline_width_spectrum)
    baseline_mse = mean((gamma_baseline[mask] .- cgyro_gamma_interp).^2)

    println("\n" * "=" ^ 60)
    println("PER-KY OPTIMIZATION RESULTS")
    println("=" ^ 60)
    println("Optimization time: $(round(elapsed_time, digits=2)) seconds")
    println("Optimal MSE: $(round(total_error/n_ky_opt, digits=6))")
    println("Baseline MSE: $(round(baseline_mse, digits=6))")
    println("Error reduction: $(round(100*(1 - (total_error/n_ky_opt)/baseline_mse), digits=1))%")

    println("\nOptimal WIDTH values by ky:")
    for i in 1:n_ky_opt
        println("  ky = $(round(ky_grid_in[i], digits=3)): WIDTH = $(round(optimal_widths[i], digits=4))")
    end

    # Create a simple result object for compatibility
    result = (converged=true, iterations=elapsed_time, minimum=total_error/n_ky_opt)

    return width_spectrum_optimal, optimal_widths, ky_grid_in, result, input_tjlf, mask, cgyro_gamma_interp, cgyro_freq_interp
end

@everywhere function optimize_width_spectrum(
    shot::Int;
    ky_filter_mode::String="less_than",  # Options: "less_than", "greater_than", "both", "full"
    ky_threshold_low::Union{Float64,Nothing}=nothing,
    ky_threshold_high::Union{Float64,Nothing}=nothing
)
    """
    Optimize WIDTH_SPECTRUM vector for a specific shot

    Args:
        shot: Shot index
        ky_filter_mode: How to filter ky values
            - "less_than": Only ky < ky_threshold_high (default: 0.8)
            - "greater_than": Only ky > ky_threshold_low
            - "both": ky_threshold_low < ky < ky_threshold_high
            - "full": Use all ky values in CGYRO range (no additional filtering)
        ky_threshold_low: Lower ky threshold (for "greater_than" and "both" modes)
        ky_threshold_high: Upper ky threshold (for "less_than" and "both" modes)
    """

    println("=" ^ 60)
    println("OPTIMIZING WIDTH_SPECTRUM FOR SHOT $shot")
    println("Shot ID: $(data["id"][shot])")
    println("ky filter mode: $ky_filter_mode")
    println("=" ^ 60)

    # Setup input file
    tmpdir = mktempdir()
    filepath = joinpath(tmpdir, "input.tglf")
    setup_input_file(shot, filepath)

    # Load and setup TJLF
    input_tglf = TurbulentTransport.load(InputTGLF(), filepath)
    TurbulentTransport.apply_presets!(input_tglf)
    TurbulentTransport.run_tjlf(input_tglf)

    input_tjlf = readInput(filepath)
    satParams = get_sat_params(input_tjlf)
    input_tjlf.KY_SPECTRUM .= get_ky_spectrum(input_tjlf, satParams.grad_r0)

    # Get CGYRO data
    ky_cgyro = Float64.(data["cgyro_growthrate_spectra"][shot][1])
    cgyro_gamma = Float64.(data["cgyro_growthrate_spectra"][shot][2])

    # Setup optimization domain
    ky_grid = input_tjlf.KY_SPECTRUM

    # Start with CGYRO data range as base filter
    mask = (ky_grid .>= minimum(ky_cgyro)) .& (ky_grid .<= maximum(ky_cgyro))

    # Apply additional ky filtering based on mode
    if ky_filter_mode == "less_than"
        threshold = isnothing(ky_threshold_high) ? 0.8 : ky_threshold_high
        mask = mask .& (ky_grid .< threshold)
        println("ky filter: ky < $threshold")
    elseif ky_filter_mode == "greater_than"
        if isnothing(ky_threshold_low)
            error("ky_threshold_low must be specified for 'greater_than' mode")
        end
        mask = mask .& (ky_grid .> ky_threshold_low)
        println("ky filter: ky > $ky_threshold_low")
    elseif ky_filter_mode == "both"
        low = isnothing(ky_threshold_low) ? 0.0 : ky_threshold_low
        high = isnothing(ky_threshold_high) ? 0.8 : ky_threshold_high
        mask = mask .& (ky_grid .> low) .& (ky_grid .< high)
        println("ky filter: $low < ky < $high")
    elseif ky_filter_mode == "full"
        println("ky filter: full CGYRO range (no additional filtering)")
    else
        error("Unknown ky_filter_mode: $ky_filter_mode. Use 'less_than', 'greater_than', 'both', or 'full'")
    end

    ky_grid_in = ky_grid[mask]
    n_ky_opt = length(ky_grid_in)

    println("Total ky points: $(length(ky_grid))")
    println("Optimization ky points: $n_ky_opt")
    println("Optimization ky range: $(round(minimum(ky_grid_in), digits=3)) to $(round(maximum(ky_grid_in), digits=3))")

    # Interpolate CGYRO gamma onto TJLF grid
    itp = interpolate((ky_cgyro,), cgyro_gamma, Gridded(Linear()))
    cgyro_gamma_interp = itp.(ky_grid_in)

    # Get baseline WIDTH value to use as initial guess
    baseline_width = input_tjlf.WIDTH
    println("Baseline WIDTH: $baseline_width")

    # Initialize WIDTH_SPECTRUM - we need a vector for ALL ky points, not just optimization points
    n_ky_total = length(ky_grid)

    # Initial guess: use baseline WIDTH for all ky points
    width_spectrum_init = fill(baseline_width, n_ky_total)

    # We'll only optimize the width values for the ky points in our optimization region
    # The rest will remain at baseline values

    # Optimization bounds for the width values we're optimizing
    lower_bounds = fill(0.1 * baseline_width, n_ky_opt)
    upper_bounds = fill(10.0 * baseline_width, n_ky_opt)

    # Initial guess for optimization variables (only the ky points we're optimizing)
    x0 = fill(baseline_width, n_ky_opt)

    function objective(width_opt_values)
        # Create full WIDTH_SPECTRUM vector
        width_spectrum_full = copy(width_spectrum_init)

        # Update only the optimization region
        opt_indices = findall(mask)
        width_spectrum_full[opt_indices] = width_opt_values

        # Run TJLF with this WIDTH_SPECTRUM
        try
            gamma_tjlf, freq_tjlf = run_tjlf_with_width_spectrum(input_tjlf, width_spectrum_full)
            gamma_tjlf_opt = gamma_tjlf[mask]

            # Calculate MSE
            mse = mean((gamma_tjlf_opt .- cgyro_gamma_interp).^2)
            return mse
        catch e
            println("Error in TJLF run: $e")
            return 1e10  # Return large error if TJLF fails
        end
    end

    # Run baseline to get initial error
    println("\nCalculating baseline error...")
    baseline_error = objective(x0)
    println("Baseline MSE: $(round(baseline_error, digits=6))")

    # Optimize
    println("\nStarting optimization...")
    result = optimize(
        objective,
        lower_bounds, upper_bounds, x0, Fminbox(),
        Optim.Options(iterations=50, show_trace=true)
    )

    optimal_widths = Optim.minimizer(result)
    optimal_error = Optim.minimum(result)

    println("\n" * "=" ^ 60)
    println("OPTIMIZATION RESULTS")
    println("=" ^ 60)
    println("Optimal MSE: $(round(optimal_error, digits=6))")
    println("Baseline MSE: $(round(baseline_error, digits=6))")
    println("Error reduction: $(round(100*(1 - optimal_error/baseline_error), digits=1))%")

    println("\nOptimal WIDTH_SPECTRUM values (optimization region only):")
    for i in 1:n_ky_opt
        println("ky = $(round(ky_grid_in[i], digits=3)): WIDTH = $(round(optimal_widths[i], digits=4))")
    end

    # Create full optimal WIDTH_SPECTRUM
    width_spectrum_optimal = copy(width_spectrum_init)
    opt_indices = findall(mask)
    width_spectrum_optimal[opt_indices] = optimal_widths

    return width_spectrum_optimal, optimal_widths, ky_grid_in, result, input_tjlf, mask, cgyro_gamma_interp
end

function save_combined_results(results_dict::Dict, output_dir::String="width_spec_opt_grid")
    """
    Save combined results from all shots to a single JSON file in nested dictionary format

    Args:
        results_dict: Dictionary mapping shot_id to optimization results tuple
        output_dir: Output directory path

    Output format matches qlnn_training_subset_merged_result_dict_dmn36p5_dict.json:
        - Top-level keys are field names
        - Each field is a vector indexed by shot
        - Most fields contain nested vectors (one per shot, containing values per ky point)
    """
    mkpath(output_dir)

    # Initialize nested dictionary structure
    combined_data = Dict(
        "id" => [],
        "ky_grid" => [],
        "width_baseline" => [],
        "width_optimized_spectrum" => [],
        "gamma_cgyro" => [],
        "gamma_baseline" => [],
        "gamma_optimized" => [],
        "freq_cgyro" => [],
        "freq_baseline" => [],
        "freq_optimized" => [],
        "error_gamma_baseline" => [],
        "error_gamma_optimized" => []
    )

    # Process shots in sorted order for consistent indexing
    for shot in sort(collect(keys(results_dict)))
        result_tuple = results_dict[shot]
        width_spectrum_optimal, optimal_widths, ky_grid_in, result, input_tjlf, mask, cgyro_gamma_interp, cgyro_freq_interp = result_tuple

        # Skip shots that were not optimized (insufficient data)
        if isempty(optimal_widths)
            println("  Skipping Shot $shot in output (insufficient CGYRO data)")
            continue
        end

        # Get CGYRO data
        ky_cgyro = Float64.(data["cgyro_growthrate_spectra"][shot][1])
        cgyro_gamma = Float64.(data["cgyro_growthrate_spectra"][shot][2])
        cgyro_freq = Float64.(data["cgyro_growthrate_spectra"][shot][3])

        # Skip if insufficient CGYRO data for interpolation
        if length(ky_cgyro) < 2
            println("  Skipping Shot $shot in output (insufficient CGYRO data for interpolation)")
            continue
        end

        # Run TJLF with baseline and optimized widths
        ky_grid = input_tjlf.KY_SPECTRUM
        baseline_width_spectrum = fill(input_tjlf.WIDTH, length(ky_grid))

        gamma_baseline, freq_baseline = run_tjlf_with_width_spectrum(input_tjlf, baseline_width_spectrum)
        gamma_optimized, freq_optimized = run_tjlf_with_width_spectrum(input_tjlf, width_spectrum_optimal)

        # Interpolate CGYRO data onto full TJLF grid
        itp_gamma = interpolate((ky_cgyro,), cgyro_gamma, Gridded(Linear()))
        itp_freq = interpolate((ky_cgyro,), cgyro_freq, Gridded(Linear()))

        # Build vectors for this shot (only include ky points within CGYRO range)
        ky_vec = Float64[]
        width_baseline_vec = Float64[]
        width_optimized_vec = Float64[]
        gamma_cgyro_vec = Float64[]
        gamma_baseline_vec = Float64[]
        gamma_optimized_vec = Float64[]
        freq_cgyro_vec = Float64[]
        freq_baseline_vec = Float64[]
        freq_optimized_vec = Float64[]
        error_gamma_baseline_vec = Float64[]
        error_gamma_optimized_vec = Float64[]

        for (idx, ky) in enumerate(ky_grid)
            # Only include ky points within CGYRO range
            if ky >= minimum(ky_cgyro) && ky <= maximum(ky_cgyro)
                cgyro_gamma_val = itp_gamma(ky)
                cgyro_freq_val = itp_freq(ky)

                push!(ky_vec, ky)
                push!(width_baseline_vec, input_tjlf.WIDTH)
                push!(width_optimized_vec, width_spectrum_optimal[idx])
                push!(gamma_cgyro_vec, cgyro_gamma_val)
                push!(gamma_baseline_vec, gamma_baseline[idx])
                push!(gamma_optimized_vec, gamma_optimized[idx])
                push!(freq_cgyro_vec, cgyro_freq_val)
                push!(freq_baseline_vec, freq_baseline[idx])
                push!(freq_optimized_vec, freq_optimized[idx])
                push!(error_gamma_baseline_vec, (gamma_baseline[idx] - cgyro_gamma_val)^2)
                push!(error_gamma_optimized_vec, (gamma_optimized[idx] - cgyro_gamma_val)^2)
            end
        end

        # Append shot data to combined structure
        push!(combined_data["id"], data["id"][shot])
        push!(combined_data["ky_grid"], ky_vec)
        push!(combined_data["width_baseline"], width_baseline_vec)
        push!(combined_data["width_optimized_spectrum"], width_optimized_vec)
        push!(combined_data["gamma_cgyro"], gamma_cgyro_vec)
        push!(combined_data["gamma_baseline"], gamma_baseline_vec)
        push!(combined_data["gamma_optimized"], gamma_optimized_vec)
        push!(combined_data["freq_cgyro"], freq_cgyro_vec)
        push!(combined_data["freq_baseline"], freq_baseline_vec)
        push!(combined_data["freq_optimized"], freq_optimized_vec)
        push!(combined_data["error_gamma_baseline"], error_gamma_baseline_vec)
        push!(combined_data["error_gamma_optimized"], error_gamma_optimized_vec)
    end

    # Save to JSON
    json_path = joinpath(output_dir, "combined_results.json")
    open(json_path, "w") do f
        JSON.print(f, combined_data, 4)
    end

    println("Saved combined results to: $json_path")
    println("  Format: Nested dictionary with $(length(combined_data["id"])) shots")
    println("  Structure matches qlnn_training_subset_merged_result_dict_dmn36p5_dict.json")
    return json_path
end

function save_optimization_results(shot::Int, width_spectrum_optimal, optimal_widths, ky_grid_in,
                                   input_tjlf, mask, cgyro_gamma_interp, result, output_prefix="width_optimization")
    """Save optimization results to JSON and CSV files in organized folder structure"""

    # Create organized folder structure: output_prefix/shot_XX/
    shot_dir = joinpath(output_prefix, "shot_$(shot)")
    mkpath(shot_dir)  # Create directory if it doesn't exist

    # Get CGYRO data
    ky_cgyro = Float64.(data["cgyro_growthrate_spectra"][shot][1])
    cgyro_gamma = Float64.(data["cgyro_growthrate_spectra"][shot][2])
    ky_grid = input_tjlf.KY_SPECTRUM

    # Run baseline and optimized TJLF
    baseline_width_spectrum = fill(input_tjlf.WIDTH, length(ky_grid))
    gamma_baseline, freq_baseline = run_tjlf_with_width_spectrum(input_tjlf, baseline_width_spectrum)
    gamma_optimized, freq_optimized = run_tjlf_with_width_spectrum(input_tjlf, width_spectrum_optimal)

    # Calculate errors
    baseline_mse = mean((gamma_baseline[mask] .- cgyro_gamma_interp).^2)
    optimal_mse = mean((gamma_optimized[mask] .- cgyro_gamma_interp).^2)
    error_reduction = 100 * (1 - optimal_mse / baseline_mse)

    # --- Save JSON file with all results ---
    json_filename = joinpath(shot_dir, "results.json")
    results_dict = Dict(
        "shot_index" => shot,
        "shot_id" => data["id"][shot],
        "optimization_info" => Dict(
            "baseline_width" => input_tjlf.WIDTH,
            "baseline_mse" => baseline_mse,
            "optimal_mse" => optimal_mse,
            "error_reduction_percent" => error_reduction,
            "converged" => result.converged,
            "iterations" => result.iterations
        ),
        "ky_grid_full" => collect(ky_grid),
        "width_spectrum_optimal" => collect(width_spectrum_optimal),
        "optimization_region" => Dict(
            "ky_grid_optimized" => collect(ky_grid_in),
            "optimal_widths" => collect(optimal_widths),
            "mask" => collect(mask)
        ),
        "cgyro_data" => Dict(
            "ky" => collect(ky_cgyro),
            "gamma" => collect(cgyro_gamma)
        ),
        "tjlf_results" => Dict(
            "ky_grid" => collect(ky_grid),
            "gamma_baseline" => collect(gamma_baseline),
            "gamma_optimized" => collect(gamma_optimized)
        )
    )

    open(json_filename, "w") do f
        JSON.print(f, results_dict, 4)
    end
    println("  Saved: $json_filename")

    # --- Save CSV file with WIDTH_SPECTRUM ---
    csv_filename = joinpath(shot_dir, "width_spectrum.csv")
    csv_data = hcat(ky_grid, width_spectrum_optimal, fill(input_tjlf.WIDTH, length(ky_grid)))
    csv_header = ["ky", "optimized_width", "baseline_width"]

    open(csv_filename, "w") do f
        writedlm(f, [csv_header], ',')
        writedlm(f, csv_data, ',')
    end
    println("  Saved: $csv_filename")

    # --- Save CSV file with growth rates comparison ---
    csv_gamma_filename = joinpath(shot_dir, "gamma_comparison.csv")
    csv_gamma_data = hcat(ky_grid, gamma_baseline, gamma_optimized)
    csv_gamma_header = ["ky", "gamma_baseline", "gamma_optimized"]

    open(csv_gamma_filename, "w") do f
        writedlm(f, [csv_gamma_header], ',')
        writedlm(f, csv_gamma_data, ',')
    end
    println("  Saved: $csv_gamma_filename")

    # --- Save CSV file with CGYRO data ---
    csv_cgyro_filename = joinpath(shot_dir, "cgyro_data.csv")
    csv_cgyro_data = hcat(ky_cgyro, cgyro_gamma)
    csv_cgyro_header = ["ky", "gamma"]

    open(csv_cgyro_filename, "w") do f
        writedlm(f, [csv_cgyro_header], ',')
        writedlm(f, csv_cgyro_data, ',')
    end
    println("  Saved: $csv_cgyro_filename")

    return shot_dir, json_filename, csv_filename
end

function plot_gamma_comparison(shot::Int, width_spectrum_optimal, ky_grid_in, input_tjlf)
    """Plot gamma (growth rate) comparison between CGYRO, baseline TJLF, and optimized TJLF"""

    # Get CGYRO data for plotting
    ky_cgyro = Float64.(data["cgyro_growthrate_spectra"][shot][1])
    cgyro_gamma = Float64.(data["cgyro_growthrate_spectra"][shot][2])
    ky_grid = input_tjlf.KY_SPECTRUM

    # Run baseline TJLF
    baseline_width_spectrum = fill(input_tjlf.WIDTH, length(ky_grid))
    gamma_baseline, freq_baseline = run_tjlf_with_width_spectrum(input_tjlf, baseline_width_spectrum)

    # Run optimized TJLF
    gamma_optimized, freq_optimized = run_tjlf_with_width_spectrum(input_tjlf, width_spectrum_optimal)

    # Create gamma comparison plot
    plt = plot(
        ky_cgyro, cgyro_gamma,
        label = "CGYRO data",
        marker = :circle,
        lw = 2,
        legend = :bottomright,
        fontfamily = "Computer Modern",
        xscale = :log10,
        yscale = :log10,
        ylims = (1e-4, Inf),
        xlabel = L"k_y",
        ylabel = L"\gamma",
        title = "Growth Rate Comparison - Shot $shot",
        size = (800, 600)
    )

    # Plot baseline TJLF
    plot!(plt, ky_grid, gamma_baseline,
        label = "TJLF baseline (WIDTH=$(round(input_tjlf.WIDTH, digits=2)))",
        linestyle = :dash,
        alpha = 0.7,
        lw = 2
    )

    # Plot optimized TJLF
    plot!(plt, ky_grid, gamma_optimized,
        label = "TJLF optimized WIDTH_SPECTRUM",
        marker = :diamond,
        lw = 2,
        color = :red
    )

    # Highlight optimization region
    if !isempty(ky_grid_in)
        vline!(plt, [minimum(ky_grid_in)], linestyle = :dot, color = :gray, alpha = 0.5, label = "")
        vline!(plt, [maximum(ky_grid_in)], linestyle = :dot, color = :gray, alpha = 0.5, label = "")
    end

    return plt
end

function plot_frequency_comparison(shot::Int, width_spectrum_optimal, ky_grid_in, input_tjlf)
    """Plot frequency comparison between CGYRO, baseline TJLF, and optimized TJLF

    Includes automatic sign-flip detection: if TJLF and CGYRO frequencies have opposite
    signs at high ky (last 2-3 points), TJLF frequencies are flipped for better agreement.
    """

    # Get CGYRO data for plotting
    ky_cgyro = Float64.(data["cgyro_growthrate_spectra"][shot][1])
    cgyro_freq = Float64.(data["cgyro_growthrate_spectra"][shot][3])
    ky_grid = input_tjlf.KY_SPECTRUM

    # Run baseline TJLF
    baseline_width_spectrum = fill(input_tjlf.WIDTH, length(ky_grid))
    gamma_baseline, freq_baseline = run_tjlf_with_width_spectrum(input_tjlf, baseline_width_spectrum)

    # Run optimized TJLF
    gamma_optimized, freq_optimized = run_tjlf_with_width_spectrum(input_tjlf, width_spectrum_optimal)

    # Interpolate CGYRO frequencies onto TJLF grid for sign comparison
    itp_freq = interpolate((ky_cgyro,), cgyro_freq, Gridded(Linear()))

    # Check sign agreement at high ky (last 2-3 points in the ky_grid)
    # Only check points within CGYRO range
    n_check = min(3, length(ky_grid))
    high_ky_indices = []
    for i in (length(ky_grid) - n_check + 1):length(ky_grid)
        if ky_grid[i] >= minimum(ky_cgyro) && ky_grid[i] <= maximum(ky_cgyro)
            push!(high_ky_indices, i)
        end
    end

    sign_flip_needed = false
    if !isempty(high_ky_indices) && length(high_ky_indices) >= 2
        # Compare signs at high ky points
        cgyro_signs = [sign(itp_freq(ky_grid[i])) for i in high_ky_indices]
        tjlf_baseline_signs = [sign(freq_baseline[i]) for i in high_ky_indices]
        tjlf_opt_signs = [sign(freq_optimized[i]) for i in high_ky_indices]

        # If majority of TJLF signs disagree with CGYRO, flip the sign
        baseline_disagree = sum(cgyro_signs .!= tjlf_baseline_signs) > length(high_ky_indices) / 2
        opt_disagree = sum(cgyro_signs .!= tjlf_opt_signs) > length(high_ky_indices) / 2

        if baseline_disagree || opt_disagree
            sign_flip_needed = true
            freq_baseline = -freq_baseline
            freq_optimized = -freq_optimized
            println("  Note: Sign flip applied to TJLF frequencies for Shot $shot (high-ky sign mismatch detected)")
        end
    end

    # Create frequency comparison plot
    # Note: Frequency can be positive or negative (direction of rotation)
    # So we use log scale only for x-axis (ky), not y-axis
    title_suffix = sign_flip_needed ? " (TJLF sign flipped)" : ""
    plt = plot(
        ky_cgyro, cgyro_freq,
        label = "CGYRO data",
        marker = :circle,
        lw = 2,
        legend = :best,
        fontfamily = "Computer Modern",
        xscale = :log10,
        xlabel = L"k_y",
        ylabel = L"\omega",
        title = "Frequency Comparison - Shot $shot" * title_suffix,
        size = (800, 600)
    )

    # Plot baseline TJLF
    plot!(plt, ky_grid, freq_baseline,
        label = "TJLF baseline (WIDTH=$(round(input_tjlf.WIDTH, digits=2)))",
        linestyle = :dash,
        alpha = 0.7,
        lw = 2
    )

    # Plot optimized TJLF
    plot!(plt, ky_grid, freq_optimized,
        label = "TJLF optimized WIDTH_SPECTRUM",
        marker = :diamond,
        lw = 2,
        color = :red
    )

    # Add horizontal line at zero for reference
    hline!(plt, [0.0], linestyle = :dot, color = :black, alpha = 0.3, label = "")

    # Highlight optimization region
    if !isempty(ky_grid_in)
        vline!(plt, [minimum(ky_grid_in)], linestyle = :dot, color = :gray, alpha = 0.5, label = "")
        vline!(plt, [maximum(ky_grid_in)], linestyle = :dot, color = :gray, alpha = 0.5, label = "")
    end

    return plt
end

# Main execution (single shot)
function main_fn(;
    shot::Int=11,
    ky_filter_mode::String="less_than",
    ky_threshold_low::Union{Float64,Nothing}=nothing,
    ky_threshold_high::Union{Float64,Nothing}=nothing,
    parallel_mode::String="joint"  # "joint" or "per_ky"
)
    """
    Run optimization for a single shot

    Args:
        shot: Shot index (default: 11)
        ky_filter_mode: How to filter ky values (default: "less_than")
        ky_threshold_low: Lower ky threshold
        ky_threshold_high: Upper ky threshold
        parallel_mode: "joint" (optimize all ky jointly) or "per_ky" (optimize each ky independently in parallel)
    """

    # Choose optimization function based on parallel_mode
    if parallel_mode == "per_ky"
        width_spectrum_optimal, optimal_widths, ky_grid_in, result, input_tjlf, mask, cgyro_gamma_interp =
            optimize_width_spectrum_per_ky(shot;
                ky_filter_mode=ky_filter_mode,
                ky_threshold_low=ky_threshold_low,
                ky_threshold_high=ky_threshold_high
            )
    elseif parallel_mode == "joint"
        width_spectrum_optimal, optimal_widths, ky_grid_in, result, input_tjlf, mask, cgyro_gamma_interp =
            optimize_width_spectrum(shot;
                ky_filter_mode=ky_filter_mode,
                ky_threshold_low=ky_threshold_low,
                ky_threshold_high=ky_threshold_high
            )
    else
        error("Unknown parallel_mode: $parallel_mode. Use 'joint' or 'per_ky'")
    end

    # Save data files
    println("\nSaving results to data files...")
    shot_dir, _, _ = save_optimization_results(shot, width_spectrum_optimal, optimal_widths, ky_grid_in,
                                                input_tjlf, mask, cgyro_gamma_interp, result)

    # Plot results
    plt = plot_width_spectrum_results(shot, width_spectrum_optimal, optimal_widths, ky_grid_in, input_tjlf, mask, cgyro_gamma_interp)
    display(plt)

    # Save plot to shot directory
    plot_filename = joinpath(shot_dir, "optimization_plot.png")
    savefig(plot_filename)
    println("  Saved: $plot_filename")

    return width_spectrum_optimal, result
end

# Main execution (parallel - multiple shots)
function main_parallel(
    shots::Vector{Int};
    ky_filter_mode::String="less_than",
    ky_threshold_low::Union{Float64,Nothing}=nothing,
    ky_threshold_high::Union{Float64,Nothing}=nothing,
    parallel_mode::String="per_shot",  # "per_shot", "per_ky", or "both"
    lambda::Float64=0.01  # Regularization for per_ky mode
)
    """
    Run optimization for multiple shots in parallel

    Args:
        shots: Vector of shot indices to optimize
        ky_filter_mode: How to filter ky values (default: "less_than")
        ky_threshold_low: Lower ky threshold
        ky_threshold_high: Upper ky threshold
        parallel_mode: Parallelization strategy
            - "per_shot": Parallelize across shots using joint optimization (default, original behavior)
            - "per_ky": Use per-ky optimization for a single shot (ignores multiple shots, uses first shot only)
            - "both": Parallelize across shots AND within each shot (per-ky) - requires many workers

    Returns:
        Dictionary mapping shot index to optimization results
    """
    println("=" ^ 80)
    println("PARALLEL WIDTH_SPECTRUM OPTIMIZATION")
    println("Number of workers: $(nworkers())")
    println("Shots to optimize: $shots")
    println("ky filter mode: $ky_filter_mode")
    println("Parallelization mode: $parallel_mode")
    println("=" ^ 80)

    # Select optimization function based on parallel_mode
    if parallel_mode == "per_shot"
        optimize_func = shot -> optimize_width_spectrum(shot;
            ky_filter_mode=ky_filter_mode,
            ky_threshold_low=ky_threshold_low,
            ky_threshold_high=ky_threshold_high
        )
    elseif parallel_mode == "per_ky"
        if length(shots) > 1
            println("WARNING: per_ky mode with multiple shots - only processing first shot: $(shots[1])")
            println("For multiple shots with per-ky optimization, use parallel_mode='both'")
            shots = [shots[1]]
        end
        optimize_func = shot -> optimize_width_spectrum_per_ky(shot;
            ky_filter_mode=ky_filter_mode,
            ky_threshold_low=ky_threshold_low,
            ky_threshold_high=ky_threshold_high,
            lambda=lambda
        )
    elseif parallel_mode == "both"
        # This uses per-ky optimization for each shot, and parallelizes across shots
        # FULL NESTED PARALLELISM: Parallelizes both across shots AND within each shot
        println("Using FULL nested parallelization: per-shot (parallel) + per-ky (parallel)")
        println("Note: Requires many workers (~80 for 4 shots × 20 ky points)")
        optimize_func = shot -> optimize_width_spectrum_per_ky(shot;
            ky_filter_mode=ky_filter_mode,
            ky_threshold_low=ky_threshold_low,
            ky_threshold_high=ky_threshold_high,
            lambda=lambda,
            use_parallel=true  # Enable full nested parallelism
        )
    else
        error("Unknown parallel_mode: $parallel_mode. Use 'per_shot', 'per_ky', or 'both'")
    end

    # Run optimizations in parallel using pmap
    println("\nStarting parallel optimizations...")
    start_time = time()

    # Use pmap with retry and error handling for robustness
    results = pmap(optimize_func, shots;
                   retry_delays=ExponentialBackOff(n=3, first_delay=1.0, max_delay=10.0),
                   on_error=ex -> (println("Warning: Shot optimization failed with error: $ex"); nothing))

    # Check for failures
    if any(isnothing, results)
        error("Some shot optimizations failed. Check worker connections and reduce n_workers.")
    end

    elapsed_time = time() - start_time
    println("\n" * "=" ^ 80)
    println("ALL OPTIMIZATIONS COMPLETE")
    println("Total time: $(round(elapsed_time, digits=2)) seconds")
    println("Average time per shot: $(round(elapsed_time/length(shots), digits=2)) seconds")
    println("=" ^ 80)

    # Create dictionary for easy access
    results_dict = Dict{Int, Any}()
    for (i, shot) in enumerate(shots)
        results_dict[shot] = results[i]
    end

    # Save combined results to single JSON file
    println("\nSaving combined results...")
    output_dir = "width_spec_opt_grid"
    save_combined_results(results_dict, output_dir)

    # Generate plots for all shots
    println("\nGenerating plots for all shots...")
    shot_figs_dir = joinpath(output_dir, "shot_figs")
    mkpath(shot_figs_dir)

    for shot in shots
        width_spectrum_optimal, optimal_widths, ky_grid_in, result, input_tjlf, mask, cgyro_gamma_interp, cgyro_freq_interp = results_dict[shot]

        # Skip plotting for shots with insufficient data
        if isempty(optimal_widths)
            println("\nShot $shot: Skipped (insufficient CGYRO data)")
            continue
        end

        println("\nShot $shot:")

        # Generate and save gamma comparison plot
        plt_gamma = plot_gamma_comparison(shot, width_spectrum_optimal, ky_grid_in, input_tjlf)
        gamma_filename = joinpath(shot_figs_dir, "shot_$(shot)_gamma.png")
        savefig(plt_gamma, gamma_filename)
        println("  Saved: $gamma_filename")

        # Generate and save frequency comparison plot
        plt_freq = plot_frequency_comparison(shot, width_spectrum_optimal, ky_grid_in, input_tjlf)
        freq_filename = joinpath(shot_figs_dir, "shot_$(shot)_frequency.png")
        savefig(plt_freq, freq_filename)
        println("  Saved: $freq_filename")
    end

    println("\n" * "=" ^ 80)
    println("All results saved to: $output_dir")
    println("  - Combined data: combined_results.json")
    println("  - Figures: shot_figs/")
    println("=" ^ 80)

    return results_dict
end

println("Functions defined successfully!")
println("Worker count: $(nworkers())")
println("SCRIPT VERSION: 2024-10-22-v4-grid (restructured output with frequency data)")
println("=" ^ 60)

# Verify that workers have the updated function
println("\nVerifying workers have loaded the updated code...")
test_result = @fetchfrom 2 begin
    hasmethod(optimize_single_ky, Tuple{Int, Any, Vector{Float64}, Float64, Float64, Float64}) ||
    hasmethod(optimize_single_ky, (Int, Any, Vector{Float64}, Float64, Float64, Float64))
end
println("Worker 2 has optimize_single_ky function: $test_result")
println("If false, workers haven't loaded the new code - restart Julia!")
println("=" ^ 60)

# ============================================================================
# EXECUTION EXAMPLES
# ============================================================================

shots_to_optimize = collect(1:20)
# Alternative options:
# shots_to_optimize = [1,15,8,10,11,12]  # Specific shots
# shots_to_optimize = 1:20  # Range (also works, but collect() is clearer)
# shots_to_optimize = collect(1:10)  # First 10 shots
# shots_to_optimize = collect(11:20)  # Last 10 shots

width_spectrum_results = main_parallel(shots_to_optimize,
    parallel_mode="both",
    ky_filter_mode="full",
    lambda=0.0
)
