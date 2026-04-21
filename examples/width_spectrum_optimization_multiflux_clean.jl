# ============================================================================
# WIDTH_SPECTRUM Optimization with 4-Component QL Weights Cost Function
# VERSION: 2025-11-26-v4-normalized
#
# Per-ky optimization using raw QL_weights (not flux_out with intensities)
# Proper normalization: loss = asinh((TJLF - CGYRO) / STD)^2
# Auto-computes normalization constants from dataset at startup
#
# IMPORTANT: If you edit @everywhere functions, restart Julia to reload workers
# ============================================================================

# ============================================================================
# COST FUNCTION WEIGHT CONFIGURATION
# ============================================================================
# Set gamma weight (0.0 to 1.0). Remaining weight distributed to QL weights
# Example: GAMMA_WEIGHT = 0.25 means 25% gamma, 75% QL weights (25% each: Γ_e, Q_e, Q_ions)
const GAMMA_WEIGHT = 0.95
const QL_WEIGHT = 1.0 - GAMMA_WEIGHT
const QL_COMPONENT_WEIGHT = QL_WEIGHT / 3.0  # Divided among 3 QL components
# ============================================================================

using Distributed

# Add worker processes
if nprocs() == 1
    n_workers = 120
    println("Adding $n_workers worker processes...")
    addprocs(n_workers,
             exeflags="--project=/global/homes/t/trifolio/.julia/dev/TurbulentTransport",
             topology=:all_to_all)
    println("Workers added\n")
end

# Load packages
using Pkg
Pkg.activate("/global/homes/t/trifolio/.julia/dev/TurbulentTransport")

println("Loading packages...")
using TurbulentTransport, JSON, Plots, LaTeXStrings, TJLF
using Interpolations, Optim, BlackBoxOptim, Statistics, DelimitedFiles
println("Main process ready!\n")

# Load packages on workers
if nworkers() > 0
    println("Loading packages on $(nworkers()) workers...")
    @everywhere begin
        # Note: Project environment already activated via exeflags in addprocs()
        # No need to call Pkg.activate() here - it causes pidfile race conditions
        using TurbulentTransport, JSON, Plots, LaTeXStrings, TJLF
        using Interpolations, Optim, BlackBoxOptim, Statistics, DelimitedFiles
    end
    println("Workers ready!\n")

    # Share weight configuration with workers
    @eval @everywhere const GAMMA_WEIGHT = $GAMMA_WEIGHT
    @eval @everywhere const QL_WEIGHT = $QL_WEIGHT
    @eval @everywhere const QL_COMPONENT_WEIGHT = $QL_COMPONENT_WEIGHT
end

# Load data
@everywhere data = JSON.parsefile("qlnn_training_subset_merged_result_dict_dmn36p5_dict.json")
# @everywhere data = JSON.parsefile("qlnn_training_subset_merged_result_dict_dmn40p2_dict_test_500.json")

# ============================================================================
# COMPUTE NORMALIZATION STATISTICS FROM DATASET
# ============================================================================
println("Computing normalization statistics from dataset...")
all_gamma = Float64[]
all_gamma_electron = Float64[]
all_q_electron = Float64[]
all_q_ions = Float64[]

n_shots = length(data["id"])
for shot in 1:n_shots
    # Get CGYRO gamma
    cgyro_gamma = Float64.(data["cgyro_growthrate_spectra"][shot][2])
    append!(all_gamma, cgyro_gamma)

    # Get CGYRO fluxes
    cgyro_fluxes = data["cgyro_ql_fluxes"][shot]
    n_ky = length(cgyro_fluxes)

    for ky_idx in 1:n_ky
        # Extract QL weights: species [0=D, 1=C, 2=e], field [0,1,2], type [0=particle, 1=energy]
        gamma_e = sum(cgyro_fluxes[ky_idx][3][field+1][1] for field in 0:2)
        q_e = sum(cgyro_fluxes[ky_idx][3][field+1][2] for field in 0:2)
        q_ions = sum(cgyro_fluxes[ky_idx][species+1][field+1][2]
                     for species in 0:1, field in 0:2)

        push!(all_gamma_electron, gamma_e)
        push!(all_q_electron, q_e)
        push!(all_q_ions, q_ions)
    end
end

GAMMA_STD = std(all_gamma)
GAMMA_ELECTRON_STD = std(all_gamma_electron)
Q_ELECTRON_STD = std(all_q_electron)
Q_IONS_STD = std(all_q_ions)

println("Normalization constants (std):")
println("  GAMMA_STD = $(round(GAMMA_STD, digits=6))")
println("  GAMMA_ELECTRON_STD = $(round(GAMMA_ELECTRON_STD, digits=6))")
println("  Q_ELECTRON_STD = $(round(Q_ELECTRON_STD, digits=6))")
println("  Q_IONS_STD = $(round(Q_IONS_STD, digits=6))")
println()

# Share with workers (use global vars, not const, to avoid redeclaration errors)
@everywhere GAMMA_STD = $GAMMA_STD
@everywhere GAMMA_ELECTRON_STD = $GAMMA_ELECTRON_STD
@everywhere Q_ELECTRON_STD = $Q_ELECTRON_STD
@everywhere Q_IONS_STD = $Q_IONS_STD
# ============================================================================

# Template input.tglf file
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
C_B = 0.315
SIG_B = 0.34
BOUNCE_COEFF = 3.0
C_NORM = 1.82770384
C_EXP = 1.39786897
C_COEFF = 0.36017009
C_ETG = 1.25
""";

@everywhere function setup_input_file(shot::Int, filepath::String)
    """Setup input.tglf file for a specific shot"""
    input_lines = split(input_tglf_lines, '\n')

    param_names = String[]
    for line in input_lines
        m = match(r"^([A-Za-z0-9_]+)\s*=", line)
        m !== nothing && push!(param_names, m.captures[1])
    end

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

    for i in 1:length(input_lines)
        m = match(r"^([A-Za-z0-9_]+)\s*=\s*(.*)", input_lines[i])
        if m !== nothing
            pname = m.captures[1]
            template_val_str = strip(m.captures[2])
            if haskey(params_to_update, pname)
                val = params_to_update[pname]
                if val === true
                    val = ".true."
                elseif val === false
                    val = ".false."
                else
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

    open(filepath, "w") do f
        write(f, join(input_lines, "\n"))
    end

    return filepath
end

@everywhere function run_tjlf_with_width_spectrum(input_tjlf, width_spectrum::Vector{Float64})
    """Run TJLF with custom WIDTH_SPECTRUM, return (gamma, frequency, QL_weights)"""
    input = deepcopy(input_tjlf)
    input.WIDTH_SPECTRUM = width_spectrum

    satParams = get_sat_params(input)
    input.KY_SPECTRUM .= get_ky_spectrum(input, satParams.grad_r0)
    outputHermite = gauss_hermite(input)
    QL_weights, eigenvalue = tjlf_TM(input, satParams, outputHermite)

    gamma = eigenvalue[1, :, 1]
    frequency = eigenvalue[1, :, 2]
    return (gamma, frequency, QL_weights)
end

@everywhere function optimize_single_ky(
    ky_index::Int,
    input_tjlf,
    ky_grid::Vector{Float64},
    cgyro_gamma_target::Float64,
    cgyro_gamma_electron_target::Float64,
    cgyro_q_electron_target::Float64,
    cgyro_q_ions_target::Float64,
    baseline_width::Float64,
    lambda::Float64=0.01;
    _version::Int=8
)
    """Optimize WIDTH using configurable gamma/QL weight split with proper normalization"""

    n_ky_total = length(ky_grid)
    best_width = Ref(baseline_width)
    best_error = Ref(Inf)

    function objective_1d(width_value)
        width_spectrum = fill(baseline_width, n_ky_total)
        width_spectrum[ky_index] = width_value

        try
            gamma_tjlf, freq_tjlf, QL_weights = run_tjlf_with_width_spectrum(input_tjlf, width_spectrum)
            gamma_this_ky = gamma_tjlf[ky_index]

            # Penalize TJLF gamma=0 cases (clamped negative values)
            if gamma_this_ky == 0.0 && cgyro_gamma_target > 0.0
                return 1e8
            end

            # Extract per-ky QL weights: QL_weights[field, species, mode, ky, type]
            # species: 1=electron, 2=D, 3=Carbon; type: 1=particle, 2=energy
            gamma_electron_ky = sum(QL_weights[:, 1, 1, ky_index, 1])
            q_electron_ky = sum(QL_weights[:, 1, 1, ky_index, 2])
            q_ions_ky = sum(QL_weights[:, 2:3, 1, ky_index, 2])

            # Calculate normalized errors: asinh(diff / STD)
            # This ensures all components contribute equally regardless of scale
            gamma_diff = gamma_this_ky - cgyro_gamma_target
            gamma_error = asinh(gamma_diff / max(GAMMA_STD, 1e-10))

            gamma_electron_diff = gamma_electron_ky - cgyro_gamma_electron_target
            gamma_electron_error = asinh(gamma_electron_diff / max(GAMMA_ELECTRON_STD, 1e-10))

            q_electron_diff = q_electron_ky - cgyro_q_electron_target
            q_electron_error = asinh(q_electron_diff / max(Q_ELECTRON_STD, 1e-10))

            q_ions_diff = q_ions_ky - cgyro_q_ions_target
            q_ions_error = asinh(q_ions_diff / max(Q_IONS_STD, 1e-10))

            # 4-component loss using configurable weights:
            component_loss = GAMMA_WEIGHT * gamma_error^2 +
                           QL_COMPONENT_WEIGHT * (gamma_electron_error^2 + q_electron_error^2 + q_ions_error^2)

            regularization = lambda * ((width_value - baseline_width) / baseline_width)^2
            total_error = component_loss + regularization

            if total_error < best_error[]
                best_error[] = total_error
                best_width[] = width_value
            end

            return total_error
        catch e
            println("ERROR in TJLF run for ky_index=$ky_index: $e")
            return 1e10
        end
    end

    # Grid search optimization (2 phases)
    lower = 0.1 * baseline_width
    upper = 10.0 * baseline_width

    # Phase 1: Coarse grid
    n_grid = 20
    for width_val in range(lower, upper, length=n_grid)
        objective_1d(width_val)
    end

    # Phase 2: Fine grid around best point
    if best_error[] < Inf
        refine_lower = max(lower, best_width[] * 0.8)
        refine_upper = min(upper, best_width[] * 1.2)
        for width_val in range(refine_lower, refine_upper, length=n_grid)
            objective_1d(width_val)
        end
    end

    return (Float64(best_width[]), Float64(best_error[]))
end

@everywhere function optimize_width_spectrum_per_ky(
    shot::Int;
    ky_filter_mode::String="less_than",
    ky_threshold_low::Union{Float64,Nothing}=nothing,
    ky_threshold_high::Union{Float64,Nothing}=nothing,
    lambda::Float64=0.01,
    use_parallel::Bool=true
)
    """Optimize WIDTH_SPECTRUM per-ky using 4-component cost function"""

    println("=" ^ 60)
    println("PER-KY OPTIMIZATION FOR SHOT $shot (4-COMPONENT MULTIFLUX)")
    println("Shot ID: $(data["id"][shot])")
    println("=" ^ 60)

    # Setup TJLF
    tmpdir = mktempdir()
    filepath = joinpath(tmpdir, "input.tglf")
    setup_input_file(shot, filepath)

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

    if length(ky_cgyro) < 2
        println("WARNING: Insufficient CGYRO data, skipping shot")
        baseline_width = input_tjlf.WIDTH
        width_spectrum_optimal = fill(baseline_width, length(input_tjlf.KY_SPECTRUM))
        result = (converged=false, iterations=0, minimum=Inf)
        return width_spectrum_optimal, Float64[], Float64[], result, input_tjlf, Bool[], Float64[], Float64[]
    end

    # Extract CGYRO fluxes: [shot][ky][species][field][type]
    cgyro_fluxes = data["cgyro_ql_fluxes"][shot]
    n_ky_cgyro = length(ky_cgyro)
    cgyro_gamma_electron = zeros(Float64, n_ky_cgyro)
    cgyro_q_electron = zeros(Float64, n_ky_cgyro)
    cgyro_q_ions = zeros(Float64, n_ky_cgyro)

    for ky_idx in 1:n_ky_cgyro
        # CGYRO indexing: species [0=D, 1=C, 2=e], field [0,1,2], type [0=particle, 1=energy]
        cgyro_gamma_electron[ky_idx] = sum(cgyro_fluxes[ky_idx][3][field+1][1] for field in 0:2)
        cgyro_q_electron[ky_idx] = sum(cgyro_fluxes[ky_idx][3][field+1][2] for field in 0:2)
        cgyro_q_ions[ky_idx] = sum(cgyro_fluxes[ky_idx][species+1][field+1][2]
                                   for species in 0:1, field in 0:2)
    end

    # Setup optimization domain
    ky_grid = input_tjlf.KY_SPECTRUM
    mask = (ky_grid .>= minimum(ky_cgyro)) .& (ky_grid .<= maximum(ky_cgyro))

    # Apply ky filtering
    if ky_filter_mode == "less_than"
        threshold = isnothing(ky_threshold_high) ? 0.8 : ky_threshold_high
        mask = mask .& (ky_grid .< threshold)
    elseif ky_filter_mode == "greater_than"
        isnothing(ky_threshold_low) && error("ky_threshold_low required for 'greater_than'")
        mask = mask .& (ky_grid .> ky_threshold_low)
    elseif ky_filter_mode == "both"
        low = isnothing(ky_threshold_low) ? 0.0 : ky_threshold_low
        high = isnothing(ky_threshold_high) ? 0.8 : ky_threshold_high
        mask = mask .& (ky_grid .> low) .& (ky_grid .< high)
    elseif ky_filter_mode != "full"
        error("Unknown ky_filter_mode: $ky_filter_mode")
    end

    ky_grid_in = ky_grid[mask]
    opt_indices = findall(mask)
    n_ky_opt = length(ky_grid_in)

    if n_ky_opt == 0
        println("WARNING: No ky points to optimize after filtering")
        baseline_width = input_tjlf.WIDTH
        width_spectrum_optimal = fill(baseline_width, length(ky_grid))
        result = (converged=false, iterations=0, minimum=Inf)
        return width_spectrum_optimal, Float64[], Float64[], result, input_tjlf, mask, Float64[], Float64[]
    end

    println("Optimizing $n_ky_opt ky points in range [$(round(minimum(ky_grid_in), digits=3)), $(round(maximum(ky_grid_in), digits=3))]")

    # Interpolate CGYRO data onto TJLF grid
    itp_gamma = interpolate((ky_cgyro,), cgyro_gamma, Gridded(Linear()))
    itp_freq = interpolate((ky_cgyro,), cgyro_freq, Gridded(Linear()))
    itp_gamma_electron = interpolate((ky_cgyro,), cgyro_gamma_electron, Gridded(Linear()))
    itp_q_electron = interpolate((ky_cgyro,), cgyro_q_electron, Gridded(Linear()))
    itp_q_ions = interpolate((ky_cgyro,), cgyro_q_ions, Gridded(Linear()))

    cgyro_gamma_interp = itp_gamma.(ky_grid_in)
    cgyro_freq_interp = itp_freq.(ky_grid_in)
    cgyro_gamma_electron_interp = itp_gamma_electron.(ky_grid_in)
    cgyro_q_electron_interp = itp_q_electron.(ky_grid_in)
    cgyro_q_ions_interp = itp_q_ions.(ky_grid_in)

    baseline_width = input_tjlf.WIDTH
    println("Baseline WIDTH: $baseline_width\n")

    # Prepare optimization arguments
    opt_args = [(idx, input_tjlf, ky_grid,
                 cgyro_gamma_interp[i],
                 cgyro_gamma_electron_interp[i],
                 cgyro_q_electron_interp[i],
                 cgyro_q_ions_interp[i],
                 baseline_width, lambda)
                for (i, idx) in enumerate(opt_indices)]

    # Run optimization
    start_time = time()
    if use_parallel
        results = pmap(args -> optimize_single_ky(args...), opt_args;
                       retry_delays=ExponentialBackOff(n=3, first_delay=1.0, max_delay=10.0),
                       on_error=ex -> (println("Warning: Worker failed: $ex"); nothing))
        any(isnothing, results) && error("Some ky optimizations failed")
    else
        results = map(args -> optimize_single_ky(args...), opt_args)
    end
    elapsed_time = time() - start_time

    # Extract results
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
    total_error = sum(optimal_errors)

    # Calculate baseline error
    baseline_width_spectrum = fill(baseline_width, length(ky_grid))
    gamma_baseline, freq_baseline, QL_weights_baseline = run_tjlf_with_width_spectrum(input_tjlf, baseline_width_spectrum)

    baseline_loss = 0.0
    for (i, idx) in enumerate(opt_indices)
        # Use same normalization as optimization
        gamma_diff = gamma_baseline[idx] - cgyro_gamma_interp[i]
        gamma_error = asinh(gamma_diff / max(GAMMA_STD, 1e-10))

        gamma_electron_baseline = sum(QL_weights_baseline[:, 1, 1, idx, 1])
        gamma_electron_diff = gamma_electron_baseline - cgyro_gamma_electron_interp[i]
        gamma_electron_error = asinh(gamma_electron_diff / max(GAMMA_ELECTRON_STD, 1e-10))

        q_electron_baseline = sum(QL_weights_baseline[:, 1, 1, idx, 2])
        q_electron_diff = q_electron_baseline - cgyro_q_electron_interp[i]
        q_electron_error = asinh(q_electron_diff / max(Q_ELECTRON_STD, 1e-10))

        q_ions_baseline = sum(QL_weights_baseline[:, 2:3, 1, idx, 2])
        q_ions_diff = q_ions_baseline - cgyro_q_ions_interp[i]
        q_ions_error = asinh(q_ions_diff / max(Q_IONS_STD, 1e-10))

        baseline_loss += GAMMA_WEIGHT * gamma_error^2 +
                        QL_COMPONENT_WEIGHT * (gamma_electron_error^2 + q_electron_error^2 + q_ions_error^2)
    end
    baseline_mse = baseline_loss / n_ky_opt

    println("\n" * "=" ^ 60)
    println("RESULTS")
    println("Time: $(round(elapsed_time, digits=2))s")
    println("Optimal loss: $(round(total_error/n_ky_opt, digits=6))")
    println("Baseline loss: $(round(baseline_mse, digits=6))")
    println("Error reduction: $(round(100*(1 - (total_error/n_ky_opt)/baseline_mse), digits=1))%")
    println("=" ^ 60)

    result = (converged=true, iterations=elapsed_time, minimum=total_error/n_ky_opt)
    return width_spectrum_optimal, optimal_widths, ky_grid_in, result, input_tjlf, mask, cgyro_gamma_interp, cgyro_freq_interp
end

function save_combined_results(results_dict::Dict, output_dir::String="width_spec_opt_multiflux")
    """Save optimization results to JSON file"""
    mkpath(output_dir)

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

    for shot in sort(collect(keys(results_dict)))
        result_tuple = results_dict[shot]
        width_spectrum_optimal, optimal_widths, ky_grid_in, result, input_tjlf, mask, cgyro_gamma_interp, cgyro_freq_interp = result_tuple

        isempty(optimal_widths) && continue

        ky_cgyro = Float64.(data["cgyro_growthrate_spectra"][shot][1])
        cgyro_gamma = Float64.(data["cgyro_growthrate_spectra"][shot][2])
        cgyro_freq = Float64.(data["cgyro_growthrate_spectra"][shot][3])

        length(ky_cgyro) < 2 && continue

        ky_grid = input_tjlf.KY_SPECTRUM
        baseline_width_spectrum = fill(input_tjlf.WIDTH, length(ky_grid))

        gamma_baseline, freq_baseline, _ = run_tjlf_with_width_spectrum(input_tjlf, baseline_width_spectrum)
        gamma_optimized, freq_optimized, _ = run_tjlf_with_width_spectrum(input_tjlf, width_spectrum_optimal)

        itp_gamma = interpolate((ky_cgyro,), cgyro_gamma, Gridded(Linear()))
        itp_freq = interpolate((ky_cgyro,), cgyro_freq, Gridded(Linear()))

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

    json_path = joinpath(output_dir, "combined_results.json")
    open(json_path, "w") do f
        JSON.print(f, combined_data, 4)
    end

    println("Saved results: $json_path ($(length(combined_data["id"])) shots)")
    return json_path
end

function plot_gamma_comparison(shot::Int, width_spectrum_optimal, ky_grid_in, input_tjlf)
    """Plot gamma comparison: CGYRO vs baseline vs optimized TJLF"""

    ky_cgyro = Float64.(data["cgyro_growthrate_spectra"][shot][1])
    cgyro_gamma = Float64.(data["cgyro_growthrate_spectra"][shot][2])
    ky_grid = input_tjlf.KY_SPECTRUM

    baseline_width_spectrum = fill(input_tjlf.WIDTH, length(ky_grid))
    gamma_baseline, _, _ = run_tjlf_with_width_spectrum(input_tjlf, baseline_width_spectrum)
    gamma_optimized, _, _ = run_tjlf_with_width_spectrum(input_tjlf, width_spectrum_optimal)

    plt = plot(
        ky_cgyro, cgyro_gamma,
        label = "CGYRO",
        marker = :circle,
        lw = 2,
        legend = :bottomright,
        fontfamily = "Computer Modern",
        xscale = :log10,
        yscale = :log10,
        ylims = (1e-4, Inf),
        xlabel = L"k_y",
        ylabel = L"\gamma",
        title = "Growth Rate - Shot $shot (4-component)",
        size = (800, 600)
    )

    plot!(plt, ky_grid, gamma_baseline,
        label = "TJLF baseline",
        linestyle = :dash,
        alpha = 0.7,
        lw = 2
    )

    plot!(plt, ky_grid, gamma_optimized,
        label = "TJLF optimized",
        marker = :diamond,
        lw = 2,
        color = :red
    )

    if !isempty(ky_grid_in)
        vline!(plt, [minimum(ky_grid_in)], linestyle = :dot, color = :gray, alpha = 0.5, label = "")
        vline!(plt, [maximum(ky_grid_in)], linestyle = :dot, color = :gray, alpha = 0.5, label = "")
    end

    vline!(plt, [0.17], linestyle = :dash, color = :green, alpha = 0.5, label = "")
    vline!(plt, [0.8], linestyle = :dash, color = :green, alpha = 0.5, label = "")

    return plt
end

function plot_frequency_comparison(shot::Int, width_spectrum_optimal, ky_grid_in, input_tjlf)
    """Plot frequency comparison with automatic sign-flip detection"""

    ky_cgyro = Float64.(data["cgyro_growthrate_spectra"][shot][1])
    cgyro_freq = Float64.(data["cgyro_growthrate_spectra"][shot][3])
    ky_grid = input_tjlf.KY_SPECTRUM

    baseline_width_spectrum = fill(input_tjlf.WIDTH, length(ky_grid))
    _, freq_baseline, _ = run_tjlf_with_width_spectrum(input_tjlf, baseline_width_spectrum)
    _, freq_optimized, _ = run_tjlf_with_width_spectrum(input_tjlf, width_spectrum_optimal)

    # Check sign agreement at high ky
    itp_freq = interpolate((ky_cgyro,), cgyro_freq, Gridded(Linear()))
    n_check = min(3, length(ky_grid))
    high_ky_indices = []
    for i in (length(ky_grid) - n_check + 1):length(ky_grid)
        if ky_grid[i] >= minimum(ky_cgyro) && ky_grid[i] <= maximum(ky_cgyro)
            push!(high_ky_indices, i)
        end
    end

    sign_flip_needed = false
    if !isempty(high_ky_indices) && length(high_ky_indices) >= 2
        cgyro_signs = [sign(itp_freq(ky_grid[i])) for i in high_ky_indices]
        tjlf_baseline_signs = [sign(freq_baseline[i]) for i in high_ky_indices]
        tjlf_opt_signs = [sign(freq_optimized[i]) for i in high_ky_indices]

        baseline_disagree = sum(cgyro_signs .!= tjlf_baseline_signs) > length(high_ky_indices) / 2
        opt_disagree = sum(cgyro_signs .!= tjlf_opt_signs) > length(high_ky_indices) / 2

        if baseline_disagree || opt_disagree
            sign_flip_needed = true
            freq_baseline = -freq_baseline
            freq_optimized = -freq_optimized
            println("  Sign flip applied to TJLF frequencies for Shot $shot")
        end
    end

    title_suffix = sign_flip_needed ? " (TJLF sign flipped)" : ""
    plt = plot(
        ky_cgyro, cgyro_freq,
        label = "CGYRO",
        marker = :circle,
        lw = 2,
        legend = :best,
        fontfamily = "Computer Modern",
        xscale = :log10,
        xlabel = L"k_y",
        ylabel = L"\omega",
        title = "Frequency - Shot $shot" * title_suffix,
        size = (800, 600)
    )

    plot!(plt, ky_grid, freq_baseline,
        label = "TJLF baseline",
        linestyle = :dash,
        alpha = 0.7,
        lw = 2
    )

    plot!(plt, ky_grid, freq_optimized,
        label = "TJLF optimized",
        marker = :diamond,
        lw = 2,
        color = :red
    )

    hline!(plt, [0.0], linestyle = :dot, color = :black, alpha = 0.3, label = "")

    if !isempty(ky_grid_in)
        vline!(plt, [minimum(ky_grid_in)], linestyle = :dot, color = :gray, alpha = 0.5, label = "")
        vline!(plt, [maximum(ky_grid_in)], linestyle = :dot, color = :gray, alpha = 0.5, label = "")
    end

    return plt
end

function plot_ql_weights_comparison(shot::Int, width_spectrum_optimal, ky_grid_in, input_tjlf)
    """Plot QL weights comparison: CGYRO vs baseline vs optimized TJLF for all 3 components"""

    # Get CGYRO flux data
    ky_cgyro = Float64.(data["cgyro_growthrate_spectra"][shot][1])
    cgyro_fluxes = data["cgyro_ql_fluxes"][shot]
    n_ky_cgyro = length(ky_cgyro)

    # Extract CGYRO QL weights
    cgyro_gamma_electron = zeros(Float64, n_ky_cgyro)
    cgyro_q_electron = zeros(Float64, n_ky_cgyro)
    cgyro_q_ions = zeros(Float64, n_ky_cgyro)

    for ky_idx in 1:n_ky_cgyro
        cgyro_gamma_electron[ky_idx] = sum(cgyro_fluxes[ky_idx][3][field+1][1] for field in 0:2)
        cgyro_q_electron[ky_idx] = sum(cgyro_fluxes[ky_idx][3][field+1][2] for field in 0:2)
        cgyro_q_ions[ky_idx] = sum(cgyro_fluxes[ky_idx][species+1][field+1][2]
                                   for species in 0:1, field in 0:2)
    end

    # Run TJLF baseline and optimized
    ky_grid = input_tjlf.KY_SPECTRUM
    baseline_width_spectrum = fill(input_tjlf.WIDTH, length(ky_grid))
    _, _, QL_weights_baseline = run_tjlf_with_width_spectrum(input_tjlf, baseline_width_spectrum)
    _, _, QL_weights_optimized = run_tjlf_with_width_spectrum(input_tjlf, width_spectrum_optimal)

    # Extract TJLF QL weights
    gamma_electron_baseline = [sum(QL_weights_baseline[:, 1, 1, ky_idx, 1]) for ky_idx in 1:length(ky_grid)]
    q_electron_baseline = [sum(QL_weights_baseline[:, 1, 1, ky_idx, 2]) for ky_idx in 1:length(ky_grid)]
    q_ions_baseline = [sum(QL_weights_baseline[:, 2:3, 1, ky_idx, 2]) for ky_idx in 1:length(ky_grid)]

    gamma_electron_optimized = [sum(QL_weights_optimized[:, 1, 1, ky_idx, 1]) for ky_idx in 1:length(ky_grid)]
    q_electron_optimized = [sum(QL_weights_optimized[:, 1, 1, ky_idx, 2]) for ky_idx in 1:length(ky_grid)]
    q_ions_optimized = [sum(QL_weights_optimized[:, 2:3, 1, ky_idx, 2]) for ky_idx in 1:length(ky_grid)]

    # Create 3-panel plot
    p1 = plot(
        ky_cgyro, cgyro_gamma_electron,
        label = "CGYRO",
        marker = :circle,
        lw = 2,
        legend = :topright,
        fontfamily = "Computer Modern",
        xscale = :log10,
        yscale = :log10,
        ylims = (1e-4, Inf),
        xlabel = L"k_y",
        ylabel = L"\Gamma_{e}",
        title = "Electron Particle QL Weight",
        size = (800, 400)
    )
    plot!(p1, ky_grid, gamma_electron_baseline, label = "TJLF baseline", linestyle = :dash, alpha = 0.7, lw = 2)
    plot!(p1, ky_grid, gamma_electron_optimized, label = "TJLF optimized", marker = :diamond, lw = 2, color = :red)
    if !isempty(ky_grid_in)
        vline!(p1, [minimum(ky_grid_in)], linestyle = :dot, color = :gray, alpha = 0.5, label = "")
        vline!(p1, [maximum(ky_grid_in)], linestyle = :dot, color = :gray, alpha = 0.5, label = "")
    end

    p2 = plot(
        ky_cgyro, cgyro_q_electron,
        label = "CGYRO",
        marker = :circle,
        lw = 2,
        legend = :topright,
        fontfamily = "Computer Modern",
        xscale = :log10,
        yscale = :log10,
        ylims = (1e-4, Inf),
        xlabel = L"k_y",
        ylabel = L"Q_{e}",
        title = "Electron Heat QL Weight",
        size = (800, 400)
    )
    plot!(p2, ky_grid, q_electron_baseline, label = "TJLF baseline", linestyle = :dash, alpha = 0.7, lw = 2)
    plot!(p2, ky_grid, q_electron_optimized, label = "TJLF optimized", marker = :diamond, lw = 2, color = :red)
    if !isempty(ky_grid_in)
        vline!(p2, [minimum(ky_grid_in)], linestyle = :dot, color = :gray, alpha = 0.5, label = "")
        vline!(p2, [maximum(ky_grid_in)], linestyle = :dot, color = :gray, alpha = 0.5, label = "")
    end

    p3 = plot(
        ky_cgyro, cgyro_q_ions,
        label = "CGYRO",
        marker = :circle,
        lw = 2,
        legend = :topright,
        fontfamily = "Computer Modern",
        xscale = :log10,
        yscale = :log10,
        ylims = (1e-4, Inf),
        xlabel = L"k_y",
        ylabel = L"Q_{ions}",
        title = "Ion Heat QL Weight",
        size = (800, 400)
    )
    plot!(p3, ky_grid, q_ions_baseline, label = "TJLF baseline", linestyle = :dash, alpha = 0.7, lw = 2)
    plot!(p3, ky_grid, q_ions_optimized, label = "TJLF optimized", marker = :diamond, lw = 2, color = :red)
    if !isempty(ky_grid_in)
        vline!(p3, [minimum(ky_grid_in)], linestyle = :dot, color = :gray, alpha = 0.5, label = "")
        vline!(p3, [maximum(ky_grid_in)], linestyle = :dot, color = :gray, alpha = 0.5, label = "")
    end

    # Combine into single figure
    plt = plot(p1, p2, p3, layout = (3, 1), size = (800, 1200))
    return plt
end

function main_parallel(
    shots::Vector{Int};
    ky_filter_mode::String="less_than",
    ky_threshold_low::Union{Float64,Nothing}=nothing,
    ky_threshold_high::Union{Float64,Nothing}=nothing,
    lambda::Float64=0.01
)
    """Run 4-component QL weights optimization for multiple shots in parallel"""

    gamma_pct = round(GAMMA_WEIGHT * 100, digits=1)
    ql_pct = round(QL_WEIGHT * 100, digits=1)
    ql_component_pct = round(QL_COMPONENT_WEIGHT * 100, digits=2)

    println("=" ^ 80)
    println("PARALLEL WIDTH_SPECTRUM OPTIMIZATION (4-COMPONENT QL WEIGHTS)")
    println("Workers: $(nworkers()) | Shots: $shots")
    println("Cost: $(gamma_pct)% gamma + $(ql_pct)% QL weights ($(ql_component_pct)% each: Γ_e, Q_e, Q_ions)")
    println("=" ^ 80)

    optimize_func = shot -> optimize_width_spectrum_per_ky(shot;
        ky_filter_mode=ky_filter_mode,
        ky_threshold_low=ky_threshold_low,
        ky_threshold_high=ky_threshold_high,
        lambda=lambda,
        use_parallel=true
    )

    start_time = time()
    results = pmap(optimize_func, shots;
                   retry_delays=ExponentialBackOff(n=3, first_delay=1.0, max_delay=10.0),
                   on_error=ex -> (println("Warning: Shot failed: $ex"); nothing))

    any(isnothing, results) && error("Some shot optimizations failed")

    elapsed_time = time() - start_time
    println("\n" * "=" ^ 80)
    println("ALL COMPLETE | Total: $(round(elapsed_time, digits=2))s | Avg: $(round(elapsed_time/length(shots), digits=2))s/shot")
    println("=" ^ 80)

    results_dict = Dict(shot => results[i] for (i, shot) in enumerate(shots))

    # Save results
    output_dir = "width_spec_opt_multiflux_n95"
    println("\nSaving results...")
    save_combined_results(results_dict, output_dir)

    # Generate plots
    println("\nGenerating plots...")
    shot_figs_dir = joinpath(output_dir, "shot_figs")
    mkpath(shot_figs_dir)

    for shot in shots
        width_spectrum_optimal, optimal_widths, ky_grid_in, result, input_tjlf, mask, cgyro_gamma_interp, cgyro_freq_interp = results_dict[shot]

        if isempty(optimal_widths)
            println("Shot $shot: Skipped (insufficient data)")
            continue
        end

        println("Shot $shot:")

        plt_gamma = plot_gamma_comparison(shot, width_spectrum_optimal, ky_grid_in, input_tjlf)
        gamma_filename = joinpath(shot_figs_dir, "shot_$(shot)_gamma.png")
        savefig(plt_gamma, gamma_filename)
        println("  Saved: $gamma_filename")

        plt_freq = plot_frequency_comparison(shot, width_spectrum_optimal, ky_grid_in, input_tjlf)
        freq_filename = joinpath(shot_figs_dir, "shot_$(shot)_frequency.png")
        savefig(plt_freq, freq_filename)
        println("  Saved: $freq_filename")

        plt_ql = plot_ql_weights_comparison(shot, width_spectrum_optimal, ky_grid_in, input_tjlf)
        ql_filename = joinpath(shot_figs_dir, "shot_$(shot)_ql_weights.png")
        savefig(plt_ql, ql_filename)
        println("  Saved: $ql_filename")
    end

    println("\n" * "=" ^ 80)
    println("Results saved to: $output_dir")
    println("=" ^ 80)

    return results_dict
end

gamma_pct = round(GAMMA_WEIGHT * 100, digits=1)
ql_pct = round(QL_WEIGHT * 100, digits=1)
ql_component_pct = round(QL_COMPONENT_WEIGHT * 100, digits=2)

println("Functions loaded successfully!")
println("Workers: $(nworkers()) | Version: 2024-11-26-v4-normalized")
println("Cost: $(gamma_pct)% gamma + $(ql_pct)% QL weights ($(ql_component_pct)% each: Γ_e, Q_e, Q_ions)")
println("Normalization: asinh((TJLF - CGYRO) / STD) - dataset-specific")
println("=" ^ 60)

# ============================================================================
# EXECUTION
# ============================================================================

# shots_to_optimize = [8,10,11,12]  # Test with single shot
shots_to_optimize = collect(1:20)  # All shots
# shots_to_optimize = [1,8,10,11,12]  # Specific shots
# shots_to_optimize = vcat(1:7,9,13:20)

width_spectrum_results = main_parallel(shots_to_optimize,
    ky_filter_mode="full",
    lambda=0.0
)
