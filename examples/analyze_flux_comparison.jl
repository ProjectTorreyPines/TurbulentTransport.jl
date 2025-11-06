#!/usr/bin/env julia
# ============================================================================
# Total Flux Comparison Analysis Script
# ============================================================================
# Loads optimized WIDTH_SPECTRUM from combined_results.json,
# runs TJLF with both baseline and optimized widths,
# and compares total heat fluxes vs CGYRO
#
# Heat flux is the total energy flux summed over all species and fields
# - Q_total: total heat flux (electron + ions + impurities)
# ============================================================================

using Pkg
Pkg.activate("/global/homes/t/trifolio/.julia/dev/TurbulentTransport")

using TurbulentTransport
using JSON
using Plots
using LaTeXStrings
using TJLF
using Statistics

println("="^80)
println("TOTAL FLUX COMPARISON ANALYSIS")
println("="^80)

# Accept results directory as command-line argument
results_dir = length(ARGS) >= 1 ? ARGS[1] : "width_spec_opt_multiflux_80g/"
if !endswith(results_dir, "/")
    results_dir = results_dir * "/"
end

# ============================================================================
# USER PLOT TITLE CUSTOMIZATION
# ============================================================================
# Add custom text to the plot title (leave empty for default titles only)
plot_title_suffix = " 80% gamma"  # e.g., " - 25% gamma" or " - Final Results"
# ============================================================================

# Load the optimized results
results_file = results_dir * "combined_results.json"
println("\nLoading optimized results from: $results_file")
if !isfile(results_file)
    error("Results file not found: $results_file")
end
opt_results = JSON.parsefile(results_file)

# Load original CGYRO data
data_file = "qlnn_training_subset_merged_result_dict_dmn36p5_dict.json"
println("Loading CGYRO data from: $data_file")
data = JSON.parsefile(data_file)

# Template input.tglf for setting up TJLF runs
input_tglf_lines = """
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

function setup_input_file(shot::Int, filepath::String)
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

function run_tjlf_with_width_spectrum(input_tjlf, width_spectrum::Vector{Float64})
    """Run TJLF with custom WIDTH_SPECTRUM and return total heat flux

    Returns: Q_total (total heat flux for all species)
    - Q_total: total heat flux summed over all fields and species
    """

    # Create a copy to avoid modifying the original
    input = deepcopy(input_tjlf)

    # Set the WIDTH_SPECTRUM
    input.WIDTH_SPECTRUM = width_spectrum

    satParams = get_sat_params(input)
    input.KY_SPECTRUM .= get_ky_spectrum(input, satParams.grad_r0)
    outputHermite = gauss_hermite(input)
    QL_weights, eigenvalue = tjlf_TM(input, satParams, outputHermite)

    # Compute fluxes using sum_ky_spectrum
    # QL_flux_out has shape (field, species, type)
    # type: 1=particle, 2=energy (heat), 3=toroidal stress, 4=parallel stress, 5=exchange
    QL_flux_out, flux_out = sum_ky_spectrum(input, satParams, eigenvalue[:, :, 1], QL_weights)

    # Q_total: sum over all fields and all species for heat flux (type=2)
    Q_total = sum(QL_flux_out[:, :, 2])

    return Q_total
end

function compute_cgyro_flux_totals(shot_idx::Int)
    """Compute total heat flux from CGYRO data

    Returns: Q_total (total heat flux)
    """

    # CGYRO fluxes stored in QLGYRO_SAT3_flux_Q_azf-1
    # This is a 3-element array: [electron, ion, carbon]
    cgyro_flux_Q = data["QLGYRO_SAT3_flux_Q_azf-1"][shot_idx]

    # Sum all species to get total heat flux
    Q_total = sum(cgyro_flux_Q)

    return Q_total
end

# Main analysis loop
println("\n" * "="^80)
println("ANALYZING FLUX COMPARISON FOR ALL SHOTS")
println("="^80)

n_shots = length(opt_results["id"])
println("Number of shots in optimized results: $n_shots")

# Storage for results
shot_ids = String[]
cgyro_Q = Float64[]
baseline_Q = Float64[]
optimized_Q = Float64[]

# Process each shot
for shot_idx in 1:n_shots
    shot_id = opt_results["id"][shot_idx]
    println("\n" * "-"^60)
    println("Shot $shot_idx (ID: $shot_id)")
    println("-"^60)

    # Setup TJLF input for this shot
    tmpdir = mktempdir()
    filepath = joinpath(tmpdir, "input.tglf")
    setup_input_file(shot_idx, filepath)

    # Load and setup TJLF
    input_tglf = TurbulentTransport.load(InputTGLF(), filepath)
    TurbulentTransport.apply_presets!(input_tglf)
    TurbulentTransport.run_tjlf(input_tglf)

    input_tjlf = readInput(filepath)
    satParams = get_sat_params(input_tjlf)
    input_tjlf.KY_SPECTRUM .= get_ky_spectrum(input_tjlf, satParams.grad_r0)

    # Get CGYRO flux totals
    Q_cgyro = compute_cgyro_flux_totals(shot_idx)

    # Get optimized WIDTH_SPECTRUM from results
    width_optimized = Float64.(opt_results["width_optimized_spectrum"][shot_idx])
    width_baseline = Float64.(opt_results["width_baseline"][shot_idx])

    # Make sure widths match the ky grid size
    ky_grid = input_tjlf.KY_SPECTRUM
    n_ky = length(ky_grid)

    # If optimized width is shorter, pad with baseline values
    if length(width_optimized) < n_ky
        println("  WARNING: width_optimized has $(length(width_optimized)) values, ky_grid has $n_ky points")
        println("           Padding with baseline width")
        width_opt_full = fill(input_tjlf.WIDTH, n_ky)
        width_opt_full[1:length(width_optimized)] = width_optimized
        width_optimized = width_opt_full
    end

    if length(width_baseline) < n_ky
        width_base_full = fill(input_tjlf.WIDTH, n_ky)
        width_base_full[1:length(width_baseline)] = width_baseline
        width_baseline = width_base_full
    end

    # Run TJLF with baseline WIDTH_SPECTRUM
    println("  Running TJLF with baseline WIDTH...")
    Q_base = run_tjlf_with_width_spectrum(input_tjlf, width_baseline)

    # Run TJLF with optimized WIDTH_SPECTRUM
    println("  Running TJLF with optimized WIDTH...")
    Q_opt = run_tjlf_with_width_spectrum(input_tjlf, width_optimized)

    # Store results
    push!(shot_ids, shot_id)
    push!(cgyro_Q, Q_cgyro)
    push!(baseline_Q, Q_base)
    push!(optimized_Q, Q_opt)

    # Print comparison
    Q_base_error = 100 * abs(Q_base - Q_cgyro) / abs(Q_cgyro)
    Q_opt_error = 100 * abs(Q_opt - Q_cgyro) / abs(Q_cgyro)

    println("\n  Total Heat Flux (Q_total):")
    println("    CGYRO:     $(round(Q_cgyro, digits=6))")
    println("    Baseline:  $(round(Q_base, digits=6))  (error: $(round(Q_base_error, digits=1))%)")
    println("    Optimized: $(round(Q_opt, digits=6))  (error: $(round(Q_opt_error, digits=1))%)")

    # Clean up temp directory
    rm(tmpdir, recursive=true)
end

println("\n" * "="^80)
println("CREATING COMPARISON PLOT")
println("="^80)

# Calculate errors only for non-zero CGYRO values
function safe_relative_error(predicted, target, max_error=1000.0)
    errors = Float64[]
    for i in 1:length(predicted)
        if target[i] != 0.0
            error_pct = 100 * abs(predicted[i] - target[i]) / abs(target[i])
            if error_pct <= max_error
                push!(errors, error_pct)
            end
        end
    end
    return errors
end

# Filter out zeros and extreme outliers (>1000% error) for log-log plotting
max_error_pct = 1000.0

# Heat flux filtering
Q_valid_baseline = (cgyro_Q .!= 0) .& (baseline_Q .!= 0) .&
                   (100 * abs.(baseline_Q .- cgyro_Q) ./ abs.(cgyro_Q) .<= max_error_pct)
Q_valid_optimized = (cgyro_Q .!= 0) .& (optimized_Q .!= 0) .&
                    (100 * abs.(optimized_Q .- cgyro_Q) ./ abs.(cgyro_Q) .<= max_error_pct)

# Find common valid points (valid in BOTH baseline and optimized)
Q_valid_common = Q_valid_baseline .& Q_valid_optimized

Q_cgyro_base_filtered = abs.(cgyro_Q[Q_valid_baseline])
Q_base_filtered = abs.(baseline_Q[Q_valid_baseline])
Q_cgyro_opt_filtered = abs.(cgyro_Q[Q_valid_optimized])
Q_opt_filtered = abs.(optimized_Q[Q_valid_optimized])

n_Q_filtered_base = sum(.!Q_valid_baseline)
n_Q_filtered_opt = sum(.!Q_valid_optimized)

# Calculate errors on common valid points for fair comparison
Q_base_errors_common = Float64[]
Q_opt_errors_common = Float64[]
for i in 1:length(cgyro_Q)
    if Q_valid_common[i]
        push!(Q_base_errors_common, 100 * abs(baseline_Q[i] - cgyro_Q[i]) / abs(cgyro_Q[i]))
        push!(Q_opt_errors_common, 100 * abs(optimized_Q[i] - cgyro_Q[i]) / abs(cgyro_Q[i]))
    end
end

# Also calculate errors for all valid points (for summary statistics)
Q_base_errors = safe_relative_error(baseline_Q, cgyro_Q)
Q_opt_errors = safe_relative_error(optimized_Q, cgyro_Q)

println("Heat flux: Filtered $n_Q_filtered_base baseline points and $n_Q_filtered_opt optimized points (zeros + outliers >$(max_error_pct)%)")
println("           Common valid points for improvement calculation: $(length(Q_base_errors_common))")

# Create scatter plot: TJLF vs CGYRO (diagonal = perfect agreement)
p = plot(size=(800, 700), margin=5Plots.mm)

# Find min/max for diagonal line (heat flux)
if !isempty(Q_cgyro_base_filtered) && !isempty(Q_cgyro_opt_filtered)
    Q_min = min(minimum(Q_cgyro_base_filtered), minimum(Q_base_filtered),
                minimum(Q_cgyro_opt_filtered), minimum(Q_opt_filtered))
    Q_max = max(maximum(Q_cgyro_base_filtered), maximum(Q_base_filtered),
                maximum(Q_cgyro_opt_filtered), maximum(Q_opt_filtered))

    # Heat flux scatter plot - LOG-LOG scale
    # Create title with optional user suffix
    Q_title = "Total Heat Flux: TJLF vs CGYRO (log-log)"
    if !isempty(plot_title_suffix)
        Q_title = Q_title * " - " * plot_title_suffix
    end

    plot!(p, [Q_min, Q_max], [Q_min, Q_max],
          label="Perfect Agreement", color=:black, ls=:dash, lw=2, alpha=0.5,
          xscale=:log10, yscale=:log10)
    scatter!(p, Q_cgyro_base_filtered, Q_base_filtered,
          label="TJLF Baseline ($(length(Q_base_filtered)) pts)", marker=:square, ms=6, alpha=0.7,
          xlabel="CGYRO Total Heat Flux (Q)", ylabel="TJLF Total Heat Flux (Q)",
          title=Q_title,
          legend=:topleft, grid=true)
    scatter!(p, Q_cgyro_opt_filtered, Q_opt_filtered,
          label="TJLF Optimized ($(length(Q_opt_filtered)) pts)", marker=:diamond, ms=7, color=:red)

    # Add improvement statistic to plot (using common valid points)
    # if !isempty(Q_base_errors_common) && !isempty(Q_opt_errors_common)
    #     Q_improvement = mean(Q_base_errors_common) - mean(Q_opt_errors_common)
    #     improvement_text = "Overall Improvement: $(round(Q_improvement, digits=1))%"
    #     n_common = length(Q_base_errors_common)
    #     n_filtered = length(cgyro_Q) - n_common
    #     annotate!(p, :bottomright, text("$improvement_text\n(on $n_common common points)\n\n$n_filtered points excluded\n(zeros + outliers >1000%)", 9, :darkblue, :right))
    # elseif n_Q_filtered_base > 0 || n_Q_filtered_opt > 0
    #     annotate!(p, :bottomright, text("$(max(n_Q_filtered_base, n_Q_filtered_opt)) points excluded\n(zeros + outliers >1000%)", 8, :gray))
    # end
else
    plot!(p, title="Heat Flux: No valid data (all zeros)")
end

# Save plot
output_file = results_dir * "flux_comparison_scatter.png"
savefig(p, output_file)
println("\nSaved plot to: $output_file")

# Calculate and print summary statistics
println("\n" * "="^80)
println("SUMMARY STATISTICS")
println("="^80)

# Count filtered points
n_Q_zeros = sum(cgyro_Q .== 0)
n_Q_outliers = length(cgyro_Q) - n_Q_zeros - length(Q_base_errors)

println("\nTotal Heat Flux (Q) Errors:")
println("  Excluded: $n_Q_zeros zeros, $n_Q_outliers outliers (>1000% error)")
println("  Valid points for statistics: $(length(Q_base_errors))/$(length(cgyro_Q))")
if !isempty(Q_base_errors)
    println("  Baseline  - Mean: $(round(mean(Q_base_errors), digits=1))%, Median: $(round(median(Q_base_errors), digits=1))%, Max: $(round(maximum(Q_base_errors), digits=1))%")
    println("  Optimized - Mean: $(round(mean(Q_opt_errors), digits=1))%, Median: $(round(median(Q_opt_errors), digits=1))%, Max: $(round(maximum(Q_opt_errors), digits=1))%")
    println("  Improvement (on all valid points): $(round(mean(Q_base_errors) - mean(Q_opt_errors), digits=1))%")
else
    println("  No valid data (all values filtered out)")
end

# Additional statistics on common valid points
if !isempty(Q_base_errors_common) && !isempty(Q_opt_errors_common)
    println("\n  Common Valid Points (valid in BOTH baseline and optimized): $(length(Q_base_errors_common))/$(length(cgyro_Q))")
    println("  Baseline  - Mean: $(round(mean(Q_base_errors_common), digits=1))%, Median: $(round(median(Q_base_errors_common), digits=1))%")
    println("  Optimized - Mean: $(round(mean(Q_opt_errors_common), digits=1))%, Median: $(round(median(Q_opt_errors_common), digits=1))%")
    println("  Improvement (on common points): $(round(mean(Q_base_errors_common) - mean(Q_opt_errors_common), digits=1))%")
end

# Save results to JSON
println("\n" * "="^80)
println("SAVING RESULTS TO JSON")
println("="^80)

# Create per-shot error arrays
Q_base_errors_per_shot = []
Q_opt_errors_per_shot = []

for i in 1:length(cgyro_Q)
    if cgyro_Q[i] != 0.0
        push!(Q_base_errors_per_shot, 100 * abs(baseline_Q[i] - cgyro_Q[i]) / abs(cgyro_Q[i]))
        push!(Q_opt_errors_per_shot, 100 * abs(optimized_Q[i] - cgyro_Q[i]) / abs(cgyro_Q[i]))
    else
        push!(Q_base_errors_per_shot, nothing)
        push!(Q_opt_errors_per_shot, nothing)
    end
end

results_summary = Dict(
    "shots" => shot_ids,
    "heat_flux" => Dict(
        "cgyro" => cgyro_Q,
        "tjlf_baseline" => baseline_Q,
        "tjlf_optimized" => optimized_Q,
        "baseline_error_percent" => Q_base_errors_per_shot,
        "optimized_error_percent" => Q_opt_errors_per_shot
    ),
    "summary_statistics" => Dict(
        "heat_flux" => Dict(
            "n_valid_shots" => length(Q_base_errors),
            "n_zero_cgyro" => length(cgyro_Q) - length(Q_base_errors),
            "baseline_mean_error" => isempty(Q_base_errors) ? nothing : mean(Q_base_errors),
            "baseline_median_error" => isempty(Q_base_errors) ? nothing : median(Q_base_errors),
            "baseline_max_error" => isempty(Q_base_errors) ? nothing : maximum(Q_base_errors),
            "optimized_mean_error" => isempty(Q_opt_errors) ? nothing : mean(Q_opt_errors),
            "optimized_median_error" => isempty(Q_opt_errors) ? nothing : median(Q_opt_errors),
            "optimized_max_error" => isempty(Q_opt_errors) ? nothing : maximum(Q_opt_errors),
            "improvement" => isempty(Q_base_errors) ? nothing : mean(Q_base_errors) - mean(Q_opt_errors),
            "n_common_valid_shots" => length(Q_base_errors_common),
            "common_baseline_mean_error" => isempty(Q_base_errors_common) ? nothing : mean(Q_base_errors_common),
            "common_optimized_mean_error" => isempty(Q_opt_errors_common) ? nothing : mean(Q_opt_errors_common),
            "common_improvement" => isempty(Q_base_errors_common) ? nothing : mean(Q_base_errors_common) - mean(Q_opt_errors_common)
        )
    )
)

json_output = results_dir * "flux_comparison_results.json"
open(json_output, "w") do f
    JSON.print(f, results_summary, 4)
end
println("Saved results to: $json_output")

println("\n" * "="^80)
println("ANALYSIS COMPLETE")
println("="^80)
