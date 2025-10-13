using TurbulentTransport
using JSON
using Plots
using LaTeXStrings
using TJLF
using Interpolations
using Optim
using Statistics

# Load data
data = JSON.parsefile("qlnn_training_subset_merged_result_dict_dmn36p5_dict.json")

# Template input.tglf content
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
C_B = 0.315
SIG_B = 0.34
BOUNCE_COEFF = 3.0
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
    """Run TJLF with a custom WIDTH_SPECTRUM vector"""

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

    return eigenvalue[1, :, 1]  # Return gamma values
end

function optimize_width_spectrum(shot::Int)
    """Optimize WIDTH_SPECTRUM vector for a specific shot"""

    println("=" ^ 60)
    println("OPTIMIZING WIDTH_SPECTRUM FOR SHOT $shot")
    println("Shot ID: $(data["id"][shot])")
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
    mask = (ky_grid .>= minimum(ky_cgyro)) .& (ky_grid .<= maximum(ky_cgyro)) .& (ky_grid .< 0.8)
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
            gamma_tjlf = run_tjlf_with_width_spectrum(input_tjlf, width_spectrum_full)
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

function plot_width_spectrum_results(shot::Int, width_spectrum_optimal, optimal_widths, ky_grid_in, input_tjlf, mask, cgyro_gamma_interp)
    """Plot the optimization results"""

    # Get CGYRO data for plotting
    ky_cgyro = Float64.(data["cgyro_growthrate_spectra"][shot][1])
    cgyro_gamma = Float64.(data["cgyro_growthrate_spectra"][shot][2])
    ky_grid = input_tjlf.KY_SPECTRUM

    # Run baseline TJLF
    baseline_width_spectrum = fill(input_tjlf.WIDTH, length(ky_grid))
    gamma_baseline = run_tjlf_with_width_spectrum(input_tjlf, baseline_width_spectrum)

    # Run optimized TJLF
    gamma_optimized = run_tjlf_with_width_spectrum(input_tjlf, width_spectrum_optimal)

    # Create main comparison plot
    plt1 = plot(
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
        ylabel = L"\\gamma",
        title = "WIDTH_SPECTRUM Optimization Results - Shot $shot"
    )

    # Plot baseline TJLF
    plot!(plt1, ky_grid, gamma_baseline,
        label = "TJLF baseline (WIDTH=$(round(input_tjlf.WIDTH, digits=2)))",
        linestyle = :dash,
        alpha = 0.7,
        lw = 2
    )

    # Plot optimized TJLF
    plot!(plt1, ky_grid, gamma_optimized,
        label = "TJLF optimized WIDTH_SPECTRUM",
        marker = :diamond,
        lw = 2,
        color = :red
    )

    # Highlight optimization region
    vline!(plt1, [minimum(ky_grid_in)], linestyle = :dot, color = :gray, alpha = 0.5, label = "")
    vline!(plt1, [maximum(ky_grid_in)], linestyle = :dot, color = :gray, alpha = 0.5, label = "")

    # Create WIDTH_SPECTRUM plot
    plt2 = plot(
        ky_grid, width_spectrum_optimal,
        label = "Optimized WIDTH_SPECTRUM",
        marker = :circle,
        lw = 2,
        xlabel = L"k_y",
        ylabel = "WIDTH value",
        title = "Optimized WIDTH_SPECTRUM vs ky",
        xscale = :log10,
        legend = :best
    )

    # Show baseline WIDTH as horizontal line
    hline!(plt2, [input_tjlf.WIDTH], label = "Baseline WIDTH", linestyle = :dash, alpha = 0.7)

    # Highlight optimization region
    vline!(plt2, [minimum(ky_grid_in)], linestyle = :dot, color = :gray, alpha = 0.5, label = "")
    vline!(plt2, [maximum(ky_grid_in)], linestyle = :dot, color = :gray, alpha = 0.5, label = "")

    # Combine plots
    combined_plot = plot(plt1, plt2, layout = (2, 1), size = (800, 800))

    return combined_plot
end

# Main execution
function main_fn()
    shot = 11  # Change this to test different shots

    # Run optimization
    width_spectrum_optimal, optimal_widths, ky_grid_in, result, input_tjlf, mask, cgyro_gamma_interp = optimize_width_spectrum(shot)

    # Plot results
    plt = plot_width_spectrum_results(shot, width_spectrum_optimal, optimal_widths, ky_grid_in, input_tjlf, mask, cgyro_gamma_interp)
    display(plt)

    # Save plot
    savefig("width_spectrum_optimization_shot_$shot.png")
    println("\nPlot saved as: width_spectrum_optimization_shot_$shot.png")

    return width_spectrum_optimal, result
end

println("Functions defined: ", names(Main))

# Run the optimization
width_spectrum_result = main_fn()
