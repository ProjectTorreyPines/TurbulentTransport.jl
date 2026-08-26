function Base.show(io::IO, ::MIME"text/plain", input_tglf::InputTGLF)
    for fname in sort!(collect(fieldnames(typeof(input_tglf))))
        value = getfield(input_tglf, fname)
        species_idx = tryparse(Int, split(string(fname), "_")[end])
        if value !== missing && (species_idx === nothing || species_idx <= input_tglf.NS)
            println(io, "$fname = $(value)")
        end
    end
end

"""
    InputTGLF(dd::IMAS.dd, rho::AbstractVector{Float64}, sat::Symbol=:sat0, electromagnetic::Bool=false, lump_ions::Bool=true)

Evaluate TGLF input parameters at given radii
"""
function InputTGLF(dd::IMAS.dd, rho::AbstractVector{Float64}, sat::Symbol=:sat0, electromagnetic::Bool=false, lump_ions::Bool=true; MXH_modes::Int=1)
    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]
    gridpoint_cp = [argmin_abs(cp1d.grid.rho_tor_norm, ρ) for ρ in rho]
    return InputTGLF(eqt, cp1d, gridpoint_cp, sat, electromagnetic, lump_ions; MXH_modes, wall=dd.wall)
end

"""
    InputTGLF(dd::IMAS.dd, gridpoint_cp::AbstractVector{Int}, sat::Symbol=:sat0, electromagnetic::Bool=false, lump_ions::Bool=true)

Evaluate TGLF input parameters at given core profiles grid indexes
"""
function InputTGLF(dd::IMAS.dd, gridpoint_cp::AbstractVector{Int}, sat::Symbol=:sat0, electromagnetic::Bool=false, lump_ions::Bool=true; MXH_modes::Int=1)
    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]
    return InputTGLF(eqt, cp1d, gridpoint_cp, sat, electromagnetic, lump_ions; MXH_modes, wall=dd.wall)
end

function InputTGLF(
    eqt::IMAS.equilibrium__time_slice{T},
    cp1d::IMAS.core_profiles__profiles_1d{T},
    gridpoint_cp::AbstractVector{Int},
    sat::Symbol,
    electromagnetic::Bool,
    lump_ions::Bool;
    MXH_modes::Int=1,
    wall=nothing) where {T<:Real}

    e = IMAS.cgs.e # statcoul
    k = IMAS.cgs.k # erg/eV
    mp = IMAS.cgs.mp # g
    me = IMAS.cgs.me # g
    md = IMAS.cgs.md # g
    m_to_cm = IMAS.cgs.m_to_cm
    m³_to_cm³ = IMAS.cgs.m³_to_cm³
    T_to_Gauss = IMAS.cgs.T_to_Gauss

    eqt1d = eqt.profiles_1d
    rho_cp = cp1d.grid.rho_tor_norm
    rho_eq = eqt1d.rho_tor_norm

    if lump_ions
        ions = IMAS.lump_ions_as_bulk_and_impurity(cp1d)
    else
        ions = cp1d.ion
    end

    Rmaj = IMAS.interp1d(rho_eq, m_to_cm * 0.5 * (eqt1d.r_outboard .+ eqt1d.r_inboard)).(rho_cp)

    rmin = GACODE.r_min_core_profiles(eqt1d, rho_cp)

    q_profile = IMAS.interp1d(rho_eq, eqt1d.q).(rho_cp)

    if MXH_modes > 1
        coeffs = IMAS.rec_to_mxh_coeffs(eqt, wall; MXH_modes)
        ns = size(coeffs, 2)
        rho_eq_s = rho_eq[1:ns]
        L = MXH_modes
        n_shape = min(L, 6)

        # Override Rmaj with MXH R0
        Rmaj = IMAS.interp1d(rho_eq_s, m_to_cm .* coeffs[1, :]).(rho_cp)

        kappa  = IMAS.interp1d(rho_eq_s, coeffs[4, :]).(rho_cp)
        s1_mxh = IMAS.interp1d(rho_eq_s, coeffs[6+L, :]).(rho_cp)  # asin(δ) = s[1]
        s2_mxh = IMAS.interp1d(rho_eq_s, coeffs[7+L, :]).(rho_cp)  # -ζ = s[2]

        delta = sin.(s1_mxh)
        zeta  = .-s2_mxh

        skappa = rmin .* IMAS.gradient(rmin, kappa) ./ kappa
        sdelta = rmin .* IMAS.gradient(rmin, s1_mxh)   # r * d(asin(δ))/dr
        szeta  = .-rmin .* IMAS.gradient(rmin, s2_mxh) # r * dζ/dr

        # cos coefficients c[0..n_shape] and shears: SHAPE_COS0..SHAPE_COS{n_shape}
        c_mxh  = [IMAS.interp1d(rho_eq_s, coeffs[5+m, :]).(rho_cp) for m in 0:n_shape]
        dc_mxh = [rmin .* IMAS.gradient(rmin, c_mxh[i]) for i in eachindex(c_mxh)]

        # sin coefficients s[3..n_shape] and shears: SHAPE_SIN3..SHAPE_SIN{n_shape}
        s_mxh  = [IMAS.interp1d(rho_eq_s, coeffs[5+L+m, :]).(rho_cp) for m in 3:n_shape]
        ds_mxh = [rmin .* IMAS.gradient(rmin, s_mxh[i]) for i in eachindex(s_mxh)]
    else
        if !ismissing(eqt1d, :elongation)
            kappa = IMAS.interp1d(rho_eq, eqt1d.elongation).(rho_cp)
        else
            kappa = zero(rho_cp)
        end

        if !ismissing(eqt1d, :triangularity_lower) && !ismissing(eqt1d, :triangularity_upper)
            delta = IMAS.interp1d(rho_eq, 0.5 * (eqt1d.triangularity_lower + eqt1d.triangularity_upper)).(rho_cp)
        else
            delta = zero(rho_cp)
        end

        if !ismissing(eqt1d, :squareness_lower_inner) && !ismissing(eqt1d, :squareness_lower_outer) && !ismissing(eqt1d, :squareness_upper_inner) &&
        !ismissing(eqt1d, :squareness_upper_outer)
            tmp = 0.25 .* (eqt1d.squareness_lower_inner .+ eqt1d.squareness_lower_outer .+ eqt1d.squareness_upper_inner .+ eqt1d.squareness_upper_outer)
            zeta = IMAS.interp1d(rho_eq, tmp).(rho_cp)
        else
            zeta = zero(rho_cp)
        end
        skappa = rmin .* IMAS.gradient(rmin, kappa) ./ kappa
        sdelta = rmin .* IMAS.gradient(rmin, delta)
        szeta = rmin .* IMAS.gradient(rmin, zeta)
    end
    a = rmin[end]
    @views q = q_profile[gridpoint_cp]

    Te_full = cp1d.electrons.temperature
    dlntedr_full = .-IMAS.calc_z(rmin, Te_full, :backward)
    @views Te = Te_full[gridpoint_cp]
    @views dlntedr = dlntedr_full[gridpoint_cp]

    ne_full = cp1d.electrons.density_thermal ./ m³_to_cm³
    dlnnedr_full = .-IMAS.calc_z(rmin, ne_full, :backward)
    @views ne = ne_full[gridpoint_cp]
    @views dlnnedr = dlnnedr_full[gridpoint_cp]

    Bt = eqt.global_quantities.vacuum_toroidal_field.b0
    buitp = IMAS.interp1d(rho_eq, GACODE.bunit(eqt1d))
    bunit = @. @views buitp(rho_cp[gridpoint_cp]) * T_to_Gauss
    input_tglf = InputTGLFs{T}([InputTGLF{T}() for k in eachindex(gridpoint_cp)])

    signb = sign(Bt)
    input_tglf.SIGN_BT = signb
    input_tglf.SIGN_IT = @. signb * sign(q)

    input_tglf.NS = length(ions) + 1 # add 1 to include electrons

    # electrons first for TGLF
    input_tglf.MASS_1 = me / md
    input_tglf.TAUS_1 = 1.0
    input_tglf.AS_1 = 1.0
    input_tglf.ZS_1 = -1.0
    input_tglf.RLNS_1 = @. a * dlnnedr
    input_tglf.RLTS_1 = @. a * dlntedr

    c_s = GACODE.c_s.(Te)
    w0 = @. -cp1d.rotation_frequency_tor_sonic
    w0p = IMAS.gradient(rmin, w0; method=:backward)
    gamma_p = @. @views -Rmaj[gridpoint_cp] * w0p[gridpoint_cp]
    gamma_e = @. @views -rmin[gridpoint_cp] / q * w0p[gridpoint_cp]
    mach = @. @views Rmaj[gridpoint_cp] * w0[gridpoint_cp] / c_s
    input_tglf.VPAR_1 = @. -input_tglf.SIGN_IT * mach
    input_tglf.VPAR_SHEAR_1 = @. -input_tglf.SIGN_IT * (a / c_s) * gamma_p
    input_tglf.VEXB_SHEAR = @. -gamma_e * (a / c_s)

    for iion in eachindex(ions)
        species = iion + 1
        Ti_full = ions[iion].temperature
        dlntidr_full = .-IMAS.calc_z(rmin, Ti_full, :backward)
        @views Ti = Ti_full[gridpoint_cp]
        @views dlntidr = dlntidr_full[gridpoint_cp]

        Zi = IMAS.avgZ(ions[iion].element[1].z_n, Ti)
        setproperty!(input_tglf, Symbol("ZS_$species"), Zi)
        setproperty!(input_tglf, Symbol("MASS_$species"), ions[iion].element[1].a .* mp ./ md)

        ni_full = ions[iion].density_thermal ./ m³_to_cm³
        dlnnidr_full = .-IMAS.calc_z(rmin, ni_full, :backward)
        @views ni = ni_full[gridpoint_cp]
        @views dlnnidr = dlnnidr_full[gridpoint_cp]

        setproperty!(input_tglf, Symbol("TAUS_$species"), Ti ./ Te)
        setproperty!(input_tglf, Symbol("AS_$species"), ni ./ ne)
        setproperty!(input_tglf, Symbol("VPAR_$species"), input_tglf.VPAR_1)
        setproperty!(input_tglf, Symbol("VPAR_SHEAR_$species"), input_tglf.VPAR_SHEAR_1)
        setproperty!(input_tglf, Symbol("RLNS_$species"), a .* dlnnidr)
        setproperty!(input_tglf, Symbol("RLTS_$species"), a .* dlntidr)
    end

    input_tglf.BETAE = @. 8π * ne * k * Te / bunit^2
    input_tglf.XNUE = @. a / c_s * sqrt(ions[1].element[1].a) * e^4 * π * ne * (24.0 - log(sqrt(ne) / Te)) / (sqrt(me) * (k * Te)^1.5)
    input_tglf.ZEFF = @views cp1d.zeff[gridpoint_cp]
    rho_s = @views GACODE.rho_s(cp1d, eqt)[gridpoint_cp]
    input_tglf.DEBYE = @. 7.43e2 * sqrt(Te / ne) / rho_s
    input_tglf.RMIN_LOC = @. @views rmin[gridpoint_cp] / a
    input_tglf.RMAJ_LOC = @. @views Rmaj[gridpoint_cp] / a
    input_tglf.ZMAJ_LOC = 0
    input_tglf.DRMINDX_LOC = 1.0

    drmaj = IMAS.gradient(rmin, Rmaj)

    input_tglf.DRMAJDX_LOC = @views drmaj[gridpoint_cp]
    input_tglf.DZMAJDX_LOC = 0.0

    input_tglf.Q_LOC = @. abs(q)

    input_tglf.KAPPA_LOC = @views kappa[gridpoint_cp]

    input_tglf.S_KAPPA_LOC = @views skappa[gridpoint_cp]
    input_tglf.DELTA_LOC = @views delta[gridpoint_cp]
    input_tglf.S_DELTA_LOC = @views sdelta[gridpoint_cp]
    input_tglf.ZETA_LOC = @views zeta[gridpoint_cp]
    input_tglf.S_ZETA_LOC = @views szeta[gridpoint_cp]

    if MXH_modes > 1
        for m in 0:n_shape
            setproperty!(input_tglf, Symbol("SHAPE_COS$m"),   @views c_mxh[m+1][gridpoint_cp])
            setproperty!(input_tglf, Symbol("SHAPE_S_COS$m"), @views dc_mxh[m+1][gridpoint_cp])
        end
        for (i, m) in enumerate(3:n_shape)
            setproperty!(input_tglf, Symbol("SHAPE_SIN$m"),   @views s_mxh[i][gridpoint_cp])
            setproperty!(input_tglf, Symbol("SHAPE_S_SIN$m"), @views ds_mxh[i][gridpoint_cp])
        end
    end

    press = cp1d.pressure
    Pa_to_dyn = 10.0

    dpdr = @views IMAS.gradient(rmin, press)[gridpoint_cp] .* Pa_to_dyn
    input_tglf.P_PRIME_LOC = @. @views abs(q) / (rmin[gridpoint_cp] / a)^2 * rmin[gridpoint_cp] / bunit^2 * dpdr

    dqdr = @views IMAS.gradient(rmin, q_profile)[gridpoint_cp]
    input_tglf.Q_PRIME_LOC = @. @views q * a^2 / rmin[gridpoint_cp] * dqdr

    # saturation rules
    input_tglf.ALPHA_ZF = -1.0 # 1 = default, -1 = low ky cutoff kypeak search
    input_tglf.USE_MHD_RULE = false
    input_tglf.NMODES = input_tglf.NS .+ 2 # capture main branches: ES each species + BPER + VPAR_SHEAR
    input_tglf.NKY = 12 # 12 is default, 16 for smoother spectrum
    input_tglf.ALPHA_QUENCH = 0 # 0 = spectral shift, 1 = quench
    input_tglf.SAT_RULE = parse(Int, split(string(sat), "sat")[end][1]) # Parsing will fail for double digit sat rules
    
    if sat == :sat2 || sat == :sat3
        input_tglf.UNITS = "CGYRO"
        input_tglf.KYGRID_MODEL = 4
        input_tglf.NBASIS_MIN = 2
        input_tglf.NBASIS_MAX = 6
        input_tglf.USE_AVE_ION_GRID = true
        input_tglf.XNU_MODEL = 3
        input_tglf.WDIA_TRAPPED = 1.0
    else
        input_tglf.UNITS = "GYRO"
        if sat == :sat1
        elseif sat == :sat1geo
            input_tglf.UNITS = "CGYRO"
        elseif sat == :sat0quench
            input_tglf.ALPHA_QUENCH = 1
        end
        input_tglf.KYGRID_MODEL = 1
        input_tglf.NBASIS_MIN = 2
        input_tglf.NBASIS_MAX = 4
        input_tglf.USE_AVE_ION_GRID = false # default is false
        input_tglf.XNU_MODEL = 2
        input_tglf.WDIA_TRAPPED = 0.0
    end

    # electrostatic/electromagnetic
    if electromagnetic
        input_tglf.USE_BPER = true
        input_tglf.USE_BPAR = true
    else
        input_tglf.USE_BPER = false
        input_tglf.USE_BPAR = false
    end

    input_tglf.ALPHA_MACH = 0.0
    return input_tglf
end

"""
    gacode_preamble() -> String

Return a shell preamble that sources the GACODE environment on NERSC Perlmutter.
Returns an empty string on non-NERSC systems where GACODE is assumed to be in PATH.
"""
function gacode_preamble()
    if get(ENV, "NERSC_HOST", "") == "perlmutter"
        gacode_root = get(ENV, "GACODE_ROOT_CPU", "")
        isempty(gacode_root) && error("Environment variable GACODE_ROOT_CPU is not set. " *
            "Set it to your CPU GACODE installation path, e.g. export GACODE_ROOT_CPU=/path/to/gacode")
        gacode_platform = "PERLMUTTER_CPU"
        return """
        export GACODE_ROOT=$gacode_root
        export GACODE_PLATFORM=$gacode_platform
        . \${GACODE_ROOT}/shared/bin/gacode_setup 2>/dev/null
        . \${GACODE_ROOT}/platform/env/env.\${GACODE_PLATFORM} 2>/dev/null
        """
    end
    return ""
end

"""
    run_tglf(input_tglf::InputTGLF)

Run TGLF starting from a InputTGLF.

Returns a `FluxSolution` structure
"""
function run_tglf(input_tglf::InputTGLF)
    folder = mktempdir()

    save(input_tglf, joinpath(folder, "input.tglf"))

    preamble = gacode_preamble()
    open(joinpath(folder, "command.sh"), "w") do io
        return write(
            io,
            """
         $preamble
         if command -v timeout &> /dev/null; then
         	(time (timeout 120 tglf -n 1 -e .)) &> command.log
         else
         	(time (tglf -n 1 -e .)) &> command.log
         fi
         """)
    end

    fluxes = try
        run(Cmd(`bash command.sh`; dir=folder))

        tmp = open(joinpath(folder, "out.tglf.gbflux"), "r") do io
            return read(io, String)
        end

        parse_out_tglf_gbflux(tmp)

    catch e
        # show last 100 lines of  chease.output
        txt = open(joinpath(folder, "command.log"), "r") do io
            return split(read(io, String), "\n")
        end
        @error "ERROR running TGLF\n...\n" * join(txt[max(1, length(txt) - 100):end], "\n")
        rethrow(e)
    end

    T = Float64
    sol = GACODE.FluxSolution{T}(
        fluxes["Q/Q_GB_elec"],
        fluxes["Q/Q_GB_ions"],
        fluxes["Gam/Gam_GB_elec"],
        fluxes["Gam/Gam_GB_all_ions"],
        fluxes["Pi/Pi_GB_ions"])
    rm(folder; force=true, recursive=true)

    return sol
end

"""
    run_tglf(input_tglfs::Vector{InputTGLF{T}}) where {T<:Real}

Run TGLF starting from a vector of InputTGLFs.

NOTE: Each run is done asyncronously (ie. in separate parallel processes)

Returns a `FluxSolution` structure
"""
function run_tglf(input_tglfs::Vector{InputTGLF{T}}) where {T<:Real}
    return collect(asyncmap(input_tglf -> run_tglf(input_tglf), input_tglfs))
end

export run_tglf

"""
    apply_presets!(input_tglf::InputTGLF)

Applies the same preset rules as TGLF Fortran for consistency:

1. SAT_RULE specific settings:
   - SAT_RULE = 2 or 3:
     * XNU_MODEL = 3
     * WDIA_TRAPPED = 1.0
     * If UNITS = "GYRO", changes to "CGYRO"
     * GEOMETRY_FLAG must not be 0
   - SAT_RULE = 1:
     * XNU_MODEL = 2
   - SAT_RULE = 0:
     * UNITS = "GYRO"
     * XNU_MODEL = 2
     * NMODES > 2 is set to 4

2. Other rules:
   - If USE_BPER is true, ALPHA_MACH = 0.0
"""
function apply_presets!(input_tglf::InputTGLF)
    # Handle SAT_RULE specific settings
    if input_tglf.SAT_RULE == 2 || input_tglf.SAT_RULE == 3
        input_tglf.XNU_MODEL = 3
        input_tglf.WDIA_TRAPPED = 1.0
        if input_tglf.UNITS == "GYRO"
            input_tglf.UNITS = "CGYRO"
        end
    elseif input_tglf.SAT_RULE == 1
        input_tglf.XNU_MODEL = 2
    elseif input_tglf.SAT_RULE == 0
        input_tglf.UNITS = "GYRO"
        input_tglf.XNU_MODEL = 2
        # Adjust NMODES for SAT_RULE=0
        if input_tglf.NMODES > 2
            input_tglf.NMODES = 4
        end
    end

    # Handle USE_BPER setting
    if input_tglf.USE_BPER
        input_tglf.ALPHA_MACH = 0.0
    end
    return input_tglf
end
