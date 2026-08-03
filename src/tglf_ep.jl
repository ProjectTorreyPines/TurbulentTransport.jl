# TGLF-EP / TJLFEP "from dd" support.
#
# This file is to TGLF-EP what `tglf.jl` is to TGLF: it builds the TGLF input
# objects from an IMAS `dd`, here augmented with the fast energetic-particle (EP)
# species and the extra EP geometry/profile constants (`extraEP`) that TJLFEP's
# stability scan consumes. `InputTGLF_EP` returns `(InputTGLFs, extraEP::Dict)`.
#
# It reuses TurbulentTransport's typed `InputTGLFs{T}` wrapper (input_tglfs.jl)
# and only needs IMAS + GACODE + TJLF.InputTGLF - no TJLFEP types - so it lives
# in TurbulentTransport without creating a dependency cycle.

# `expro_bound_deriv` / `expro_log_gradients` are DUPLICATED here from
# `TJLFEP/src/tjlfep_read_inputs.jl`. They cannot be shared from one home: the
# TJLFEP copy is needed by the always-loaded gacode-file path (which loads
# neither TurbulentTransport, IMAS, nor GACODE), while `InputTGLF_EP` below needs
# the same pure 3-point Lagrange derivative. They are tiny and pure-numeric.

"""
    expro_bound_deriv(f, r) -> Vector{Float64}

Match GACODE `expro_util.f90::bound_deriv` (Lagrange derivative on uniform radial grid).
"""
function expro_bound_deriv(f::AbstractVector{<:Real}, r::AbstractVector{<:Real})
    n = length(f)
    n == length(r) || error("expro_bound_deriv: length(f)=$(length(f)) != length(r)=$(length(r))")
    df = Vector{Float64}(undef, n)
    f = Float64.(f)
    r = Float64.(r)
    for i in 1:n
        if i == 1
            ra, r1, r2, r3 = r[1], r[1], r[2], r[3]
            f1, f2, f3 = f[1], f[2], f[3]
        elseif i == n
            ra, r1, r2, r3 = r[n], r[n - 2], r[n - 1], r[n]
            f1, f2, f3 = f[n - 2], f[n - 1], f[n]
        else
            ra, r1, r2, r3 = r[i], r[i - 1], r[i], r[i + 1]
            f1, f2, f3 = f[i - 1], f[i], f[i + 1]
        end
        df[i] = ((ra - r1) + (ra - r2)) / (r3 - r1) / (r3 - r2) * f3 +
                ((ra - r1) + (ra - r3)) / (r2 - r1) / (r2 - r3) * f2 +
                ((ra - r2) + (ra - r3)) / (r1 - r2) / (r1 - r3) * f1
    end
    return df
end

"""`dlnnidr = -d(ln n)/dr` and `dlntidr = -d(ln T)/dr` as in `expro_compute_derived`."""
function expro_log_gradients(ni::AbstractVector, ti::AbstractVector, rmin::AbstractVector)
    nr = length(ni)
    length(ti) == nr == length(rmin) || error("ni/ti/rmin length mismatch")
    eps_n = 1.0e-30
    dlnnidr = expro_bound_deriv(-log.(max.(ni, eps_n)), rmin)
    dlntidr = expro_bound_deriv(-log.(max.(ti, eps_n)), rmin)
    return dlnnidr, dlntidr
end

"""
    InputTGLFEP(dd::IMAS.dd, rho::AbstractVector{Float64}, sat::Symbol=:sat0, electromagnetic::Bool=false, lump_ions::Bool=false; is_ep::Int=1)

Evaluate TGLF-EP input parameters at given radii, appending the fast energetic-particle species.
Returns `(input_tglf::InputTGLFs, extraEP::Dict)`.
"""
function InputTGLFEP(dd::IMAS.dd, rho::AbstractVector{Float64}, sat::Symbol=:sat0, electromagnetic::Bool=false, lump_ions::Bool=false; is_ep::Int=1)
    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]
    gridpoint_cp = [argmin(abs.(cp1d.grid.rho_tor_norm .- ρ)) for ρ in rho]
    return InputTGLF_EP(eqt, cp1d, gridpoint_cp, sat, electromagnetic, lump_ions; is_ep)
end

"""
    InputTGLFEP(dd::IMAS.dd, gridpoint_cp::AbstractVector{Int}, sat::Symbol=:sat0, electromagnetic::Bool=false, lump_ions::Bool=false; is_ep::Int=1)

Evaluate TGLF-EP input parameters at given core profiles grid indexes.
Returns `(input_tglf::InputTGLFs, extraEP::Dict)`.
"""
function InputTGLFEP(dd::IMAS.dd, gridpoint_cp::AbstractVector{Int}, sat::Symbol=:sat0, electromagnetic::Bool=false, lump_ions::Bool=false; is_ep::Int=1)
    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]
    return InputTGLF_EP(eqt, cp1d, gridpoint_cp, sat, electromagnetic, lump_ions; is_ep)
end

function InputTGLF_EP(
    eqt::IMAS.equilibrium__time_slice{T},
    cp1d::IMAS.core_profiles__profiles_1d{T},
    gridpoint_cp::AbstractVector{Int},
    sat::Symbol,
    electromagnetic::Bool,
    lump_ions::Bool;
    is_ep::Int=1) where {T<:Real}

    e = IMAS.cgs.e # statcoul
    k = IMAS.cgs.k # erg/eV
    mp = IMAS.cgs.mp # g
    me = IMAS.cgs.me # g
    md = IMAS.cgs.md # g
    m_to_cm = IMAS.cgs.m_to_cm
    m³_to_cm³ = IMAS.cgs.m³_to_cm³
    T_to_Gauss = IMAS.cgs.T_to_Gauss

    eq1d = eqt.profiles_1d

    if lump_ions
        ions = IMAS.lump_ions_as_bulk_and_impurity(cp1d)
    else
        ions = cp1d.ion
    end

    # Rmaj = IMAS.interp1d(eq1d.rho_tor_norm, m_to_cm * 0.5 * (eq1d.r_outboard .+ eq1d.r_inboard)).(cp1d.grid.rho_tor_norm)
    Rmaj = IMAS.interp1d(eq1d.rho_tor_norm, 0.5 * (eq1d.r_outboard .+ eq1d.r_inboard)).(cp1d.grid.rho_tor_norm)


    rmin = GACODE.r_min_core_profiles(eqt.profiles_1d, cp1d.grid.rho_tor_norm)
    rmin = rmin / m_to_cm

    # Match the input.gacode preprocessing path (tjlfep_read_inputs.jl / GACODE expro_util.f90):
    # compute ALL gradients with the GACODE 3-point Lagrange derivative `expro_bound_deriv`
    # on rmin (with the same axis fix rmin+1e-6 and -log / eps_n conventions) instead of
    # IMAS.calc_z / IMAS.gradient. This way a dd converted from an input.gacode reproduces the
    # gacode-file TJLFEP gradients (RLNS/RLTS and geometry s-factors).
    rmin_work = rmin .+ 1.0e-6
    eps_n = 1.0e-30

    q_profile = IMAS.interp1d(eq1d.rho_tor_norm, eq1d.q).(cp1d.grid.rho_tor_norm)

    if !ismissing(eq1d, :elongation)
        kappa = IMAS.interp1d(eq1d.rho_tor_norm, eq1d.elongation).(cp1d.grid.rho_tor_norm)
    else
        kappa = zero(cp1d.grid.rho_tor_norm)
    end

    if !ismissing(eq1d, :triangularity_lower) && !ismissing(eq1d, :triangularity_upper)
        delta = IMAS.interp1d(eq1d.rho_tor_norm, 0.5 * (eq1d.triangularity_lower + eq1d.triangularity_upper)).(cp1d.grid.rho_tor_norm)
    else
        delta = zero(cp1d.grid.rho_tor_norm)
    end

    if !ismissing(eq1d, :squareness_lower_inner) && !ismissing(eq1d, :squareness_lower_outer) && !ismissing(eq1d, :squareness_upper_inner) &&
       !ismissing(eq1d, :squareness_upper_outer)
        tmp = 0.25 * (eq1d.squareness_lower_inner .+ eq1d.squareness_lower_outer .+ eq1d.squareness_upper_inner .+ eq1d.squareness_upper_outer)
        zeta = IMAS.interp1d(eq1d.rho_tor_norm, tmp).(cp1d.grid.rho_tor_norm)
    else
        zeta = zero(cp1d.grid.rho_tor_norm)
    end

    a = rmin[end]
    q = q_profile[gridpoint_cp]

    Bt = eqt.global_quantities.vacuum_toroidal_field.b0
    bunit = IMAS.interp1d(eq1d.rho_tor_norm, GACODE.bunit(eqt) .* T_to_Gauss).(cp1d.grid.rho_tor_norm)[gridpoint_cp]
    input_tglf = InputTGLFs{T}([InputTGLF{T}() for k in eachindex(gridpoint_cp)])

    # signb = sign(Bt)  # original: IMAS convention (bcentr may be stored as magnitude+positive)
    # signq = sign.(q)  # original: signed q from IMAS
    # Fortran: EXPRO_ctrl_signb = -1.0 flips B sign; EXPRO_ctrl_signq = 1.0 forces q positive
    # → sign_bt = -1, sign_it = sign_bt * 1 = -1 for D3D (negative BT, negative IP)
    signb = -sign(Bt)
    signq = sign.(abs.(q))  # EXPRO_ctrl_signq = 1 → q treated as positive
    # input_tglf.SIGN_BT = signb
    input_tglf.SIGN_BT = -1.0
    # input_tglf.SIGN_IT = signb .* signq
    input_tglf.SIGN_IT = -1.0

    is_ep = is_ep + 1 # fast EP included in indexing here
    ns = length(ions) + 2  # electrons + all thermal ions + fast EP appended at last slot
    input_tglf.NS = ns

    c_s = GACODE.c_s(cp1d)[gridpoint_cp]
    w0 = -1 * cp1d.rotation_frequency_tor_sonic
    w0p = expro_bound_deriv(w0, rmin_work)
    gamma_p = -Rmaj[gridpoint_cp] .* w0p[gridpoint_cp]
    gamma_e = -rmin[gridpoint_cp] ./ q .* w0p[gridpoint_cp]
    mach = Rmaj[gridpoint_cp] .* w0[gridpoint_cp] ./ c_s
    input_tglf.VPAR_1 = -input_tglf.SIGN_IT .* mach
    input_tglf.VPAR_SHEAR_1 = -input_tglf.SIGN_IT .* (a ./ c_s) .* gamma_p
    input_tglf.VEXB_SHEAR = -gamma_e .* (a ./ c_s)

    # ===== EP Data: Calculate for ALL radii (TJLFEP needs full profiles) =====
    # Get FULL radial grid (not just selected points)
    fac = 1e19

    nr_full = length(cp1d.grid.rho_tor_norm)
    ns_full = ns
    
    # Preallocate EP arrays for all radii and species (6 max species for compatibility)
    ep_species_data = Dict{String, Vector{Float64}}()
    
    # Full radial profiles for thermal electrons
    Te_full = cp1d.electrons.temperature ./ 1000 # eV -> keV
    dlntedr_full = expro_bound_deriv(-log.(max.(Te_full, eps_n)), rmin_work)
    
    # ne_full = cp1d.electrons.density_thermal ./ m³_to_cm³
    ne_full = cp1d.electrons.density_thermal ./ fac
    dlnnedr_full = expro_bound_deriv(-log.(max.(ne_full, eps_n)), rmin_work)
    
    # Store full radial data for TJLFEP
    ep_species_data["DENS_1"] = ne_full
    ep_species_data["TEMP_1"] = Te_full
    ep_species_data["DLNNDR_1"] = dlnnedr_full
    ep_species_data["DLNTDR_1"] = dlntedr_full
    ep_species_data["ZS_1"] = fill(-1.0, nr_full)

    # Extract values at selected grid points for TGLF
    Te = Te_full[gridpoint_cp]
    dlntedr = dlntedr_full[gridpoint_cp]
    ne = ne_full[gridpoint_cp]
    dlnnedr = dlnnedr_full[gridpoint_cp]

    Ze = -1.0

    setproperty!(input_tglf, Symbol("ZS_1"), Ze)
    setproperty!(input_tglf, Symbol("MASS_1"), me / md)
    setproperty!(input_tglf, Symbol("TAUS_1"), Te ./ Te)
    setproperty!(input_tglf, Symbol("AS_1"), ne ./ ne)
    setproperty!(input_tglf, Symbol("VPAR_1"), input_tglf.VPAR_1)
    setproperty!(input_tglf, Symbol("VPAR_SHEAR_1"), input_tglf.VPAR_SHEAR_1)
    setproperty!(input_tglf, Symbol("RLNS_1"), a * dlnnedr)
    setproperty!(input_tglf, Symbol("RLTS_1"), a * dlntedr)

    # Ion species (2 through NS) - full radial profiles
    for iion in eachindex(ions)
        species = iion + 1

        z_n = ions[iion].element[1].z_n

        # thermal
        Ti_full = ions[iion].temperature ./ 1000 # eV -> keV
        dlntidr_full = expro_bound_deriv(-log.(max.(Ti_full, eps_n)), rmin_work)
        
        # ni_full = ions[iion].density_thermal ./ m³_to_cm³
        ni_full = ions[iion].density_thermal ./ fac
        dlnnidr_full = expro_bound_deriv(-log.(max.(ni_full, eps_n)), rmin_work)
        
        # Store full radial data for TJLFEP
        ep_species_data["DENS_$species"] = ni_full
        ep_species_data["TEMP_$species"] = Ti_full
        ep_species_data["DLNNDR_$species"] = dlnnidr_full
        ep_species_data["DLNTDR_$species"] = dlntidr_full
        ep_species_data["ZS_$species"] = IMAS.avgZ(Float64(z_n), Ti_full)

        # Extract values at selected grid points for TGLF
        Ti = Ti_full[gridpoint_cp]
        dlntidr = dlntidr_full[gridpoint_cp]
        ni = ni_full[gridpoint_cp]
        dlnnidr = dlnnidr_full[gridpoint_cp]

        Zi = IMAS.avgZ(Float64(ions[iion].element[1].z_n), Ti)

        setproperty!(input_tglf, Symbol("ZS_$species"), Zi)
        setproperty!(input_tglf, Symbol("MASS_$species"), (ions[iion].element[1].a)[1] * mp / md)
        setproperty!(input_tglf, Symbol("TAUS_$species"), Ti ./ Te)
        setproperty!(input_tglf, Symbol("AS_$species"), ni ./ ne)
        setproperty!(input_tglf, Symbol("VPAR_$species"), input_tglf.VPAR_1)
        setproperty!(input_tglf, Symbol("VPAR_SHEAR_$species"), input_tglf.VPAR_SHEAR_1)
        setproperty!(input_tglf, Symbol("RLNS_$species"), a * dlnnidr)
        setproperty!(input_tglf, Symbol("RLTS_$species"), a * dlntidr)

    end

    ep_slot = ns  # EP is always the last slot (length(ions) + 2)

    if is_ep == 1
        # Full radial profiles for fast electrons
        # ne_full = cp1d.electrons.density_fast ./ m³_to_cm³
        ne_full = cp1d.electrons.density_fast ./ fac
        fast_zero = ne_full .== 0
        dlnnedr_full = ifelse.(fast_zero, 0.0, expro_bound_deriv(-log.(max.(ne_full, eps_n)), rmin_work))

        # Pressure is in Pa (J/m³); T [eV] = P / (n [m⁻³] × 1.602e-19 J/eV). Note: e = IMAS.cgs.e is CGS, do NOT use here.
        Te_full = (cp1d.electrons.pressure_fast_parallel ./ 3 + 2 * cp1d.electrons.pressure_fast_perpendicular ./ 3) ./ max.(ne_full .* fac, eps(Float64)) ./ 1.602e-19
        Te_full[fast_zero] .= 0.0
        Te_full = Te_full ./ 1000 # eV -> keV
        dlntedr_full = ifelse.(fast_zero, 0.0, expro_bound_deriv(-log.(max.(Te_full, eps_n)), rmin_work))
        dlntedr_full[fast_zero] .= 0.0
        ep_species_data["DENS_$ep_slot"] = ne_full
        ep_species_data["TEMP_$ep_slot"] = Te_full
        ep_species_data["DLNNDR_$ep_slot"] = dlnnedr_full
        ep_species_data["DLNTDR_$ep_slot"] = dlntedr_full
        ep_species_data["ZS_$ep_slot"] = fill(-1.0, nr_full)

        Te_fast = Te_full[gridpoint_cp]
        dlntedr_fast = dlntedr_full[gridpoint_cp]
        nef = ne_full[gridpoint_cp]
        dlnnedr_fast = dlnnedr_full[gridpoint_cp]

        setproperty!(input_tglf, Symbol("ZS_$ep_slot"), -1.0)
        setproperty!(input_tglf, Symbol("MASS_$ep_slot"), me / md)
        setproperty!(input_tglf, Symbol("TAUS_$ep_slot"), Te_fast ./ Te)
        setproperty!(input_tglf, Symbol("AS_$ep_slot"), nef ./ ne)
        setproperty!(input_tglf, Symbol("VPAR_$ep_slot"), input_tglf.VPAR_1)
        setproperty!(input_tglf, Symbol("VPAR_SHEAR_$ep_slot"), input_tglf.VPAR_SHEAR_1)
        setproperty!(input_tglf, Symbol("RLNS_$ep_slot"), a * dlnnedr_fast)
        setproperty!(input_tglf, Symbol("RLTS_$ep_slot"), a * dlntedr_fast)
    else
        ep_ion = ions[is_ep - 1]

        ni_full = ep_ion.density_fast ./ fac
        fast_zero = ni_full .== 0
        dlnnidr_full = ifelse.(fast_zero, 0.0, expro_bound_deriv(-log.(max.(ni_full, eps_n)), rmin_work))

        Ti_full = (ep_ion.pressure_fast_parallel .+ 2 .* ep_ion.pressure_fast_perpendicular) ./ max.(ni_full .* fac, eps(Float64)) ./ 1.602e-19
        Ti_full[fast_zero] .= 0.0
        Ti_full = Ti_full ./ 1000 # eV -> keV

        # Mirror GACODE's own low-density fix (inputgacode.jl:487): replace Ti where density
        # is negligibly small (< 1e-6 × mean) with the mean Ti over valid points, so calc_z
        # sees a smooth profile instead of P/(~0·e) blow-up.
        ni_full_m3 = ep_ion.density_fast  # m⁻³
        navg = sum(ni_full_m3) / length(ni_full_m3)
        low_dens = ni_full_m3 .< 1e-6 .* navg
        valid = .!low_dens .& .!fast_zero
        Ti_full[low_dens] .= sum(Ti_full[valid]) / max(1, count(valid))

        dlntidr_full = ifelse.(fast_zero, 0.0, expro_bound_deriv(-log.(max.(Ti_full, eps_n)), rmin_work))
        dlntidr_full[fast_zero] .= 0.0
        dlntidr_full[low_dens] .= 0.0  # gradient unreliable in mean-filled region

        # Match Fortran floor on EP density gradient (TGLFEP_read_EXPRO.f90: max(...,1.0)).
        # Fortran floors EXPRO_dlnnidr at 1.0 m⁻¹, giving RLNS floor = a_meters ≈ 0.6-0.7.
        # Using 1.0 here (not 1.0/a) replicates that: RLNS = a * 1.0 = a.
        dlnnidr_full .= max.(dlnnidr_full, 1.0)

        ep_species_data["DENS_$ep_slot"] = ni_full
        ep_species_data["TEMP_$ep_slot"] = Ti_full
        ep_species_data["DLNNDR_$ep_slot"] = dlnnidr_full
        ep_species_data["DLNTDR_$ep_slot"] = dlntidr_full
        ep_species_data["ZS_$ep_slot"] = max.(IMAS.avgZ(Float64(ep_ion.element[1].z_n), Ti_full), 1.0)

        Ti_ep = Ti_full[gridpoint_cp]
        dlntidr_ep = dlntidr_full[gridpoint_cp]
        ni_ep = ni_full[gridpoint_cp]
        dlnnidr_ep = dlnnidr_full[gridpoint_cp]

        Zi_ep = max.(IMAS.avgZ(Float64(ep_ion.element[1].z_n), Ti_ep), 1.0)

        setproperty!(input_tglf, Symbol("ZS_$ep_slot"), Zi_ep)
        setproperty!(input_tglf, Symbol("MASS_$ep_slot"), (ep_ion.element[1].a)[1] * mp / md)
        setproperty!(input_tglf, Symbol("TAUS_$ep_slot"), Ti_ep ./ Te)
        setproperty!(input_tglf, Symbol("AS_$ep_slot"), ni_ep ./ ne)
        setproperty!(input_tglf, Symbol("VPAR_$ep_slot"), input_tglf.VPAR_1)
        setproperty!(input_tglf, Symbol("VPAR_SHEAR_$ep_slot"), input_tglf.VPAR_SHEAR_1)
        setproperty!(input_tglf, Symbol("RLNS_$ep_slot"), a * dlnnidr_ep)
        setproperty!(input_tglf, Symbol("RLTS_$ep_slot"), a * dlntidr_ep)
    end

    # cgs formulas (cf. tglf.jl and betae_full below): this function's working
    # units are Te [keV], ne [1e19 m^-3], a [m] — convert before use.
    ne_cgs = ne .* 1e13            # 1e19 m^-3 -> cm^-3
    Te_eV = Te .* 1e3              # keV -> eV
    input_tglf.BETAE = 8.0 * pi .* ne_cgs .* k .* Te_eV ./ bunit .^ 2
    loglam = 24.0 .- log.(sqrt.(ne_cgs) ./ Te_eV)
    input_tglf.XNUE = (a * 100) ./ c_s * sqrt.(ions[1].element[1].a) .* e^4 .* pi .* ne_cgs .* loglam ./ (sqrt.(me) .* (k .* Te_eV) .^ 1.5)
    input_tglf.ZEFF = cp1d.zeff[gridpoint_cp]
    rho_s = GACODE.rho_s(cp1d, eqt)[gridpoint_cp]
    input_tglf.DEBYE = 7.43e2 * sqrt.(Te_eV ./ ne_cgs) ./ rho_s
    input_tglf.RMIN_LOC = rmin[gridpoint_cp] ./ a
    input_tglf.RMAJ_LOC = Rmaj[gridpoint_cp] ./ a
    input_tglf.ZMAJ_LOC = 0
    input_tglf.DRMINDX_LOC = 1.0

    drmaj = expro_bound_deriv(Rmaj, rmin_work)

    input_tglf.DRMAJDX_LOC = drmaj[gridpoint_cp]
    input_tglf.DZMAJDX_LOC = 0.0

    input_tglf.Q_LOC = abs.(q)

    # omegaGAM = (1.0 ./ input_tglf.RMAJ_LOC) .* sqrt.(1.0 .+ input_tglf.TAUS_2) ./ (1.0 .+ 1.0 ./ (2.0 .* q))
    # OMEGA_TAE = sqrt.(2.0 ./ input_tglf.BETAE) ./ 2.0 ./ q ./ input_tglf.RMAJ_LOC
    RHO_STAR = rho_s
    CS = c_s
    RMIN = rmin
    NS = input_tglf.NS[1]
    

    # ===== Calculate full radial profiles for EP (before extracting selected points) =====
    # Sound speed: reproduce the input.gacode path EXACTLY (expro_util.f90 CGS formula, Te in keV).
    # GACODE.c_s uses the deuterium mass md (→ factor sqrt(2) larger); the gacode-file path
    # (tjlfep_read_inputs.jl) divides by 2*md and works in m/s. Match the file path so a dd
    # converted from an input.gacode reproduces its cs / gammaE / gammap.
    # Use thermal electron profiles (not fast electrons, which overwrote ne_full/Te_full above)
    ne_thermal_full = ep_species_data["DENS_1"]
    Te_thermal_full = ep_species_data["TEMP_1"]
    k_erg = 1.6022e-12
    mp_g  = 2.0 * 1.6726e-24
    cs_mps_full = sqrt.(k_erg .* (1e3 .* Te_thermal_full) ./ (2.0 * mp_g)) ./ 1e2   # m/s
    c_s_full = cs_mps_full .* 100.0                                                 # cm/s (readEXPRO / F_REAL CS convention)
    rho_s_full = GACODE.rho_s(cp1d, eqt)

    # OMEGA_TAE and omegaGAM need full radial betae
    # ne_thermal_full is in 10^19/m^3 → need cm^-3 (×1e13); Te_thermal_full is in keV → need eV (×1e3); combined factor 1e16
    betae_full = 8.0 * pi .* (ne_thermal_full .* 1e13) .* k .* (Te_thermal_full .* 1e3) ./ (IMAS.interp1d(eq1d.rho_tor_norm, GACODE.bunit(eqt) .* T_to_Gauss).(cp1d.grid.rho_tor_norm)) .^ 2
    
    # Full radial rotation profiles for gamma_e and gamma_p
    w0_full = -cp1d.rotation_frequency_tor_sonic
    w0p_full = expro_bound_deriv(w0_full, rmin_work)
    Rmaj_full = IMAS.interp1d(eq1d.rho_tor_norm, 0.5 * (eq1d.r_outboard .+ eq1d.r_inboard)).(cp1d.grid.rho_tor_norm)
    q_full = abs.(q_profile)  # positive q for shear / OMEGA_TAE / omegaGAM (EXPRO_ctrl_signq=1)
    # gammaE/gammap normalized to c_s/a (dimensionless). Match the gacode-file path exactly
    # (tjlfep_read_inputs.jl): gamma_e_phys = -rmin/q·w0p with SIGNED q, then ×(a/cs_mps) in m/(m/s).
    gamma_e_full = -rmin ./ q_profile .* w0p_full .* (a ./ cs_mps_full)
    gamma_p_full = -Rmaj_full .* w0p_full .* (a ./ cs_mps_full)
    
    # Calculate OMEGA_TAE and omegaGAM for full radial grid
    OMEGA_TAE_full = sqrt.(2.0 ./ betae_full) ./ (2.0 .* q_full .* Rmaj_full ./ a)
    
    # For omegaGAM, we need TAUS_2 at all radii
    # Ti_ion1_full is in eV (IMAS); Te_thermal_full is in keV → divide by 1000 to get same units
    if length(ions) >= 1
        Ti_ion1_full = ions[1].temperature
        # Fortran uses EXPRO_ctrl_signq=1 → q > 0 always → denominator 1+1/(2q) > 1.
        # Use abs(q) (= q_full) to match that convention.
        omegaGAM_full = (a ./ Rmaj_full) .* sqrt.(1.0 .+ Ti_ion1_full ./ (Te_thermal_full .* 1e3)) ./ (1.0 .+ 1.0 ./ (2.0 .* q_full))
    else
        omegaGAM_full = zeros(nr_full)
    end

    input_tglf.KAPPA_LOC = kappa[gridpoint_cp]

    skappa = (rmin_work ./ kappa) .* expro_bound_deriv(kappa, rmin_work)
    sdelta = rmin_work .* expro_bound_deriv(delta, rmin_work)
    szeta = rmin_work .* expro_bound_deriv(zeta, rmin_work)

    input_tglf.S_KAPPA_LOC = skappa[gridpoint_cp]
    input_tglf.DELTA_LOC = delta[gridpoint_cp]
    input_tglf.S_DELTA_LOC = sdelta[gridpoint_cp]
    input_tglf.ZETA_LOC = zeta[gridpoint_cp]
    input_tglf.S_ZETA_LOC = szeta[gridpoint_cp]

    # press = cp1d.pressure_thermal  # original: thermal pressure only
    # Fortran: EXPRO_ptot includes fast-ion pressure; use total pressure to match
    press = cp1d.pressure
    Pa_to_dyn = 10.0

    dpdr = expro_bound_deriv(press .* Pa_to_dyn, rmin_work)[gridpoint_cp]
    input_tglf.P_PRIME_LOC = abs.(q) ./ (rmin[gridpoint_cp] ./ a) .^ 2 .* rmin[gridpoint_cp] ./ bunit .^ 2 .* dpdr

    # Magnetic shear as in GACODE expro_util: s = rmin * d(ln|q|)/dr (log form), not (rmin/q)·dq/dr.
    shear_loc = rmin_work .* expro_bound_deriv(log.(abs.(q_profile)), rmin_work)
    s = shear_loc[gridpoint_cp]
    input_tglf.Q_PRIME_LOC = q .^ 2 .* a .^ 2 ./ rmin[gridpoint_cp] .^ 2 .* s

    # saturation rules
    input_tglf.ALPHA_ZF = 1.0 # 1 = default, -1 = low ky cutoff kypeak search
    input_tglf.USE_MHD_RULE = false
    input_tglf.NMODES = input_tglf.NS .+ 2 # capture main branches: ES each species + BPER + VPAR_SHEAR
    input_tglf.NKY = 12 # 12 is default, 16 for smoother spectrum
    input_tglf.ALPHA_QUENCH = 0 # 0 = spectral shift, 1 = quench
    input_tglf.SAT_RULE = parse(Int,split(string(sat),"sat")[end])
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
        input_tglf.USE_BPAR = false # TGLF does not have enough moments to resolve BPAR flutter
    else
        input_tglf.USE_BPER = false
        input_tglf.USE_BPAR = false
    end

    input_tglf.ALPHA_MACH = 0.0

    # ===== Compute full-radial geometry for TJLFEP profile struct =====
    bunit_full = IMAS.interp1d(eq1d.rho_tor_norm, GACODE.bunit(eqt) .* T_to_Gauss).(cp1d.grid.rho_tor_norm)
    rmin_safe = max.(rmin, rmin[2])  # avoid division by zero at magnetic axis (rmin[1] ≈ 0)
    shear_full = rmin_work .* expro_bound_deriv(log.(abs.(q_profile)), rmin_work)
    q_prime_full = abs.(q_profile) .^ 2 .* a .^ 2 ./ rmin_safe .^ 2 .* shear_full
    dpdr_full = expro_bound_deriv(press .* Pa_to_dyn, rmin_work)
    p_prime_full = abs.(q_profile) ./ (rmin_safe ./ a) .^ 2 .* rmin_safe ./ bunit_full .^ 2 .* dpdr_full
    ir_rep = max(1, nr_full ÷ 2)

    # Build full-radial ZS matrix (nr_full × NS) from ep_species_data
    zs_full = hcat([ep_species_data["ZS_$j"] for j in 1:NS]...)
    mass_vec = vcat([getproperty(input_tglf, Symbol("MASS_$j"))[1] for j in 1:NS])

    grid = collect(cp1d.grid.rho_tor_norm)

    # ===== Build extraEP dictionary with FULL radial profiles for TJLFEP =====
    extraEP = merge(ep_species_data, Dict(
        "gammaE" => gamma_e_full,
        "gammap" => gamma_p_full,
        "omegaGAM" => omegaGAM_full,
        "OMEGA_TAE" => OMEGA_TAE_full,
        "RHO_STAR" => rho_s_full ./ (a .* m_to_cm),  # rho_s_full in cm, a in m → rho_s_cm/a_cm (dimensionless ρ_s/a), matching Fortran: rho_star = EXPRO_rhos/a_cm
        "CS" => c_s_full,
        "RMIN" => rmin,
        "NR" => nr_full,
        "NS" => ns_full,
        # "SIGN_BT" => Float64(signb),
        "SIGN_BT" => -1.0,
        # "SIGN_IT" => Float64(signb * sign(q_profile[ir_rep])),  # original: uses signed q
        # Fortran: EXPRO_ctrl_signq = 1 forces q positive → SIGN_IT = SIGN_BT
        # "SIGN_IT" => Float64(signb),
        "SIGN_IT" => -1.0,
        "RMAJ" => Rmaj ./ a,
        "SHIFT" => drmaj,
        "Q" => abs.(q_profile),
        "SHEAR" => shear_full,
        "Q_PRIME" => q_prime_full,
        "P_PRIME" => p_prime_full,
        "KAPPA" => kappa,
        "S_KAPPA" => skappa,
        "DELTA" => delta,
        "S_DELTA" => sdelta,
        "ZETA" => zeta,
        "S_ZETA" => szeta,
        "BETAE" => betae_full,
        "ZEFF" => Float64.(cp1d.zeff),
        "B_UNIT" => bunit_full,
        "N_ION" => length(ions),
        "EP_SLOT" => ep_slot,
        "ZS" => zs_full,
        "MASS" => mass_vec,
        "grid" => grid
    ))

    return input_tglf, extraEP
end

# ===================================================================================
# run_tjlfep: SLURM submit/poll driver for the FUSE-dd TJLFEP scan.
#
# Mirrors qlgyro.jl's run_qlgyro: render a batch script, sbatch it, parse the job id,
# and persist a small state object that refresh_tjlfep!/tjlfep_status poll via sacct.
# The submission is a pure shell-out (no FUSE/TJLFEP dependency added here); the heavy
# lifting runs inside the rendered job (master builds the dd via FUSE, addprocs the GPU
# workers, ActorTJLFEP -> TJLFEP.runTHD over the rho scan).
#
# Defaults reproduce the validated 5-node / 20-GPU + premium/1h layout.
# ===================================================================================

const TJLFEP_CFS_SYSIMAGE = "/global/cfs/cdirs/m3739/TJLFEP/TJLFEP_gpu_generic_sysimage.so"

"""
    TGLFEPRunState

State for a submitted `run_tjlfep` SLURM job, polled by [`refresh_tjlfep!`](@ref) /
[`tjlfep_status`](@ref) and consumed by [`load_tjlfep_results`](@ref).
"""
mutable struct TGLFEPRunState
    basedir::String
    case::Symbol
    n_scan::Int
    n_basis::Int
    ngrid::Int
    gpu::Bool
    alpha_solver::Symbol
    nodes::Int
    inner::Symbol
    mps_team::Int
    job_id::String
    build_job_id::String
    batchfile::String
    logfile::String
    sysimage::String
    status::Symbol
end

_git_sha(dir::AbstractString) =
    try
        strip(read(`git -C $dir rev-parse HEAD`, String))
    catch
        ""
    end

"""
    _tjlfep_sysimage_ok(sysimage, tjlfep_root, tjlf_root) -> Bool

True iff `sysimage` exists and its `<sysimage>.sha` sidecar matches the current
TJLFEP + TJLF git HEADs (written by `batch_build_gpu_sysimage_generic.sh`).
"""
function _tjlfep_sysimage_ok(sysimage::AbstractString, tjlfep_root::AbstractString, tjlf_root::AbstractString)
    (isempty(sysimage) || !isfile(sysimage)) && return false
    sidecar = sysimage * ".sha"
    isfile(sidecar) || return false
    want = "TJLFEP=$(_git_sha(tjlfep_root))\nTJLF=$(_git_sha(tjlf_root))"
    return strip(read(sidecar, String)) == strip(want)
end

const _TJLFEP_BATCH_TEMPLATE = raw"""#!/bin/bash -l
# Auto-generated by TurbulentTransport.run_tjlfep -- FUSE-dd TJLFEP scan on the
# @NODES@-node / @NTASKS@-GPU layout. Master builds the @CASE@ dd via FUSE and
# addprocs(SlurmManager()) the GPU workers; ActorTJLFEP -> TJLFEP.runTHD scans the radii.
#SBATCH -A @ACCOUNT@
#SBATCH -q @QOS@
#SBATCH -N @NODES@
#SBATCH -n @NTASKS@
#SBATCH -t @WALLTIME@
#SBATCH -C gpu
#SBATCH -J @JOBNAME@
#SBATCH -o @BASEDIR@/tjlfep_%j.out
#SBATCH -e @BASEDIR@/tjlfep_%j.err
#SBATCH --ntasks-per-node=@TPN@
#SBATCH --cpus-per-task=@CPT@
#SBATCH --gpus-per-node=@GPN@

set -uo pipefail

module load cudatoolkit/12.9
module load julia/1.11.7
export JULIA_DEPOT_PATH="${PSCRATCH}/.julia${JULIA_DEPOT_PATH:+:${JULIA_DEPOT_PATH}}"

export TJLFEP_ROOT="@TJLFEP_ROOT@"
export FUSE_ROOT="@FUSE_ROOT@"

export CASE="@CASE@"
export SCAN_N="@SCAN_N@"
export N_BASIS="@N_BASIS@"
export NGRID="@NGRID@"
export ALPHA_SOLVER="@ALPHA_SOLVER@"
export SOLVER="@SOLVER@"
export AD_EXTEND_MODE="@EXTEND_MODE@"
export AD_WIDE_KDESC="@WIDE_KDESC@"
export AD_FAITHFUL_CONFIRM="@FAITHFUL_CONFIRM@"
export INNER="@INNER@"
export MPS_TEAM="@MPS_TEAM@"
export TJLFEP_OUT_DIR="@BASEDIR@"
export JULIA_WORKER_THREADS="@WORKER_THREADS@"
export JULIA_CUDA_USE_COMPAT=false

# Shared GPU sysimage: master AND workers load the SAME baked image so the serialized
# pmap-closure gensyms match (no UndefVarError #NNN#NNN, no per-worker eigensolve JIT).
# Empty => JIT fallback (master + workers both JIT from dev source, closure-safe).
export TJLFEP_GPU_SYSIMAGE="@SYSIMAGE@"

# --- NVIDIA MPS: one control daemon per node so an MPS team of workers can share a GPU ---
if [ "${INNER}" = "mps_team" ]; then
  export CUDA_MPS_PIPE_DIRECTORY="/tmp/nvidia-mps.${SLURM_JOB_ID}"
  export CUDA_MPS_LOG_DIRECTORY="/tmp/nvidia-log.${SLURM_JOB_ID}"
  srun --ntasks="${SLURM_NNODES}" --ntasks-per-node=1 bash -lc '
    mkdir -p "${CUDA_MPS_PIPE_DIRECTORY}" "${CUDA_MPS_LOG_DIRECTORY}"
    nvidia-cuda-mps-control -d || true
  ' || true
fi

cd "${TJLFEP_ROOT}/build"
echo "=== run_tjlfep case=${CASE} job=${SLURM_JOB_ID} nodes=${SLURM_NNODES} tasks=${SLURM_NTASKS} ==="
echo "N_BASIS=${N_BASIS} SCAN_N=${SCAN_N} NGRID=${NGRID} ALPHA_SOLVER=${ALPHA_SOLVER} SOLVER=${SOLVER} EXTEND_MODE=${AD_EXTEND_MODE} WIDE_KDESC=${AD_WIDE_KDESC} FAITHFUL_CONFIRM=${AD_FAITHFUL_CONFIRM} INNER=${INNER} MPS_TEAM=${MPS_TEAM}"
nvidia-smi -L 2>/dev/null | head -4 || true

SYSFLAG=""
if [ -n "${TJLFEP_GPU_SYSIMAGE}" ] && [ -f "${TJLFEP_GPU_SYSIMAGE}" ]; then
  SYSFLAG="--sysimage=${TJLFEP_GPU_SYSIMAGE}"
  echo "master sysimage = ${TJLFEP_GPU_SYSIMAGE}"
else
  echo "master sysimage = <none, JIT from dev source>"
fi

# Master runs directly (NOT under srun) so addprocs(SlurmManager()) can claim all task slots.
stdbuf -oL -eL julia --startup-file=no ${SYSFLAG} --project="${FUSE_ROOT}" \
    "${TJLFEP_ROOT}/build/timing/run_iter_fuse_scan20_gpu.jl"

# --- tear down MPS ---
if [ "${INNER}" = "mps_team" ]; then
  srun --ntasks="${SLURM_NNODES}" --ntasks-per-node=1 bash -lc 'echo quit | nvidia-cuda-mps-control || true' || true
fi

echo "=== done; see TIMING_RESULT markers above ==="
"""

# SPMD layout (inner=:mps_team): the verified DIII-D gacode path, but driving the FUSE dd.
# Phase 1 builds the dd once (master); phase 2 is `srun -n SCAN_N` through mps-scan-wrapper.sh
# (1 radius : 1 GPU + an MPS team sharing that GPU); phase 3 merges via the normal ActorTJLFEP
# (TJLFEP_PRECOMPUTED_DIR short-circuits runTHD's pmap). No single-master/pmap topology, so the
# per-radius team `addprocs` always runs from a real task master.
const _TJLFEP_SPMD_BATCH_TEMPLATE = raw"""#!/bin/bash -l
# Auto-generated by TurbulentTransport.run_tjlfep (inner=mps_team SPMD) -- FUSE-dd TJLFEP scan
# on the @NODES@-node / @NTASKS@-GPU layout, 1 radius : 1 GPU + an MPS team per GPU.
#SBATCH -A @ACCOUNT@
#SBATCH -q @QOS@
#SBATCH -N @NODES@
#SBATCH -n @NTASKS@
#SBATCH -t @WALLTIME@
#SBATCH -C gpu
#SBATCH -J @JOBNAME@
#SBATCH -o @BASEDIR@/tjlfep_%j.out
#SBATCH -e @BASEDIR@/tjlfep_%j.err
#SBATCH --ntasks-per-node=@TPN@
#SBATCH --cpus-per-task=@CPT@
#SBATCH --gpus-per-node=@GPN@

set -uo pipefail

module load cudatoolkit/12.9
module load julia/1.11.7
export JULIA_DEPOT_PATH="${PSCRATCH}/.julia${JULIA_DEPOT_PATH:+:${JULIA_DEPOT_PATH}}"

export TJLFEP_ROOT="@TJLFEP_ROOT@"
export FUSE_ROOT="@FUSE_ROOT@"

export CASE="@CASE@"
export SCAN_N="@SCAN_N@"
export N_BASIS="@N_BASIS@"
export NGRID="@NGRID@"
export ALPHA_SOLVER="@ALPHA_SOLVER@"
export SOLVER="@SOLVER@"
export AD_EXTEND_MODE="@EXTEND_MODE@"
export AD_WIDE_KDESC="@WIDE_KDESC@"
export AD_FAITHFUL_CONFIRM="@FAITHFUL_CONFIRM@"
export INNER="@INNER@"
export MPS_TEAM="@MPS_TEAM@"
export TJLFEP_OUT_DIR="@BASEDIR@"
export JULIA_WORKER_THREADS="@WORKER_THREADS@"
export JULIA_CUDA_USE_COMPAT=false
export USE_GPU=1
export TJLFEP_GPU_SYSIMAGE="@SYSIMAGE@"

cd "${TJLFEP_ROOT}/build"
echo "=== run_tjlfep SPMD case=${CASE} job=${SLURM_JOB_ID} nodes=${SLURM_NNODES} tasks=${SLURM_NTASKS} ==="
echo "N_BASIS=${N_BASIS} SCAN_N=${SCAN_N} NGRID=${NGRID} ALPHA_SOLVER=${ALPHA_SOLVER} SOLVER=${SOLVER} EXTEND_MODE=${AD_EXTEND_MODE} WIDE_KDESC=${AD_WIDE_KDESC} FAITHFUL_CONFIRM=${AD_FAITHFUL_CONFIRM} INNER=${INNER} MPS_TEAM=${MPS_TEAM}"
nvidia-smi -L 2>/dev/null | head -4 || true

SYSFLAG=""
if [ -n "${TJLFEP_GPU_SYSIMAGE}" ] && [ -f "${TJLFEP_GPU_SYSIMAGE}" ]; then
  SYSFLAG="--sysimage=${TJLFEP_GPU_SYSIMAGE}"
  echo "sysimage = ${TJLFEP_GPU_SYSIMAGE}"
else
  echo "sysimage = <none, JIT from dev source>"
fi

# --- Phase 1: build the dd once on the head node + serialize inputs ---
echo "--- phase 1: prepare (build dd) ---"
stdbuf -oL -eL julia --startup-file=no ${SYSFLAG} --project="${FUSE_ROOT}" \
    "${TJLFEP_ROOT}/build/timing/run_fuse_dd_prepare.jl" || { echo "prepare FAILED"; exit 1; }

# --- Phase 2: one radius per GPU + MPS team (mps-scan-wrapper starts the per-node MPS daemon,
#     pins GPUs, and sets SCAN_INDEX = global procid + 1) ---
echo "--- phase 2: per-radius MPS scan (srun -n ${SCAN_N}) ---"
GPUS_PER_RADIUS=1 srun -n "${SCAN_N}" --ntasks-per-node=@TPN@ --gpus-per-node=@GPN@ \
    "${TJLFEP_ROOT}/build/common/mps-scan-wrapper.sh" \
    julia --startup-file=no ${SYSFLAG} --project="${FUSE_ROOT}" \
        "${TJLFEP_ROOT}/build/timing/run_fuse_dd_mps_task.jl" || { echo "per-radius scan FAILED"; exit 1; }

# --- Phase 3: merge per-radius results through the normal ActorTJLFEP (CPU, head node) ---
echo "--- phase 3: merge (ActorTJLFEP reads precomputed radii) ---"
stdbuf -oL -eL julia --startup-file=no ${SYSFLAG} --project="${FUSE_ROOT}" \
    "${TJLFEP_ROOT}/build/timing/run_fuse_dd_merge.jl" || { echo "merge FAILED"; exit 1; }

echo "=== done; see TIMING_RESULT markers above ==="
"""

"""
    run_tjlfep(case::Symbol=:ITER; kwargs...) -> TGLFEPRunState

Render and submit a SLURM job that runs the FUSE-dd TJLFEP energetic-particle
stability scan for `case` on the 5-node / 20-GPU layout, then return a
[`TGLFEPRunState`](@ref) to poll with [`refresh_tjlfep!`](@ref).

Keyword arguments (defaults reproduce the validated 5N/20-GPU premium/1h layout):
- `n_scan=20, n_basis=32, ngrid=201, gpu=true, alpha_solver=:stiff, nodes=5`
- `solver=:ad` (critical-factor engine, matches the `ActorTJLFEP` production default;
  `:grid` for Fortran-equivalence, `:robust_ad` for the most accurate AD path, `:truth`
  for the narrow-width EP-driven onset). `extend_mode=:locate` (`:wide` for the fast
  single-pass width-aware mode), `wide_kdesc=2` (`:wide` multistart breadth),
  `faithful_confirm=true`. Each defaults from its env var (`SOLVER`, `AD_EXTEND_MODE`,
  `AD_WIDE_KDESC`, `AD_FAITHFUL_CONFIRM`) so a preset shell env is honored; the rendered
  batch script then re-exports them so the scan/SPMD-task scripts pick them up.
- `inner=:threads` (proven 1-radius-per-GPU baseline, single-master + `pmap` topology) or
  `:mps_team` with `mps_team=8`. `:mps_team` renders the **SPMD layout** (the verified DIII-D
  gacode path): a 3-phase job — phase 1 builds the dd once, phase 2 is `srun -n n_scan` through
  `mps-scan-wrapper.sh` (1 radius : 1 GPU, each with an `mps_team`-sized pool of MPS clients
  sharing that GPU for the within-radius kw-scan), phase 3 merges via the normal `ActorTJLFEP`
  (`TJLFEP_PRECOMPUTED_DIR` makes `runTHD` load the per-radius results). The per-radius team
  `addprocs` runs from a genuine task master, so it is not the nested-from-`pmap`-worker case.
- `walltime="01:00:00", qos="premium", account="m3739_g", worker_threads=8`
- `basedir=""`: defaults to `default_results_dir("TJLFEP")`
  (`/global/cfs/cdirs/m3739/results/FUSE/TJLFEP/tjlfep_<user>_<timestamp>`).
- `sysimage`: shared CFS GPU sysimage (default `$TJLFEP_CFS_SYSIMAGE`). When present
  and SHA-matched to the current TJLFEP+TJLF source it is used by master AND workers;
  otherwise `rebuild_sysimage` controls a (re)build.
- `rebuild_sysimage=:if_stale` (`:always`/`:never`): if the sysimage is missing/stale,
  submit `build_script` and chain the scan with `--dependency=afterok`. `:never`
  (or no build script) with a stale image => JIT fallback (master+workers both JIT).
- `submit=true`: set `false` to only render the batch script (dry run).
- `wait_for_completion=false, poll_interval=60`.
"""
function run_tjlfep(case::Symbol=:ITER;
    n_scan::Int=20,
    n_basis::Int=32,
    ngrid::Int=201,
    gpu::Bool=true,
    alpha_solver::Symbol=:stiff,
    solver::Symbol=Symbol(get(ENV, "SOLVER", "ad")),
    extend_mode::Symbol=Symbol(get(ENV, "AD_EXTEND_MODE", "locate")),
    wide_kdesc::Int=parse(Int, get(ENV, "AD_WIDE_KDESC", "2")),
    faithful_confirm::Bool=get(ENV, "AD_FAITHFUL_CONFIRM", "1") != "0",
    nodes::Int=5,
    inner::Symbol=:threads,
    mps_team::Int=8,
    walltime::String="01:00:00",
    qos::String="premium",
    account::String="m3739_g",
    worker_threads::Int=8,
    basedir::String="",
    fuse_root::String=get(ENV, "FUSE_ROOT", normpath(pkgdir(@__MODULE__), "..", "FUSE")),
    tjlfep_root::String=get(ENV, "TJLFEP_ROOT", normpath(pkgdir(@__MODULE__), "..", "TJLFEP")),
    tjlf_root::String=pkgdir(TJLF),
    sysimage::String=TJLFEP_CFS_SYSIMAGE,
    build_script::String="",
    rebuild_sysimage::Symbol=:if_stale,
    jobname::String="tjlfep_$(lowercase(string(case)))",
    submit::Bool=true,
    wait_for_completion::Bool=false,
    poll_interval::Int=60)

    isempty(basedir) && (basedir = default_results_dir("TJLFEP"))
    mkpath(basedir)

    ntasks = n_scan
    ntasks_per_node = cld(n_scan, nodes)
    gpus_per_node = min(ntasks_per_node, 4)
    cpus_per_task = 32

    # Sysimage staleness policy (Section 7): default to the shared CFS image, rebuild on change.
    use_sysimage = ""
    build_job_id = ""
    if gpu && !isempty(sysimage)
        ok = _tjlfep_sysimage_ok(sysimage, tjlfep_root, tjlf_root)
        if ok && rebuild_sysimage != :always
            use_sysimage = sysimage
        elseif rebuild_sysimage in (:if_stale, :always)
            bs = isempty(build_script) ?
                 joinpath(tjlfep_root, "build", "sysimage", "batch_build_gpu_sysimage_generic.sh") :
                 build_script
            if isfile(bs)
                out = read(`sbatch $bs`, String)
                build_job_id = parse_slurm_jobid(out)
                use_sysimage = sysimage  # available once the afterok build completes
                @info "run_tjlfep: submitted sysimage (re)build" build_job_id script = bs target = sysimage
            else
                @warn "run_tjlfep: sysimage stale/missing and build script not found; JIT fallback" script = bs
            end
        else
            @warn "run_tjlfep: sysimage stale/missing and rebuild_sysimage=:never; JIT fallback (master+workers JIT)"
        end
    end

    template = inner === :mps_team ? _TJLFEP_SPMD_BATCH_TEMPLATE : _TJLFEP_BATCH_TEMPLATE
    script = replace(template,
        "@ACCOUNT@" => account,
        "@QOS@" => qos,
        "@NODES@" => string(nodes),
        "@NTASKS@" => string(ntasks),
        "@WALLTIME@" => walltime,
        "@JOBNAME@" => jobname,
        "@BASEDIR@" => basedir,
        "@TPN@" => string(ntasks_per_node),
        "@CPT@" => string(cpus_per_task),
        "@GPN@" => string(gpus_per_node),
        "@TJLFEP_ROOT@" => tjlfep_root,
        "@FUSE_ROOT@" => fuse_root,
        "@CASE@" => string(case),
        "@SCAN_N@" => string(n_scan),
        "@N_BASIS@" => string(n_basis),
        "@NGRID@" => string(ngrid),
        "@ALPHA_SOLVER@" => string(alpha_solver),
        "@SOLVER@" => string(solver),
        "@EXTEND_MODE@" => string(extend_mode),
        "@WIDE_KDESC@" => string(wide_kdesc),
        "@FAITHFUL_CONFIRM@" => faithful_confirm ? "1" : "0",
        "@INNER@" => string(inner),
        "@MPS_TEAM@" => string(mps_team),
        "@WORKER_THREADS@" => string(worker_threads),
        "@SYSIMAGE@" => use_sysimage)

    batchfile = joinpath(basedir, "submit_tjlfep.sh")
    open(io -> write(io, script), batchfile, "w")
    chmod(batchfile, 0o755)

    job_id = ""
    if submit
        cmd = isempty(build_job_id) ? `sbatch $batchfile` :
              `sbatch --dependency=afterok:$build_job_id $batchfile`
        out = read(cmd, String)
        job_id = parse_slurm_jobid(out)
        @info "run_tjlfep: submitted scan" job_id case nodes n_scan n_basis solver extend_mode inner basedir
    else
        @info "run_tjlfep: rendered batch script (submit=false)" batchfile
    end

    logfile = isempty(job_id) ? "" : joinpath(basedir, "tjlfep_$(job_id).out")
    state = TGLFEPRunState(basedir, case, n_scan, n_basis, ngrid, gpu, alpha_solver,
        nodes, inner, mps_team, job_id, build_job_id, batchfile, logfile, use_sysimage,
        submit ? :submitted : :rendered)
    save_state(state, joinpath(basedir, ".tjlfep_state.jls"))

    if submit && wait_for_completion
        while true
            refresh_tjlfep!(state)
            state.status in (:completed, :failed, :unknown) && break
            sleep(poll_interval)
        end
    end

    return state
end

"""
    refresh_tjlfep!(state::TGLFEPRunState) -> TGLFEPRunState

Poll `sacct` for the scan job and update `state.status`. Once the job leaves the
queue, resolves `state.logfile` (the `%j` output) if not already known. Saves state.
"""
function refresh_tjlfep!(state::TGLFEPRunState)
    if isempty(state.job_id)
        state.status = :rendered
        return state
    end
    state.status = check_slurm_status(state.job_id)
    if isempty(state.logfile) || !isfile(state.logfile)
        cand = joinpath(state.basedir, "tjlfep_$(state.job_id).out")
        isfile(cand) && (state.logfile = cand)
    end
    save_state(state, joinpath(state.basedir, ".tjlfep_state.jls"))
    return state
end

"""
    tjlfep_status(state::TGLFEPRunState) -> String

Human-readable one-line status. When the job has finished, appends the last
`TIMING_RESULT` total-job line from the log if available.
"""
function tjlfep_status(state::TGLFEPRunState)
    base = "TJLFEP $(state.case) job=$(state.job_id) [$(state.status)] " *
           "nodes=$(state.nodes) n_scan=$(state.n_scan) n_basis=$(state.n_basis) inner=$(state.inner)"
    if state.status in (:completed, :failed) && !isempty(state.logfile) && isfile(state.logfile)
        timing = ""
        for ln in eachline(state.logfile)
            occursin("TIMING_RESULT", ln) && occursin("phase=total_job", ln) && (timing = strip(ln))
        end
        !isempty(timing) && (base *= "\n  " * timing)
    end
    return base
end

"""
    load_tjlfep_results(state::TGLFEPRunState) -> NamedTuple

Load the per-run results the driver serialized to `basedir/tjlfep_results.jls`
(`(; rho_scan, SFmin, width, kymark, n_EP, p_EP)`), plus the path to `dd_out.json`
if present. Returns `nothing` if results are not yet written.
"""
function load_tjlfep_results(state::TGLFEPRunState)
    results = joinpath(state.basedir, "tjlfep_results.jls")
    isfile(results) || return nothing
    res = load_state(results)
    dd_json = joinpath(state.basedir, "dd_out.json")
    return (; res..., dd_json = isfile(dd_json) ? dd_json : nothing)
end
