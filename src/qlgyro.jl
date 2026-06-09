using Serialization
import Dates

Base.@kwdef mutable struct InputQLGYRO
    N_PARALLEL::Union{Int,Missing} = missing
    N_RUNS::Union{Int,Missing} = missing
    GAMMA_E::Union{Float64,Missing} = missing
    CODE::Union{Int,Missing} = missing
    NKY::Union{Int,Missing} = missing
    KYGRID_MODEL::Union{Int,Missing} = missing
    KY::Union{Float64,Missing} = missing
    XNU_MODEL::Union{Int,Missing} = missing
    AUTO_BOX_SIZE::Union{Int,Missing} = missing
    KX_MAX_BOX::Union{Float64,Missing} = missing
    SAT_RULE::Union{Int,Missing} = missing
    N_PX0::Union{Int,Missing} = missing
    PX0GRID_MODEL::Union{Int,Missing} = missing
    RESTART_MODE::Union{Int,Missing} = missing
end

"""
    run_qlgyro(input_qlgyro::InputQLGYRO, input_cgyro::InputCGYRO)

Run QLGYRO starting from a InputQLGYRO and InputCGYRO

Returns a `FluxSolution` structure
"""
function run_qlgyro(input_qlgyro::InputQLGYRO, input_cgyro::InputCGYRO)
    folder = mktempdir()

    save(input_cgyro, joinpath(folder, "input.cgyro"))
    save(input_qlgyro, joinpath(folder, "input.qlgyro"))

    n_parallel = input_qlgyro.N_PARALLEL
    preamble = gacode_preamble()
    open(joinpath(folder, "command.sh"), "w") do io
        return write(
            io,
            """
         $preamble
         	(time (qlgyro -n $n_parallel -e .)) &> command.log
         """)
    end

    fluxes = try
        run(Cmd(`bash command.sh`; dir=folder))

        tmp = open(joinpath(folder, "out.qlgyro.gbflux"), "r") do io
            return read(io, String)
        end
        fluxes = parse_out_tglf_gbflux(tmp)

    catch e
        # show last 100 lines of  chease.output
        txt = open(joinpath(folder, "command.log"), "r") do io
            return split(read(io, String), "\n")
        end
        @error "ERROR running QLGYRO\n...\n" * join(txt[max(1, length(txt) - 100):end], "\n")
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
    run_qlgyro(input_qlgyros::Vector{InputQLGYRO}, input_cgyros::Vector{InputCGYRO})

Run QLGYRO starting from a vectors of InputQLGYRO and InputCGYRO

NOTE: Each run is done sequentially, one after the other

Returns a vector of `FluxSolution` structures
"""
function run_qlgyro(input_qlgyros::Vector{InputQLGYRO}, input_cgyros::Vector{InputCGYRO})
    @assert length(input_qlgyros) == length(input_cgyros)
    return collect(map((input_qlgyro, input_cgyro) -> run_qlgyro(input_qlgyro, input_cgyro), input_qlgyros, input_cgyros))
end

# ============================================================================
# Modular QLGYRO: linear CGYRO + TJLF saturation rules
# ============================================================================

function _require_env(key::String)
    val = get(ENV, key, "")
    isempty(val) && error("Environment variable $key is not set. " *
        "Set it to your GACODE installation path, e.g. export $key=/path/to/gacode")
    return val
end

function gacode_env(gpu::Bool)
    root = gpu ? _require_env("GACODE_ROOT_GPU") : _require_env("GACODE_ROOT_CPU")
    platform = gpu ? "PERLMUTTER_GPU" : "PERLMUTTER_CPU"
    return root, platform
end

"""
    tglf_to_cgyro(input_tglf::InputTGLF) -> InputCGYRO

Convert InputTGLF parameters to InputCGYRO.

Species ordering: TGLF has electrons as species 1, ions as 2..NS.
CGYRO has ions first (1..NS-1), electrons last (NS).
"""
function tglf_to_cgyro(input_tglf::InputTGLF)
    ic = InputCGYRO()

    ns = input_tglf.NS
    ic.N_SPECIES = ns

    # Geometry (Miller local equilibrium)
    ic.EQUILIBRIUM_MODEL = 2  # Miller Extended Harmonic
    ic.RMIN = input_tglf.RMIN_LOC
    ic.RMAJ = input_tglf.RMAJ_LOC
    ic.Q = input_tglf.Q_LOC
    ic.KAPPA = input_tglf.KAPPA_LOC
    ic.DELTA = input_tglf.DELTA_LOC
    ic.ZETA = input_tglf.ZETA_LOC
    ic.S_KAPPA = input_tglf.S_KAPPA_LOC
    ic.S_DELTA = input_tglf.S_DELTA_LOC
    ic.S_ZETA = input_tglf.S_ZETA_LOC
    ic.ZMAG = input_tglf.ZMAJ_LOC
    ic.DZMAG = input_tglf.DZMAJDX_LOC
    ic.SHIFT = input_tglf.DRMAJDX_LOC

    # Higher-order MXH shaping coefficients (copy fields that exist in both structs)
    cgyro_fields = fieldnames(typeof(ic))
    for m in 0:3
        for prefix in ("SHAPE_COS", "SHAPE_S_COS")
            sym = Symbol("$(prefix)$m")
            sym in cgyro_fields || continue
            val = getproperty(input_tglf, sym)
            ismissing(val) || setproperty!(ic, sym, val)
        end
    end
    for prefix in ("SHAPE_SIN", "SHAPE_S_SIN")
        sym = Symbol("$(prefix)3")
        sym in cgyro_fields || continue
        val = getproperty(input_tglf, sym)
        ismissing(val) || setproperty!(ic, sym, val)
    end

    # Magnetic shear conversion:
    # TGLF: Q_PRIME_LOC = q * (a²/r) * dq/dr (from tglf.jl IMAS constructor)
    # CGYRO: S = (r/q) * dq/dr
    # Therefore: S = RMIN_LOC² * Q_PRIME_LOC / Q_LOC²
    ic.S = input_tglf.RMIN_LOC^2 * input_tglf.Q_PRIME_LOC / input_tglf.Q_LOC^2

    # Signs
    ic.BTCCW = Float64(input_tglf.SIGN_BT)
    ic.IPCCW = Float64(input_tglf.SIGN_IT)

    # Electromagnetic fields
    n_field = 1
    if input_tglf.USE_BPER
        n_field = 2
    end
    if input_tglf.USE_BPAR
        n_field = 3
    end
    ic.N_FIELD = n_field
    ic.BETAE_UNIT = input_tglf.BETAE

    # Debye length (CGYRO calls it LAMBDA_STAR)
    ic.LAMBDA_STAR = input_tglf.DEBYE

    # Collision model
    ic.NU_EE = input_tglf.XNUE
    ic.COLLISION_MODEL = 4  # Sugama (standard for CGYRO linear)

    # Rotation
    # TGLF VEXB_SHEAR -> CGYRO GAMMA_E
    # These are normalized differently. For LOCAL model:
    # TGLF: VEXB_SHEAR = -(r/q)*(d omega_0/d r) * (a/c_s) (in units_in=GYRO)
    # CGYRO: GAMMA_E = -(r/q)*(d omega_0/d r) * (a/c_s)
    # So GAMMA_E = VEXB_SHEAR for GYRO units
    if !ismissing(input_tglf.VEXB_SHEAR) && input_tglf.VEXB_SHEAR != 0.0
        ic.GAMMA_E = input_tglf.VEXB_SHEAR
    end

    # Species mapping: TGLF species 1 = electrons, 2..NS = ions
    # CGYRO species 1..NS-1 = ions, NS = electrons
    for i_tglf in 2:ns
        i_cgyro = i_tglf - 1  # ion index in CGYRO
        setproperty!(ic, Symbol("Z_$i_cgyro"), getfield(input_tglf, Symbol("ZS_$i_tglf")))
        setproperty!(ic, Symbol("MASS_$i_cgyro"), getfield(input_tglf, Symbol("MASS_$i_tglf")))
        # TGLF: AS = n_s/n_e, CGYRO: DENS = n_s/n_e (both normalized to electron density)
        setproperty!(ic, Symbol("DENS_$i_cgyro"), getfield(input_tglf, Symbol("AS_$i_tglf")))
        # TGLF: TAUS = T_s/T_e, CGYRO: TEMP = T_s/T_e
        setproperty!(ic, Symbol("TEMP_$i_cgyro"), getfield(input_tglf, Symbol("TAUS_$i_tglf")))
        # Gradients: TGLF RLNS = -a/n * dn/dr, CGYRO DLNNDR = -a/n * dn/dr
        setproperty!(ic, Symbol("DLNNDR_$i_cgyro"), getfield(input_tglf, Symbol("RLNS_$i_tglf")))
        setproperty!(ic, Symbol("DLNTDR_$i_cgyro"), getfield(input_tglf, Symbol("RLTS_$i_tglf")))
    end

    # Electrons last in CGYRO
    i_e_cgyro = ns
    setproperty!(ic, Symbol("Z_$i_e_cgyro"), getfield(input_tglf, :ZS_1))
    setproperty!(ic, Symbol("MASS_$i_e_cgyro"), getfield(input_tglf, :MASS_1))
    setproperty!(ic, Symbol("DENS_$i_e_cgyro"), getfield(input_tglf, :AS_1))
    setproperty!(ic, Symbol("TEMP_$i_e_cgyro"), getfield(input_tglf, :TAUS_1))
    setproperty!(ic, Symbol("DLNNDR_$i_e_cgyro"), getfield(input_tglf, :RLNS_1))
    setproperty!(ic, Symbol("DLNTDR_$i_e_cgyro"), getfield(input_tglf, :RLTS_1))

    # Rotation model (2 = Sonic, matching nstx_20/d3d_40 conventions)
    ic.ROTATION_MODEL = 2

    # Linear run settings (defaults from nstx_20/generate_ky_scan.py)
    ic.NONLINEAR_FLAG = 0
    ic.N_TOROIDAL = 1
    ic.N_RADIAL = 16
    ic.BOX_SIZE = 1
    ic.N_THETA = 24
    ic.N_XI = 24
    ic.N_ENERGY = 8
    ic.DELTA_T = 0.005
    ic.DELTA_T_METHOD = 1   # Cash-Karp adaptive timestepping
    ic.ERROR_TOL = 0.001
    ic.MAX_TIME = 100000.0
    ic.PRINT_STEP = 100.0
    ic.FREQ_TOL = 0.01
    ic.FIELD_PRINT_FLAG = 1

    return ic
end

"""
    build_input_tjlf(input_cgyro::InputCGYRO, input_qlgyro::InputQLGYRO;
                     ky_values::Union{Nothing,Vector{Float64}}=nothing,
                     nmodes::Int=1) -> InputTJLF

Construct an `InputTJLF` directly from a `(InputCGYRO, InputQLGYRO)` pair. Performs:
- Miller geometry name translation (RMIN→RMIN_LOC, etc.) and `S → Q_PRIME_LOC` conversion
- Species reordering (CGYRO ions-first → TJLF electrons-first)
- EM-field translation (`N_FIELD` → `USE_BPER`/`USE_BPAR`), unit copies (BETAE_UNIT, NU_EE, LAMBDA_STAR)
- Pulls `KYGRID_MODEL`/`NKY`/`KY`/`SAT_RULE`/`XNU_MODEL` from `input_qlgyro`
- Sets sensible TJLF defaults for parameters not present in CGYRO/QLGYRO (WIDTH, ALPHA_ZF, etc.)

If `ky_values` is provided, `KY_SPECTRUM` is pre-populated and `KYGRID_MODEL` is forced to 0,
matching the modular QLGYRO flow where ky points come from the CGYRO scan.
"""
function build_input_tjlf(input_cgyro::InputCGYRO, input_qlgyro::InputQLGYRO;
                          ky_values::Union{Nothing,Vector{Float64}}=nothing,
                          nmodes::Int=1)
    ns = input_cgyro.N_SPECIES
    ns === missing && error("InputCGYRO.N_SPECIES must be set")

    nky = ky_values === nothing ? Int(input_qlgyro.NKY) : length(ky_values)
    input_tjlf = InputTJLF{Float64}(ns, nky)

    # ---- Miller geometry (CGYRO names → TJLF *_LOC names) ----
    input_tjlf.RMIN_LOC = input_cgyro.RMIN
    input_tjlf.RMAJ_LOC = input_cgyro.RMAJ
    input_tjlf.Q_LOC = input_cgyro.Q
    input_tjlf.KAPPA_LOC = input_cgyro.KAPPA
    input_tjlf.DELTA_LOC = input_cgyro.DELTA
    input_tjlf.ZETA_LOC = input_cgyro.ZETA
    input_tjlf.S_KAPPA_LOC = input_cgyro.S_KAPPA
    input_tjlf.S_DELTA_LOC = input_cgyro.S_DELTA
    input_tjlf.S_ZETA_LOC = input_cgyro.S_ZETA
    input_tjlf.ZMAJ_LOC = input_cgyro.ZMAG
    input_tjlf.DZMAJDX_LOC = input_cgyro.DZMAG
    input_tjlf.DRMAJDX_LOC = input_cgyro.SHIFT
    input_tjlf.DRMINDX_LOC = 1.0

    # Magnetic shear: CGYRO S = (r/q) dq/dr; TGLF Q_PRIME_LOC = q (a²/r) dq/dr
    # Therefore Q_PRIME_LOC = (Q² / RMIN²) * S
    if !ismissing(input_cgyro.S) && !ismissing(input_cgyro.Q) && !ismissing(input_cgyro.RMIN)
        input_tjlf.Q_PRIME_LOC = input_cgyro.Q^2 / input_cgyro.RMIN^2 * input_cgyro.S
    end
    input_tjlf.P_PRIME_LOC = 0.0
    input_tjlf.BETA_LOC = 0.0
    input_tjlf.KX0_LOC = 0.0

    # MXH shape coefficients (only fields present in InputCGYRO)
    tjlf_fields = fieldnames(typeof(input_tjlf))
    for m in 0:3
        for prefix in ("SHAPE_COS", "SHAPE_S_COS")
            sym = Symbol("$(prefix)$m")
            sym in tjlf_fields || continue
            sym in fieldnames(InputCGYRO) || continue
            val = getproperty(input_cgyro, sym)
            ismissing(val) || setproperty!(input_tjlf, sym, Float64(val))
        end
    end
    for prefix in ("SHAPE_SIN", "SHAPE_S_SIN")
        sym = Symbol("$(prefix)3")
        sym in tjlf_fields || continue
        sym in fieldnames(InputCGYRO) || continue
        val = getproperty(input_cgyro, sym)
        ismissing(val) || setproperty!(input_tjlf, sym, Float64(val))
    end

    # ---- Signs ----
    input_tjlf.SIGN_BT = ismissing(input_cgyro.BTCCW) ? 1 : Int(sign(input_cgyro.BTCCW))
    input_tjlf.SIGN_IT = ismissing(input_cgyro.IPCCW) ? 1 : Int(sign(input_cgyro.IPCCW))

    # ---- EM fields ----
    nfield = ismissing(input_cgyro.N_FIELD) ? 1 : input_cgyro.N_FIELD
    input_tjlf.USE_BPER = nfield >= 2
    input_tjlf.USE_BPAR = nfield >= 3
    input_tjlf.BETAE = ismissing(input_cgyro.BETAE_UNIT) ? 0.0 : input_cgyro.BETAE_UNIT

    # ---- Collisions, Debye ----
    input_tjlf.XNUE = ismissing(input_cgyro.NU_EE) ? 0.0 : input_cgyro.NU_EE
    input_tjlf.DEBYE = ismissing(input_cgyro.LAMBDA_STAR) ? 0.0 : input_cgyro.LAMBDA_STAR

    # ---- Rotation: CGYRO GAMMA_E = TGLF VEXB_SHEAR (GYRO units) ----
    gamma_e = !ismissing(input_qlgyro.GAMMA_E) ? input_qlgyro.GAMMA_E :
              (ismissing(input_cgyro.GAMMA_E) ? 0.0 : input_cgyro.GAMMA_E)
    input_tjlf.VEXB_SHEAR = Float64(gamma_e)

    # ---- Species: CGYRO ions(1..NS-1) + electrons(NS) → TJLF electrons(1) + ions(2..NS) ----
    # Electrons (CGYRO species NS → TJLF species 1)
    e_idx = ns
    input_tjlf.ZS[1] = Float64(getproperty(input_cgyro, Symbol("Z_$e_idx")))
    input_tjlf.MASS[1] = Float64(getproperty(input_cgyro, Symbol("MASS_$e_idx")))
    input_tjlf.AS[1] = Float64(getproperty(input_cgyro, Symbol("DENS_$e_idx")))
    input_tjlf.TAUS[1] = Float64(getproperty(input_cgyro, Symbol("TEMP_$e_idx")))
    input_tjlf.RLNS[1] = Float64(getproperty(input_cgyro, Symbol("DLNNDR_$e_idx")))
    input_tjlf.RLTS[1] = Float64(getproperty(input_cgyro, Symbol("DLNTDR_$e_idx")))
    input_tjlf.VPAR[1] = 0.0
    input_tjlf.VPAR_SHEAR[1] = 0.0

    # Ions (CGYRO species 1..NS-1 → TJLF species 2..NS)
    for ic in 1:(ns - 1)
        it = ic + 1
        input_tjlf.ZS[it] = Float64(getproperty(input_cgyro, Symbol("Z_$ic")))
        input_tjlf.MASS[it] = Float64(getproperty(input_cgyro, Symbol("MASS_$ic")))
        input_tjlf.AS[it] = Float64(getproperty(input_cgyro, Symbol("DENS_$ic")))
        input_tjlf.TAUS[it] = Float64(getproperty(input_cgyro, Symbol("TEMP_$ic")))
        input_tjlf.RLNS[it] = Float64(getproperty(input_cgyro, Symbol("DLNNDR_$ic")))
        input_tjlf.RLTS[it] = Float64(getproperty(input_cgyro, Symbol("DLNTDR_$ic")))
        input_tjlf.VPAR[it] = 0.0
        input_tjlf.VPAR_SHEAR[it] = 0.0
    end

    # ZEFF from quasineutrality: ZEFF = Σ_i n_i Z_i² / n_e (with n_e normalized to 1)
    zeff = 0.0
    n_e = input_tjlf.AS[1]
    if n_e > 0
        for is in 2:ns
            zeff += input_tjlf.AS[is] * input_tjlf.ZS[is]^2
        end
        zeff /= n_e
    end
    input_tjlf.ZEFF = zeff > 0 ? zeff : 1.0

    # ---- ky grid + saturation rule from InputQLGYRO ----
    input_tjlf.SAT_RULE = ismissing(input_qlgyro.SAT_RULE) ? 1 : input_qlgyro.SAT_RULE
    input_tjlf.NMODES = nmodes
    input_tjlf.NKY = ismissing(input_qlgyro.NKY) ? 12 : input_qlgyro.NKY
    if !ismissing(input_qlgyro.KY)
        input_tjlf.KY = Float64(input_qlgyro.KY)
    end
    input_tjlf.KYGRID_MODEL = ismissing(input_qlgyro.KYGRID_MODEL) ? 1 : input_qlgyro.KYGRID_MODEL

    # XNU_MODEL: respect input.qlgyro if set, else apply SAT_RULE preset
    if !ismissing(input_qlgyro.XNU_MODEL)
        input_tjlf.XNU_MODEL = input_qlgyro.XNU_MODEL
    else
        input_tjlf.XNU_MODEL = input_tjlf.SAT_RULE in (2, 3) ? 3 : 2
    end

    # If caller pre-computed the ky scan points, populate KY_SPECTRUM and force KYGRID_MODEL=0
    if ky_values !== nothing
        input_tjlf.KYGRID_MODEL = 0
        input_tjlf.NKY = nky
        input_tjlf.KY_SPECTRUM .= ky_values
    end

    # ---- TJLF defaults (mirroring InputTGLF defaults at TJLF/src/tjlf_modules.jl:142-176) ----
    input_tjlf.UNITS = input_tjlf.SAT_RULE in (2, 3) ? "CGYRO" : "GYRO"
    input_tjlf.USE_BISECTION = true
    input_tjlf.USE_INBOARD_DETRAPPED = false
    input_tjlf.USE_AVE_ION_GRID = false
    input_tjlf.USE_MHD_RULE = false
    input_tjlf.NEW_EIKONAL = true
    input_tjlf.FIND_WIDTH = true
    input_tjlf.FIND_EIGEN = true
    input_tjlf.IFLUX = true
    input_tjlf.ADIABATIC_ELEC = false
    input_tjlf.NWIDTH = 21
    input_tjlf.NXGRID = 16
    input_tjlf.NBASIS_MIN = 2
    input_tjlf.NBASIS_MAX = 4
    input_tjlf.VPAR_MODEL = 0
    input_tjlf.IBRANCH = -1

    input_tjlf.WIDTH = 1.65
    input_tjlf.WIDTH_MIN = 0.3
    input_tjlf.ETG_FACTOR = 1.25
    input_tjlf.RLNP_CUTOFF = 18.0
    input_tjlf.ALPHA_E = 1.0
    input_tjlf.ALPHA_P = 1.0
    input_tjlf.ALPHA_MACH = 0.0
    input_tjlf.ALPHA_QUENCH = 0
    input_tjlf.ALPHA_ZF = 1.0
    input_tjlf.XNU_FACTOR = 1.0
    input_tjlf.DEBYE_FACTOR = 1.0

    input_tjlf.PARK = 1.0
    input_tjlf.GHAT = 1.0
    input_tjlf.GCHAT = 1.0
    input_tjlf.WD_ZERO = 0.1
    input_tjlf.LINSKER_FACTOR = 0.0
    input_tjlf.GRADB_FACTOR = 0.0
    input_tjlf.FILTER = 2.0
    input_tjlf.THETA_TRAPPED = 0.7
    input_tjlf.DAMP_PSI = 0.0
    input_tjlf.DAMP_SIG = 0.0
    input_tjlf.WDIA_TRAPPED = input_tjlf.SAT_RULE in (2, 3) ? 1.0 : 0.0

    input_tjlf.WIDTH_SPECTRUM .= input_tjlf.WIDTH

    return input_tjlf
end

"""
    generate_ky_grid(input_tglf::InputTGLF; nky::Int=0, ky::Float64=0.0, kygrid_model::Int=-1) -> Vector{Float64}

Generate a TJLF ky grid using the specified KYGRID_MODEL.

- `kygrid_model=-1` (default): uses `input_tglf.KYGRID_MODEL`
- `kygrid_model=0`: simple linear grid, `ky[i] = i * KY / NKY` (defaults: KY=1.2, NKY=12)
- `kygrid_model=1`: APS07 — 9 low-ky linear + NKY log-spaced ETG points
- `kygrid_model=2`: IAEA08 — 15 low-ky linear + NKY log-spaced ETG points
- `kygrid_model=3`: similar to APS07 with `ky_min = KY`
- `kygrid_model=4`: 12 low-ky linear (`ky_min=0.05/ρ_i`) + NKY log-spaced ETG points (recommended for SAT2/SAT3)

Models 1-4 require species info (MASS, TAUS, ZS) and geometry (for `grad_r0`).
Set `nky` to override `input_tglf.NKY` (controls the number of ETG-scale points for models 1-4).
Set `ky` to override the max ky for model 0.
"""
function generate_ky_grid(input_tglf::InputTGLF; nky::Int=0, ky::Float64=0.0, kygrid_model::Int=-1)
    model = kygrid_model >= 0 ? kygrid_model : input_tglf.KYGRID_MODEL
    nky_in = nky > 0 ? nky : input_tglf.NKY

    if model == 0
        # Simple linear grid (same as TJLF KYGRID_MODEL=0)
        # When explicitly requesting model 0, default to KY=1.2
        ky_max = ky > 0 ? ky : (kygrid_model == 0 ? 1.2 : input_tglf.KY)
        dky = ky_max / nky_in
        return [i * dky for i in 1:nky_in]
    else
        # Delegate to TJLF's get_ky_spectrum for models 1-4
        nky_total = TJLF.get_ky_spectrum_size(nky_in, model)
        input_tjlf = InputTJLF{Float64}(input_tglf.NS, nky_total)
        TJLF.update_input_tjlf!(input_tjlf, input_tglf)
        input_tjlf.KYGRID_MODEL = model
        input_tjlf.NKY = nky_in

        # Compute geometry for grad_r0 (needed for ky scaling in non-GYRO units)
        satParams = TJLF.get_sat_params(input_tjlf)
        return TJLF.get_ky_spectrum(input_tjlf, satParams.grad_r0)
    end
end

"""
    generate_ky_grid(input_cgyro::InputCGYRO, input_qlgyro::InputQLGYRO) -> Vector{Float64}

Generate the ky scan grid using `KYGRID_MODEL`/`NKY`/`KY` from `input_qlgyro` and
geometry/species data from `input_cgyro`. Mirrors the InputTGLF version but builds
the temporary `InputTJLF` directly via `build_input_tjlf`.
"""
function generate_ky_grid(input_cgyro::InputCGYRO, input_qlgyro::InputQLGYRO)
    model = ismissing(input_qlgyro.KYGRID_MODEL) ? 1 : input_qlgyro.KYGRID_MODEL
    nky_in = ismissing(input_qlgyro.NKY) ? 12 : input_qlgyro.NKY

    if model == 0
        ky_max = ismissing(input_qlgyro.KY) ? 1.2 : input_qlgyro.KY
        dky = ky_max / nky_in
        return [i * dky for i in 1:nky_in]
    else
        nky_total = TJLF.get_ky_spectrum_size(nky_in, model)
        # Build a sized-up InputTJLF so KY_SPECTRUM has nky_total slots
        input_tjlf = build_input_tjlf(input_cgyro, input_qlgyro)
        # Resize KY_SPECTRUM/WIDTH_SPECTRUM/EIGEN_SPECTRUM via fresh construction
        input_tjlf2 = InputTJLF{Float64}(input_cgyro.N_SPECIES, nky_total)
        for f in fieldnames(typeof(input_tjlf2))
            f in (:KY_SPECTRUM, :WIDTH_SPECTRUM, :EIGEN_SPECTRUM) && continue
            setfield!(input_tjlf2, f, getfield(input_tjlf, f))
        end
        input_tjlf2.KYGRID_MODEL = model
        input_tjlf2.NKY = nky_in
        input_tjlf2.WIDTH_SPECTRUM .= input_tjlf2.WIDTH

        satParams = TJLF.get_sat_params(input_tjlf2)
        return TJLF.get_ky_spectrum(input_tjlf2, satParams.grad_r0)
    end
end

"""
    qlgyro_run_hash(input_cgyro::InputCGYRO) -> UInt64

Compute a hash of the InputCGYRO parameters for caching.
"""
function qlgyro_run_hash(input_cgyro::InputCGYRO)
    buf = IOBuffer()
    for key in fieldnames(InputCGYRO)
        val = getfield(input_cgyro, key)
        if !ismissing(val)
            print(buf, key, "=", val, "\n")
        end
    end
    return hash(String(take!(buf)))
end

"""
    QLGYRORunState

Tracks the state of a modular QLGYRO run (linear CGYRO scans).
"""
mutable struct QLGYRORunState
    basedir::String
    ky_values::Vector{Float64}
    slurm_ids::Vector{String}  # Slurm job IDs per KY
    converged::Vector{Bool}
    submitted::Vector{Bool}
    input_hash::UInt64
end

"""
    save_run_state(state::QLGYRORunState)

Serialize run state to disk for persistence across calls.
"""
function save_run_state(state::QLGYRORunState)
    Serialization.serialize(joinpath(state.basedir, ".qlgyro_state.jls"), state)
end

"""
    load_run_state(basedir::String) -> Union{QLGYRORunState, Nothing}

Load previously saved run state, or nothing if not found.
"""
function load_run_state(basedir::String)
    statefile = joinpath(basedir, ".qlgyro_state.jls")
    if isfile(statefile)
        return Serialization.deserialize(statefile)
    end
    return nothing
end

"""
    check_cgyro_convergence(rundir::String) -> Symbol

Check convergence status of a CGYRO run directory.
Returns :converged, :running, :error, :terminated, or :not_started.
"""
function check_cgyro_convergence(rundir::String)
    infofile = joinpath(rundir, "out.cgyro.info")
    timefile = joinpath(rundir, "out.cgyro.time")

    if !isfile(infofile)
        return :not_started
    end

    info = read(infofile, String)

    if occursin("Linear converged", info) || occursin("Underflow in calculation of frequency error", info)
        return :converged
    end

    # Check for timeout (high-ky runs may reach MAX_TIME without converging)
    if isfile(timefile)
        lines = readlines(timefile)
        if !isempty(lines)
            last_time = parse(Float64, split(strip(lines[end]))[1])
            if last_time >= 100.0
                return :converged  # treat timeout as "converged enough"
            end
        end
    end

    if occursin("Restart", info) && !occursin("EXIT", info)
        return :running
    end

    if occursin("ERROR", info)
        return :error
    end

    if occursin("terminated at max time", info)
        return :terminated
    end

    return :running
end

# check_slurm_status is provided by slurm_utils.jl (included before this file).

"""
    submit_cgyro_job(rundir::String; gpu::Bool=true, n_mpi::Int=32, n_omp::Int=4, walltime::String="00:15:00", repo::String="m3739_g") -> String

Submit a CGYRO job via gacode_qsub. Returns the Slurm job ID.
Set `gpu=true` (default) for GPU queue, `gpu=false` for CPU queue.
`qos` sets the Slurm QOS (e.g. "regular", "premium", "debug").
"""
function submit_cgyro_job(rundir::String; gpu::Bool=true, n_mpi::Int=32, n_omp::Int=4, walltime::String="00:15:00", repo::String="m3739_g", qos::String="regular")
    # Enforce debug QOS 30-minute walltime limit
    if qos == "debug"
        parts = split(walltime, ':')
        h, m_val = parse(Int, parts[1]), parse(Int, parts[2])
        total_min = h * 60 + m_val
        if total_min > 30
            walltime = "00:30:00"
            @warn "Debug QOS has a 30-minute limit. Clamping walltime to $walltime."
        end
    end
    root, platform = gacode_env(gpu)
    gacode_setup = """
    export GACODE_ROOT=$root
    export GACODE_PLATFORM=$platform
    . \${GACODE_ROOT}/shared/bin/gacode_setup 2>/dev/null
    . \${GACODE_ROOT}/platform/env/env.\${GACODE_PLATFORM} 2>/dev/null
    """

    # Step 1: Generate batch.src (without -s, so it doesn't submit yet)
    gen_cmd = """
    $gacode_setup
    cd $(Base.shell_escape(rundir))
    gacode_qsub -code cgyro -n $n_mpi -nomp $n_omp -queue $qos -repo $repo -p $(Base.shell_escape(rundir)) -w $walltime 2>&1
    """

    try
        read(`bash -c $gen_cmd`, String)
    catch e
        @warn "Failed to generate batch.src in $rundir: $e"
        return ""
    end

    batchfile = joinpath(rundir, "batch.src")
    if !isfile(batchfile)
        @warn "batch.src not found in $rundir after gacode_qsub"
        return ""
    end

    # Step 2: Inject GACODE environment into batch.src so the compute node
    # uses the correct GPU/CPU compiled GACODE (not whatever .bashrc.ext sets)
    gacode_env_lines = """
    export GACODE_ROOT=$root
    export GACODE_PLATFORM=$platform
    . \${GACODE_ROOT}/shared/bin/gacode_setup
    . \${GACODE_ROOT}/platform/env/env.\${GACODE_PLATFORM}
    """
    original = read(batchfile, String)
    # Insert after the last #SBATCH line and export SLURM_CPU_BIND line
    patched = replace(original, "export SLURM_CPU_BIND=\"cores\"\n" =>
                      "export SLURM_CPU_BIND=\"cores\"\n" * gacode_env_lines)
    write(batchfile, patched)

    # Step 3: Submit via sbatch
    output = try
        read(`sbatch $batchfile`, String)
    catch e
        @warn "Failed to sbatch $batchfile: $e"
        return ""
    end

    # Extract job ID from sbatch output (typically "Submitted batch job XXXXXXX")
    m = match(r"(\d{5,})", output)
    job_id = m !== nothing ? m[1] : ""

    if isempty(job_id)
        @info "gacode_qsub output: $output"
    end

    return job_id
end

"""
    is_interactive_node() -> Bool

Detect if running inside an interactive Slurm allocation (salloc session).
Returns true if SLURM_JOB_ID is set and SLURM_JOB_NAME is "interactive" or "bash".
"""
function is_interactive_node()
    job_id = get(ENV, "SLURM_JOB_ID", "")
    job_name = get(ENV, "SLURM_JOB_NAME", "")
    return !isempty(job_id) && job_name in ("interactive", "bash")
end

"""
    have_slurm() -> Bool

Return true if `sbatch` is on PATH (proxy for "we're on a slurm-managed cluster").
"""
function have_slurm()
    return Sys.which("sbatch") !== nothing
end

"""
    run_cgyro_local(rundirs::Vector{String}; n_mpi::Int=4, n_omp::Int=1) -> Nothing

Run linear CGYRO directly via the local `cgyro` wrapper script (`cgyro -n N -e <dir>`)
for each rundir, sequentially. Used when no slurm/srun is available — i.e. on a
workstation where GACODE is installed locally.

Requires `cgyro` to be on PATH (typically `\$GACODE_ROOT/cgyro/bin/cgyro`).
"""
function run_cgyro_local(rundirs::Vector{String}; n_mpi::Int=4, n_omp::Int=1)
    Sys.which("cgyro") === nothing &&
        error("`cgyro` not found on PATH; cannot run locally. Source the GACODE setup first.")

    for (i, rundir) in enumerate(rundirs)
        @info "Local CGYRO ($i/$(length(rundirs))): $(basename(rundir)) (n=$n_mpi, nomp=$n_omp)"
        logfile = joinpath(rundir, "cgyro.log")
        local_dir = basename(rundir)
        parent = dirname(rundir)
        # cgyro's `-e <dir>` resolves relative to the parent of the simdir, so cd one level up
        cmd = pipeline(
            Cmd(`cgyro -n $n_mpi -nomp $n_omp -e $local_dir`; dir=parent);
            stdout=logfile, stderr=logfile, append=false,
        )
        try
            run(cmd)
        catch e
            tail = try
                lines = readlines(logfile)
                join(lines[max(1, end - 30):end], "\n")
            catch
                "(no log)"
            end
            @error "Local CGYRO failed for $(basename(rundir)):\n$tail"
            rethrow(e)
        end
    end
    return nothing
end

"""
    interactive_node_count() -> Int

Return the number of nodes in the current interactive allocation.
"""
function interactive_node_count()
    return parse(Int, get(ENV, "SLURM_NNODES", "0"))
end

"""
    run_cgyro_interactive(rundirs::Vector{String}; gpu::Bool=true, n_mpi::Int=4, n_omp::Int=1) -> Vector{Process}

Launch CGYRO runs directly via srun on an interactive allocation.
Runs up to N jobs in parallel (one per allocated node), batching if there are more
rundirs than nodes. Blocks until all runs complete.

Typical linear CGYRO on Perlmutter GPU: n_mpi=4 (one per GPU), n_omp=1.
"""
function run_cgyro_interactive(rundirs::Vector{String}; gpu::Bool=true, n_mpi::Int=4, n_omp::Int=1)
    n_nodes = interactive_node_count()
    if n_nodes == 0
        error("Not on an interactive node (SLURM_NNODES not set)")
    end

    root, platform = gacode_env(gpu)
    gacode_setup = """
    export GACODE_ROOT=$root
    export GACODE_PLATFORM=$platform
    . \${GACODE_ROOT}/shared/bin/gacode_setup 2>/dev/null
    . \${GACODE_ROOT}/platform/env/env.\${GACODE_PLATFORM} 2>/dev/null
    export MPICH_GPU_SUPPORT_ENABLED=1
    export OMP_NUM_THREADS=$n_omp
    export OMP_STACKSIZE=400M
    """

    gpus_flag = gpu ? "--gpus-per-node=4" : ""

    # Process in batches of n_nodes
    for batch_start in 1:n_nodes:length(rundirs)
        batch_end = min(batch_start + n_nodes - 1, length(rundirs))
        batch = rundirs[batch_start:batch_end]

        @info "Interactive batch: running $(length(batch)) CGYRO jobs in parallel (nodes available: $n_nodes)"

        # Launch each run on a separate node
        processes = Process[]
        for (i, rundir) in enumerate(batch)
            node_offset = i - 1  # --relative is 0-indexed
            run_cmd = """
            $gacode_setup
            srun --relative=$node_offset --nodes=1 --ntasks=$n_mpi -c $n_omp --exact $gpus_flag \\
                \${GACODE_ROOT}/platform/exec/wrap.\${GACODE_PLATFORM} \\
                \${GACODE_ROOT}/cgyro/src/cgyro 0 \\
                > $(Base.shell_escape(joinpath(rundir, "batch.out"))) 2>&1
            """
            # Use input.cgyro from the run directory
            full_cmd = "cd $(Base.shell_escape(rundir)) && $run_cmd"
            proc = run(pipeline(`bash -c $full_cmd`); wait=false)
            push!(processes, proc)
            @info "  Node $i: KY dir $(basename(rundir)) (srun --relative=$node_offset)"
        end

        # Wait for this batch to finish
        for (i, proc) in enumerate(processes)
            wait(proc)
            status = proc.exitcode == 0 ? "success" : "exit code $(proc.exitcode)"
            @info "  Node $i ($(basename(batch[i]))): $status"
        end
    end

    return nothing
end

"""
    parse_cgyro_eigenvalue(rundir::String) -> Tuple{Float64, Float64}

Parse the converged growth rate and frequency from out.cgyro.freq.
Returns (frequency, growth_rate). Growth rate < 0 means stable.
"""
function parse_cgyro_eigenvalue(rundir::String)
    freqfile = joinpath(rundir, "out.cgyro.freq")
    if !isfile(freqfile)
        return (0.0, 0.0)
    end

    lines = readlines(freqfile)
    if isempty(lines)
        return (0.0, 0.0)
    end

    # Last line has the converged eigenvalue: freq growth_rate
    parts = split(strip(lines[end]))
    if length(parts) >= 2
        freq = parse(Float64, parts[1])
        gamma = parse(Float64, parts[2])
        return (freq, gamma)
    end

    return (0.0, 0.0)
end

"""
    parse_cgyro_qlflux(rundir::String, n_species::Int, n_field::Int) -> Array{Float32, 3}

Parse the QL flux from bin.cgyro.ky_flux.
Returns the converged (last time step) array of shape (n_species, n_field, 3_moments).
Moments: (particle, energy, momentum).
"""
function parse_cgyro_qlflux(rundir::String, n_species::Int, n_field::Int)
    fluxfile = joinpath(rundir, "bin.cgyro.ky_flux")
    if !isfile(fluxfile)
        return zeros(Float32, n_species, n_field, 3)
    end

    data = Vector{Float32}(undef, filesize(fluxfile) ÷ 4)
    read!(fluxfile, data)

    n_per_step = n_species * n_field * 3
    n_time = length(data) ÷ n_per_step
    if n_time == 0
        return zeros(Float32, n_species, n_field, 3)
    end

    # CGYRO Fortran writes gflux(0, :, 1:nflux, :, :) in column-major order:
    #   (n_species, n_flux_types, n_field, n_time) — species varies fastest
    all_data = reshape(data, n_species, 3, n_field, n_time)
    # Last time step → (n_species, 3_moments, n_field), reorder to (n_species, n_field, 3_moments)
    return permutedims(all_data[:, :, :, end], (1, 3, 2))
end

"""
    run_qlgyro(input_cgyro::InputCGYRO, input_qlgyro::InputQLGYRO;
               basedir::String="",
               gpu::Bool=true,
               n_mpi::Int=32,
               n_omp::Int=4,
               walltime::String="00:15:00",
               repo::String="m3739_g",
               qos::String="regular",
               ky_indices::Union{Nothing, Vector{Int}}=nothing,
               wait_for_completion::Bool=true,
               poll_interval::Int=60) -> Union{GACODE.FluxSolution, QLGYRORunState}

Run modular QLGYRO from an `(InputCGYRO, InputQLGYRO)` pair: linear CGYRO at each ky in
the QLGYRO ky scan + TJLF saturation rules.

Workflow:
1. Generates ky grid from `input_qlgyro` (KYGRID_MODEL/NKY/KY)
2. Creates run directories, one per ky
3. If on an interactive node (salloc), runs CGYRO directly via srun (N jobs in parallel,
   one per node). Otherwise, submits batch jobs via gacode_qsub.
4. Monitors convergence and resubmits if needed
5. Parses CGYRO growth rates and QL weights
6. Applies TJLF saturation rules (sum_ky_spectrum)
7. Returns FluxSolution

If `wait_for_completion=false`, returns the run state for later checking with `refresh_qlgyro!`.
Set `ky_indices` to run only specific ky points (1-indexed) for testing.
"""
function run_qlgyro(input_cgyro::InputCGYRO, input_qlgyro::InputQLGYRO;
                    basedir::String="",
                    gpu::Bool=true,
                    n_mpi::Int=32,
                    n_omp::Int=4,
                    walltime::String="00:15:00",
                    repo::String="m3739_g",
                    qos::String="regular",
                    ky_indices::Union{Nothing,Vector{Int}}=nothing,
                    wait_for_completion::Bool=true,
                    poll_interval::Int=60,
                    local_mode::Union{Bool,Nothing}=nothing)

    input_hash = qlgyro_run_hash(input_cgyro)

    ky_grid = generate_ky_grid(input_cgyro, input_qlgyro)

    if isempty(basedir)
        preferred = "/global/cfs/cdirs/m3739/results/FUSE/QLGYRO"
        if isdir(preferred) || (try mkpath(preferred); true catch; false end)
            tmproot = preferred
        else
            scratch = get(ENV, "PSCRATCH", get(ENV, "SCRATCH", ""))
            tmproot = isempty(scratch) ? tempdir() : scratch
        end
        username = get(ENV, "USER", get(ENV, "USERNAME", "user"))
        timestamp = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
        basedir = joinpath(tmproot, "qlgyro_$(username)_$(timestamp)")
        mkpath(basedir)
        @info "QLGYRO run directory: $basedir"
    end
    mkpath(basedir)

    state = load_run_state(basedir)
    if state !== nothing && state.input_hash == input_hash
        state.basedir = basedir
        @info "Found existing QLGYRO run state. Checking status..."
    else
        nky = length(ky_grid)
        state = QLGYRORunState(basedir, ky_grid, fill("", nky), fill(false, nky), fill(false, nky), input_hash)
    end

    save(input_cgyro, joinpath(basedir, "input.cgyro"))
    save(input_qlgyro, joinpath(basedir, "input.qlgyro"))

    active_indices = ky_indices !== nothing ? ky_indices : collect(1:length(ky_grid))

    rundirs_to_launch = String[]

    for iky in active_indices
        ky_dir = joinpath(basedir, "KY_$iky")
        mkpath(ky_dir)

        if state.converged[iky]
            @info "KY_$iky already converged, skipping."
            continue
        end

        ic_ky = deepcopy(input_cgyro)
        ic_ky.KY = ky_grid[iky]
        # Per-ky inputs embed all species data directly, so use PROFILE_MODEL=1
        # (PROFILE_MODEL=2 would require an input.gacode profiles file alongside).
        ic_ky.PROFILE_MODEL = 1
        save(ic_ky, joinpath(ky_dir, "input.cgyro"))

        status = check_cgyro_convergence(ky_dir)
        if status == :converged
            state.converged[iky] = true
            @info "KY_$iky (ky=$(round(ky_grid[iky], digits=4))): already converged"
            continue
        end

        if !isempty(state.slurm_ids[iky])
            slurm_status = check_slurm_status(state.slurm_ids[iky])
            if slurm_status in (:pending, :running)
                @info "KY_$iky (ky=$(round(ky_grid[iky], digits=4))): job $(state.slurm_ids[iky]) still $slurm_status"
                continue
            end
        end

        push!(rundirs_to_launch, ky_dir)
    end

    interactive = is_interactive_node()
    use_local = local_mode === nothing ? !interactive && !have_slurm() : local_mode

    if interactive && !isempty(rundirs_to_launch)
        n_nodes = interactive_node_count()
        @info "Interactive mode detected ($n_nodes nodes). Running $(length(rundirs_to_launch)) CGYRO jobs directly via srun."
        run_cgyro_interactive(rundirs_to_launch; gpu=gpu, n_mpi=min(n_mpi, 4), n_omp=1)

        for ky_dir in rundirs_to_launch
            iky = parse(Int, match(r"KY_(\d+)", basename(ky_dir))[1])
            state.submitted[iky] = true
            status = check_cgyro_convergence(ky_dir)
            if status == :converged
                state.converged[iky] = true
                @info "KY_$iky: converged!"
            else
                @warn "KY_$iky: finished with status=$status"
            end
        end
    elseif use_local && !isempty(rundirs_to_launch)
        @info "Local mode: running $(length(rundirs_to_launch)) CGYRO jobs sequentially via local `cgyro`."
        run_cgyro_local(rundirs_to_launch; n_mpi=n_mpi, n_omp=n_omp)

        for ky_dir in rundirs_to_launch
            iky = parse(Int, match(r"KY_(\d+)", basename(ky_dir))[1])
            state.submitted[iky] = true
            status = check_cgyro_convergence(ky_dir)
            if status == :converged
                state.converged[iky] = true
                @info "KY_$iky: converged!"
            else
                @warn "KY_$iky: finished with status=$status"
            end
        end
    else
        for ky_dir in rundirs_to_launch
            iky = parse(Int, match(r"KY_(\d+)", basename(ky_dir))[1])
            job_id = submit_cgyro_job(ky_dir; gpu=gpu, n_mpi=n_mpi, n_omp=n_omp, walltime=walltime, repo=repo, qos=qos)
            state.slurm_ids[iky] = job_id
            state.submitted[iky] = true
            @info "KY_$iky (ky=$(round(ky_grid[iky], digits=4))): submitted job $job_id"
        end
    end

    save_run_state(state)

    if !wait_for_completion
        return state
    end

    # In local/interactive modes the runs are synchronous; skip the slurm poll loop.
    if !(use_local || interactive)
        @info "Waiting for all CGYRO jobs to complete..."
        max_resubmits = 3
        resubmit_count = zeros(Int, length(ky_grid))

        while !all(state.converged[i] for i in active_indices)
            sleep(poll_interval)
            refresh_qlgyro!(state; active_indices=active_indices, gpu=gpu, n_mpi=n_mpi, n_omp=n_omp,
                            walltime=walltime, repo=repo, max_resubmits=max_resubmits,
                            resubmit_count=resubmit_count)
        end
    end

    return compute_qlgyro_fluxes(input_cgyro, input_qlgyro, state;
                                 gbflux_path=joinpath(basedir, "out.qlgyro.gbflux"))
end

"""
    run_qlgyro(input_tglf::InputTGLF; kwargs...)

Backward-compatible wrapper: builds an `InputCGYRO` (via `tglf_to_cgyro`) and an
`InputQLGYRO` (via `qlgyro_from_tglf`) from the TGLF inputs, then dispatches to
`run_qlgyro(::InputCGYRO, ::InputQLGYRO; ...)`. Pass `kygrid_model` here to override
the QLGYRO `KYGRID_MODEL` for the ky scan.
"""
function run_qlgyro(input_tglf::InputTGLF; kygrid_model::Int=-1, kwargs...)
    input_cgyro = tglf_to_cgyro(input_tglf)
    input_qlgyro = qlgyro_from_tglf(input_tglf)
    if kygrid_model >= 0
        input_qlgyro.KYGRID_MODEL = kygrid_model
    end
    return run_qlgyro(input_cgyro, input_qlgyro; kwargs...)
end

"""
    refresh_qlgyro!(state::QLGYRORunState; kwargs...) -> QLGYRORunState

Check convergence of all KY runs and resubmit failed ones.
"""
function refresh_qlgyro!(state::QLGYRORunState;
                         active_indices::Vector{Int}=collect(1:length(state.ky_values)),
                         gpu::Bool=true, n_mpi::Int=32, n_omp::Int=4,
                         walltime::String="00:15:00", repo::String="m3739_g",
                         qos::String="regular",
                         max_resubmits::Int=3,
                         resubmit_count::Vector{Int}=zeros(Int, length(state.ky_values)))
    n_converged = 0
    n_running = 0
    n_pending = 0
    n_resubmitted = 0
    n_failed = 0

    for iky in active_indices
        if state.converged[iky]
            n_converged += 1
            continue
        end

        ky_dir = joinpath(state.basedir, "KY_$iky")
        status = check_cgyro_convergence(ky_dir)

        if status == :converged
            state.converged[iky] = true
            n_converged += 1
            @info "KY_$iky: converged!"
            continue
        end

        # Check Slurm status
        if !isempty(state.slurm_ids[iky])
            slurm_status = check_slurm_status(state.slurm_ids[iky])
            if slurm_status == :running
                n_running += 1
                continue
            elseif slurm_status == :pending
                n_pending += 1
                continue
            end
        end

        # Job finished but not converged - resubmit if under limit
        if status in (:terminated, :error, :not_started, :running) && resubmit_count[iky] < max_resubmits
            resubmit_count[iky] += 1
            reason = status == :error ? "CGYRO error" :
                     status == :terminated ? "hit max time without converging" :
                     status == :not_started ? "no output (job may have failed to start)" :
                     "job finished but CGYRO still running (likely walltime exceeded)"
            job_id = submit_cgyro_job(ky_dir; gpu=gpu, n_mpi=n_mpi, n_omp=n_omp, walltime=walltime, repo=repo, qos=qos)
            state.slurm_ids[iky] = job_id
            n_resubmitted += 1
            @info "KY_$iky: resubmitting (attempt $(resubmit_count[iky])/$max_resubmits), job $job_id — reason: $reason"
        elseif resubmit_count[iky] >= max_resubmits
            @warn "KY_$iky: max resubmits reached, marking as converged (best effort)"
            state.converged[iky] = true
            n_failed += 1
        end
    end

    status_parts = ["$n_converged converged", "$n_running running", "$n_pending pending"]
    if n_resubmitted > 0
        push!(status_parts, "$n_resubmitted resubmitted")
    end
    if n_failed > 0
        push!(status_parts, "$n_failed gave up (best effort)")
    end
    @info "Status: $(join(status_parts, ", ")) out of $(length(active_indices))"
    save_run_state(state)
    return state
end

"""
    compute_qlgyro_fluxes(input_cgyro::InputCGYRO, input_qlgyro::InputQLGYRO, state::QLGYRORunState;
                          sat_rule=nothing, alpha_zf=nothing) -> GACODE.FluxSolution

Parse CGYRO outputs from each `KY_*` directory in `state` and apply TJLF saturation rules
to produce quasilinear fluxes. Builds an `InputTJLF` directly from the (CGYRO, QLGYRO) pair
via `build_input_tjlf` — no `InputTGLF` round-trip.
"""
function compute_qlgyro_fluxes(input_cgyro::InputCGYRO, input_qlgyro::InputQLGYRO, state::QLGYRORunState;
                               sat_rule::Union{Int,Nothing}=nothing, alpha_zf::Union{Real,Nothing}=nothing,
                               gbflux_path::Union{Nothing,AbstractString}=nothing)
    ns = input_cgyro.N_SPECIES
    nky = length(state.ky_values)
    nmodes = 1  # linear CGYRO gives dominant mode only
    n_field = 3  # Phi, Apar, Bpar (or fewer, padded with zeros)

    # Parse growth rates and frequencies from CGYRO
    gamma_cgyro = zeros(Float64, nky)  # growth rate per ky
    freq_cgyro = zeros(Float64, nky)   # frequency per ky

    # QL weights from CGYRO: (nky, ns_cgyro, n_field, 3_moments)
    # CGYRO species order: ion1, ion2, ..., electrons
    ql_flux_cgyro = zeros(Float64, nky, ns, n_field, 3)

    for iky in 1:nky
        ky_dir = joinpath(state.basedir, "KY_$iky")

        # Parse eigenvalue
        freq, gamma = parse_cgyro_eigenvalue(ky_dir)
        gamma_cgyro[iky] = max(gamma, 0.0)  # set negative (stable) growth rates to 0
        freq_cgyro[iky] = freq

        ql_raw = parse_cgyro_qlflux(ky_dir, ns, n_field)
        # CGYRO convention: qlflux_ky = flux * ky
        ql_flux_cgyro[iky, :, :, :] .= Float64.(ql_raw) .* state.ky_values[iky]
    end

    # Reorder species CGYRO [ion1..ion(N-1), e] -> TJLF [e, ion1..ion(N-1)]
    ql_flux_tjlf = zeros(Float64, nky, ns, n_field, 3)
    ql_flux_tjlf[:, 1, :, :] .= ql_flux_cgyro[:, end, :, :]  # electrons
    for i in 2:ns
        ql_flux_tjlf[:, i, :, :] .= ql_flux_cgyro[:, i-1, :, :]  # ions
    end

    # Build InputTJLF directly from (CGYRO, QLGYRO) with the ky scan pre-populated
    input_tjlf = build_input_tjlf(input_cgyro, input_qlgyro;
                                  ky_values=state.ky_values, nmodes=nmodes)

    # Apply optional overrides
    effective_sat_rule = something(sat_rule, input_tjlf.SAT_RULE)
    if sat_rule !== nothing
        input_tjlf.SAT_RULE = sat_rule
    end
    if alpha_zf !== nothing
        input_tjlf.ALPHA_ZF = Float64(alpha_zf)
    end

    satParams = TJLF.get_sat_params(input_tjlf)

    gamma_matrix = zeros(Float64, nmodes, nky)
    gamma_matrix[1, :] .= gamma_cgyro

    # QL_weights for TJLF: shape (nf, ns, nm, nky, ntype)
    # ntype: (1=particle, 2=energy, 3=stress_tor, 4=stress_par, 5=exchange)
    QL_weights = zeros(Float64, n_field, ns, nmodes, nky, 5)
    for iky in 1:nky
        for is in 1:ns
            for ifield in 1:n_field
                QL_weights[ifield, is, 1, iky, 1] = ql_flux_tjlf[iky, is, ifield, 1]
                QL_weights[ifield, is, 1, iky, 2] = ql_flux_tjlf[iky, is, ifield, 2]
                QL_weights[ifield, is, 1, iky, 3] = ql_flux_tjlf[iky, is, ifield, 3]
            end
        end
    end

    if effective_sat_rule in (2, 3)
        vzf_out, kymax_out, jmax_out = TJLF.get_zonal_mixing(input_tjlf, satParams, gamma_cgyro)
        QL_flux_out, _ = TJLF.sum_ky_spectrum(input_tjlf, satParams, gamma_matrix, QL_weights;
                                              vzf_out_param=vzf_out,
                                              kymax_out_param=kymax_out,
                                              jmax_out_param=jmax_out)
    else
        QL_flux_out, _ = TJLF.sum_ky_spectrum(input_tjlf, satParams, gamma_matrix, QL_weights)
    end

    if gbflux_path !== nothing
        write_gbflux(QL_flux_out, gbflux_path)
    end

    return GACODE.FluxSolution{Float64}(
        TJLF.Qe(QL_flux_out),
        TJLF.Qi(QL_flux_out),
        TJLF.Γe(QL_flux_out),
        TJLF.Γi(QL_flux_out),
        TJLF.Πi(QL_flux_out)
    )
end

"""
    write_gbflux(QL_flux_out::Array, filename::AbstractString)

Write fluxes to `out.qlgyro.gbflux` in TGLF/QLGYRO format: four whitespace-separated
rows (Γ, Q, Π, S), each with NS values (electron first, then ions). The S (exchange)
row is written as zeros since linear CGYRO does not provide it.

`QL_flux_out` shape: (n_field, ns, ntype) — sums over fields and selects ntype.
"""
function write_gbflux(QL_flux_out::AbstractArray, filename::AbstractString)
    ns = size(QL_flux_out, 2)
    open(filename, "w") do io
        # Γ (particle): ntype=1
        for is in 1:ns
            print(io, sum(QL_flux_out[:, is, 1]), " ")
        end
        println(io)
        # Q (energy): ntype=2
        for is in 1:ns
            print(io, sum(QL_flux_out[:, is, 2]), " ")
        end
        println(io)
        # Π (stress_tor): ntype=3
        for is in 1:ns
            print(io, sum(QL_flux_out[:, is, 3]), " ")
        end
        println(io)
        # S (exchange): zeros (not from linear CGYRO)
        for _ in 1:ns
            print(io, "0.0 ")
        end
        println(io)
    end
    return filename
end

# Backward-compatibility shim for any caller that still passes an InputTGLF.
function compute_qlgyro_fluxes(input_tglf::InputTGLF, state::QLGYRORunState; kwargs...)
    input_cgyro = tglf_to_cgyro(input_tglf)
    input_qlgyro = qlgyro_from_tglf(input_tglf)
    return compute_qlgyro_fluxes(input_cgyro, input_qlgyro, state; kwargs...)
end

"""
    qlgyro_from_tglf(input_tglf::InputTGLF) -> InputQLGYRO

Build an `InputQLGYRO` from the QLGYRO-relevant fields of an `InputTGLF`
(KYGRID_MODEL/NKY/KY/SAT_RULE and the rotation shear).
"""
function qlgyro_from_tglf(input_tglf::InputTGLF)
    input_qlgyro = InputQLGYRO()
    input_qlgyro.NKY = input_tglf.NKY
    input_qlgyro.KYGRID_MODEL = input_tglf.KYGRID_MODEL
    if !ismissing(input_tglf.KY)
        input_qlgyro.KY = Float64(input_tglf.KY)
    end
    input_qlgyro.SAT_RULE = input_tglf.SAT_RULE
    input_qlgyro.GAMMA_E = ismissing(input_tglf.VEXB_SHEAR) ? 0.0 : Float64(input_tglf.VEXB_SHEAR)
    return input_qlgyro
end

export run_qlgyro
