using Serialization

Base.@kwdef mutable struct InputQLGYRO
    N_PARALLEL::Union{Int,Missing} = missing
    N_RUNS::Union{Int,Missing} = missing
    GAMMA_E::Union{Float64,Missing} = missing
    CODE::Union{Int,Missing} = missing
    NKY::Union{Int,Missing} = missing
    KYGRID_MODEL::Union{Int,Missing} = missing
    KY::Union{Float64,Missing} = missing
    SAT_RULE::Union{Int,Missing} = missing
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

"""
    check_slurm_status(job_id::String) -> Symbol

Check if a Slurm job is still in queue/running.
Returns :pending, :running, :completed, or :unknown.
"""
function check_slurm_status(job_id::String)
    if isempty(job_id)
        return :unknown
    end
    try
        output = strip(read(`sacct -j $job_id --format=State --noheader -P`, String))
        lines = filter(!isempty, split(output, '\n'))
        if isempty(lines)
            return :unknown
        end
        state = strip(lines[1])
        if state == "PENDING"
            return :pending
        elseif state == "RUNNING"
            return :running
        elseif state in ("COMPLETED", "FAILED", "TIMEOUT", "CANCELLED")
            return :completed
        else
            return :unknown
        end
    catch
        return :unknown
    end
end

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

    # Reshape: (n_time, n_species, n_field, 3_moments)
    all_data = reshape(data, 3, n_field, n_species, n_time)
    # Last time step, reorder to (n_species, n_field, 3_moments)
    return permutedims(all_data[:, :, :, end], (3, 2, 1))
end

"""
    run_qlgyro(input_tglf::InputTGLF;
               basedir::String="",
               gpu::Bool=true,
               n_mpi::Int=32,
               n_omp::Int=4,
               walltime::String="00:15:00",
               repo::String="m3739_g",
               ky_indices::Union{Nothing, Vector{Int}}=nothing,
               kygrid_model::Int=-1,
               wait_for_completion::Bool=true,
               poll_interval::Int=60) -> Union{GACODE.FluxSolution, QLGYRORunState}

Run modular QLGYRO: linear CGYRO at each ky + TJLF saturation rules.

Workflow:
1. Converts InputTGLF to InputCGYRO
2. Generates ky grid (default: input_tglf.KYGRID_MODEL; use `kygrid_model=0` for simple linear grid with KY=1.2)
3. Creates run directories, one per ky
4. If on an interactive node (salloc), runs CGYRO directly via srun (N jobs in parallel, one per node).
   Otherwise, submits batch jobs via gacode_qsub.
5. Monitors convergence and resubmits if needed
6. Parses CGYRO growth rates and QL weights
7. Applies TJLF saturation rules (sum_ky_spectrum)
8. Returns FluxSolution

If `wait_for_completion=false`, returns the run state for later checking with `refresh_qlgyro!`.
Set `ky_indices` to run only specific ky points (1-indexed) for testing.
"""
function run_qlgyro(input_tglf::InputTGLF;
                    basedir::String="",
                    gpu::Bool=true,
                    n_mpi::Int=32,
                    n_omp::Int=4,
                    walltime::String="00:15:00",
                    repo::String="m3739_g",
                    qos::String="regular",
                    ky_indices::Union{Nothing, Vector{Int}}=nothing,
                    kygrid_model::Int=-1,
                    wait_for_completion::Bool=true,
                    poll_interval::Int=60)

    # Step 1: Convert InputTGLF -> InputCGYRO
    input_cgyro = tglf_to_cgyro(input_tglf)
    input_hash = qlgyro_run_hash(input_cgyro)

    # Step 2: Generate ky grid
    ky_grid = generate_ky_grid(input_tglf; kygrid_model=kygrid_model)

    # Step 3: Setup run directory (must be on a shared filesystem, not /tmp which is node-local)
    if isempty(basedir)
        scratch = get(ENV, "PSCRATCH", get(ENV, "SCRATCH", ""))
        tmproot = isempty(scratch) ? tempdir() : scratch
        basedir = mktempdir(tmproot; prefix="qlgyro_", cleanup=false)
        @info "QLGYRO run directory: $basedir"
    end
    mkpath(basedir)

    # Check for existing state
    state = load_run_state(basedir)
    if state !== nothing && state.input_hash == input_hash
        @info "Found existing QLGYRO run state. Checking status..."
    else
        # Initialize new state
        nky = length(ky_grid)
        state = QLGYRORunState(basedir, ky_grid, fill("", nky), fill(false, nky), fill(false, nky), input_hash)
    end

    # Save base input.cgyro for reference
    save(input_cgyro, joinpath(basedir, "input.cgyro"))

    # Also save input.tglf for the saturation rules
    save(input_tglf, joinpath(basedir, "input.tglf"))

    # Step 4: Create KY directories and prepare input files
    active_indices = ky_indices !== nothing ? ky_indices : collect(1:length(ky_grid))

    # Collect directories that need running (not converged, not already in-flight)
    rundirs_to_launch = String[]

    for iky in active_indices
        ky_dir = joinpath(basedir, "KY_$iky")
        mkpath(ky_dir)

        # Skip if already converged
        if state.converged[iky]
            @info "KY_$iky already converged, skipping."
            continue
        end

        # Write input.cgyro with this KY value
        ic_ky = deepcopy(input_cgyro)
        ic_ky.KY = ky_grid[iky]
        save(ic_ky, joinpath(ky_dir, "input.cgyro"))

        # Check if there's an existing converged run
        status = check_cgyro_convergence(ky_dir)
        if status == :converged
            state.converged[iky] = true
            @info "KY_$iky (ky=$(round(ky_grid[iky], digits=4))): already converged"
            continue
        end

        # Check if job is still running (batch mode only)
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

    if interactive && !isempty(rundirs_to_launch)
        # Interactive mode: run CGYRO directly via srun on allocated nodes
        n_nodes = interactive_node_count()
        @info "Interactive mode detected ($n_nodes nodes). Running $(length(rundirs_to_launch)) CGYRO jobs directly via srun."
        # Linear CGYRO: 4 MPI ranks (one per GPU), 1 OpenMP thread
        run_cgyro_interactive(rundirs_to_launch; gpu=gpu, n_mpi=min(n_mpi, 4), n_omp=1)

        # Mark all launched runs as submitted and check convergence
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
        # Batch mode: submit via gacode_qsub
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

    # Step 5: Wait for completion
    @info "Waiting for all CGYRO jobs to complete..."
    max_resubmits = 3
    resubmit_count = zeros(Int, length(ky_grid))

    while !all(state.converged[i] for i in active_indices)
        sleep(poll_interval)
        refresh_qlgyro!(state; active_indices=active_indices, gpu=gpu, n_mpi=n_mpi, n_omp=n_omp,
                        walltime=walltime, repo=repo, max_resubmits=max_resubmits,
                        resubmit_count=resubmit_count)
    end

    # Step 6: Compute fluxes
    return compute_qlgyro_fluxes(input_tglf, state)
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
    compute_qlgyro_fluxes(input_tglf::InputTGLF, state::QLGYRORunState) -> GACODE.FluxSolution

Parse CGYRO outputs and apply TJLF saturation rules to compute quasilinear fluxes.
"""
function compute_qlgyro_fluxes(input_tglf::InputTGLF, state::QLGYRORunState)
    ns = input_tglf.NS  # total species (TGLF ordering: e, ion1, ion2, ...)
    nky = length(state.ky_values)
    nmodes = 1  # linear CGYRO gives dominant mode only
    n_field = 3  # Phi, Apar, Bpar (or fewer, padded with zeros)

    # Parse growth rates and frequencies from CGYRO
    gamma_cgyro = zeros(Float64, nky)  # growth rate per ky
    freq_cgyro = zeros(Float64, nky)   # frequency per ky

    # QL weights from CGYRO: (nky, ns_cgyro, n_field, 3_moments)
    # CGYRO species order: ion1, ion2, ..., electrons
    ns_cgyro = ns
    ql_flux_cgyro = zeros(Float64, nky, ns_cgyro, n_field, 3)

    for iky in 1:nky
        ky_dir = joinpath(state.basedir, "KY_$iky")

        # Parse eigenvalue
        freq, gamma = parse_cgyro_eigenvalue(ky_dir)
        gamma_cgyro[iky] = max(gamma, 0.0)  # set negative (stable) growth rates to 0
        freq_cgyro[iky] = freq

        # Parse QL flux
        ql_raw = parse_cgyro_qlflux(ky_dir, ns_cgyro, n_field)
        # ql_raw shape: (ns_cgyro, n_field, 3_moments)
        # Multiply by ky to get the QL weight (CGYRO convention: qlflux_ky = flux * ky)
        ql_flux_cgyro[iky, :, :, :] .= Float64.(ql_raw) .* state.ky_values[iky]
    end

    # Reorder species from CGYRO (ions first, electrons last) to TGLF/TJLF (electrons first, ions after)
    # CGYRO: [ion1, ion2, ..., e] -> TJLF: [e, ion1, ion2, ...]
    ql_flux_tjlf = zeros(Float64, nky, ns, n_field, 3)
    ql_flux_tjlf[:, 1, :, :] .= ql_flux_cgyro[:, end, :, :]  # electrons
    for i in 2:ns
        ql_flux_tjlf[:, i, :, :] .= ql_flux_cgyro[:, i-1, :, :]  # ions
    end

    # Build InputTJLF for saturation rules with the correct nky size
    # We construct directly with the right nky rather than resizing after
    input_tjlf = InputTJLF{Float64}(input_tglf.NS, nky)
    TJLF.update_input_tjlf!(input_tjlf, input_tglf)

    # Override ky grid settings
    # We set KYGRID_MODEL=0 because KY_SPECTRUM is pre-populated from the CGYRO runs;
    # sum_ky_spectrum uses KY_SPECTRUM directly and never reads KYGRID_MODEL or NKY.
    input_tjlf.KYGRID_MODEL = 0
    input_tjlf.NKY = nky
    input_tjlf.NMODES = nmodes
    input_tjlf.KY_SPECTRUM .= state.ky_values
    input_tjlf.WIDTH_SPECTRUM .= input_tjlf.WIDTH

    # Compute saturation parameters (geometry)
    satParams = TJLF.get_sat_params(input_tjlf)

    # Build gamma_matrix for TJLF: shape (nmodes, nky)
    gamma_matrix = zeros(Float64, nmodes, nky)
    gamma_matrix[1, :] .= gamma_cgyro

    # Build QL_weights for TJLF: shape (nf, ns, nm, nky, ntype)
    # ntype: (1=particle, 2=energy, 3=stress_tor, 4=stress_par, 5=exchange)
    # CGYRO moments: (1=particle, 2=energy, 3=momentum)
    QL_weights = zeros(Float64, n_field, ns, nmodes, nky, 5)

    for iky in 1:nky
        for is in 1:ns
            for ifield in 1:n_field
                QL_weights[ifield, is, 1, iky, 1] = ql_flux_tjlf[iky, is, ifield, 1]  # particle
                QL_weights[ifield, is, 1, iky, 2] = ql_flux_tjlf[iky, is, ifield, 2]  # energy
                QL_weights[ifield, is, 1, iky, 3] = ql_flux_tjlf[iky, is, ifield, 3]  # stress_tor
                # stress_par and exchange not available from CGYRO linear, leave as 0
            end
        end
    end

    # For ExB shear spectral shift, compute first pass eigenvalues for SAT2/3
    vexb_shear_s = input_tjlf.VEXB_SHEAR * input_tjlf.SIGN_IT
    if input_tglf.SAT_RULE in (2, 3)
        vzf_out, kymax_out, jmax_out = TJLF.get_zonal_mixing(input_tjlf, satParams, gamma_cgyro)
        QL_flux_out, flux_spectrum = TJLF.sum_ky_spectrum(input_tjlf, satParams, gamma_matrix, QL_weights;
                                                           vzf_out_param=vzf_out,
                                                           kymax_out_param=kymax_out,
                                                           jmax_out_param=jmax_out)
    else
        QL_flux_out, flux_spectrum = TJLF.sum_ky_spectrum(input_tjlf, satParams, gamma_matrix, QL_weights)
    end

    # Extract fluxes using TJLF helpers
    T = Float64
    sol = GACODE.FluxSolution{T}(
        TJLF.Qe(QL_flux_out),
        TJLF.Qi(QL_flux_out),
        TJLF.Γe(QL_flux_out),
        TJLF.Γi(QL_flux_out),
        TJLF.Πi(QL_flux_out)
    )

    return sol
end

export run_qlgyro
