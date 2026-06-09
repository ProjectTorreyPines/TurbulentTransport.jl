# Generic SLURM submit/poll helpers shared by the QLGYRO (qlgyro.jl) and TGLF-EP
# (tglf_ep.jl) "from dd" run drivers. Kept dependency-free (Serialization + base shell-out)
# so both runners reuse one implementation instead of carrying private copies.

using Serialization
import Dates

"""
    check_slurm_status(job_id::AbstractString) -> Symbol

Query `sacct` for a job's state. Returns one of
`:pending`, `:running`, `:completed`, `:failed`, or `:unknown`.
A job that finished cleanly is `:completed`; `FAILED`/`TIMEOUT`/`CANCELLED`/`NODE_FAIL`
map to `:failed` so callers can distinguish success from failure.
"""
function check_slurm_status(job_id::AbstractString)
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
        elseif state in ("RUNNING", "CONFIGURING", "COMPLETING")
            return :running
        elseif state == "COMPLETED"
            return :completed
        elseif state in ("FAILED", "TIMEOUT", "CANCELLED", "NODE_FAIL", "OUT_OF_MEMORY", "BOOT_FAIL", "DEADLINE") ||
               startswith(state, "CANCELLED")
            return :failed
        else
            return :unknown
        end
    catch
        return :unknown
    end
end

"""
    parse_slurm_jobid(sbatch_output::AbstractString) -> String

Extract the numeric job id from `sbatch` output (e.g. "Submitted batch job 12345678").
Returns "" if no id is found.
"""
function parse_slurm_jobid(sbatch_output::AbstractString)
    m = match(r"(\d{5,})", sbatch_output)
    return m !== nothing ? String(m[1]) : ""
end

"""
    default_results_dir(subdir::AbstractString;
                        runprefix::AbstractString=lowercase(subdir),
                        prefix::AbstractString="/global/cfs/cdirs/m3739/results/FUSE") -> String

Build a per-run results directory `<prefix>/<subdir>/<runprefix>_<user>_<timestamp>`,
matching the QLGYRO naming convention (e.g. `.../FUSE/TJLFEP/tjlfep_<user>_<timestamp>`).
Falls back to `\$PSCRATCH`/`tempdir()` if the CFS prefix is not writable.
"""
function default_results_dir(subdir::AbstractString;
                             runprefix::AbstractString=lowercase(subdir),
                             prefix::AbstractString="/global/cfs/cdirs/m3739/results/FUSE")
    base = joinpath(prefix, subdir)
    if !(isdir(base) || (try mkpath(base); true catch; false end))
        scratch = get(ENV, "PSCRATCH", get(ENV, "SCRATCH", ""))
        base = isempty(scratch) ? tempdir() : joinpath(scratch, "results", "FUSE", subdir)
        mkpath(base)
    end
    username = get(ENV, "USER", get(ENV, "USERNAME", "user"))
    timestamp = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
    return joinpath(base, "$(runprefix)_$(username)_$(timestamp)")
end

"""
    save_state(state, path::AbstractString)

Serialize an arbitrary run-state object to `path` (Julia `.jls`).
"""
function save_state(state, path::AbstractString)
    mkpath(dirname(path))
    Serialization.serialize(path, state)
    return path
end

"""
    load_state(path::AbstractString)

Deserialize a run-state object previously written by [`save_state`](@ref), or
`nothing` if `path` does not exist.
"""
function load_state(path::AbstractString)
    return isfile(path) ? Serialization.deserialize(path) : nothing
end
