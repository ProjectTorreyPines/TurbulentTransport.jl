#= ======================== =#
#  Model path resolution API
#= ======================== =#

using Downloads
using SHA
using Scratch

const _MODEL_SEARCH_PATHS = String[]  # Additional search paths from providers
const _LFS_POINTER_PREFIX = "version https://git-lfs.github.com/spec/v1"
const _LFS_REPO = "ProjectTorreyPines/TurbulentTransport.jl"
const _model_download_lock = ReentrantLock()

function _pkg_root()
    return normpath(joinpath(@__DIR__, ".."))
end

"""
    is_lfs_pointer(path::AbstractString) -> Bool

Return `true` when `path` points to a Git LFS stub instead of real model bytes.
"""
function is_lfs_pointer(path::AbstractString)
    isfile(path) || return false
    sz = filesize(path)
    sz > 500 && return false
    sz >= sizeof(_LFS_POINTER_PREFIX) || return false
    bytes = read(path, sizeof(_LFS_POINTER_PREFIX))
    return String(bytes) == _LFS_POINTER_PREFIX
end

# Parse a Git LFS pointer file. Returns `(oid, size)` or `nothing` if the file
# is not a recognizable pointer. `oid` is the lowercase SHA-256 hex string of
# the real content; `size` is the expected byte length.
function _lfs_pointer_info(path::AbstractString)
    is_lfs_pointer(path) || return nothing
    content = read(path, String)
    oid_match = match(r"oid sha256:([0-9a-fA-F]{64})"i, content)
    size_match = match(r"size (\d+)"i, content)
    (oid_match === nothing || size_match === nothing) && return nothing
    return (oid=lowercase(oid_match.captures[1]), size=parse(Int, size_match.captures[1]))
end

# SHA-256 of a file as lowercase hex (matches Git LFS `oid sha256:` format).
function _sha256_of_file(path::AbstractString)
    return open(path) do io
        bytes2hex(SHA.sha256(io))
    end
end

# Ordered list of GitHub refs (commit / tag / branch) to try when fetching LFS
# content via `media.githubusercontent.com`. The eventual SHA-256 check guards
# against any of these serving wrong bytes, so the order is just "most likely
# to match first".
function _candidate_refs()
    refs = String[]
    if haskey(ENV, "TURBULENTTRANSPORT_MODELS_REF")
        ref = strip(ENV["TURBULENTTRANSPORT_MODELS_REF"])
        !isempty(ref) && push!(refs, ref)
    end
    root = _pkg_root()
    if isdir(joinpath(root, ".git"))
        try
            push!(refs, readchomp(`git -C $root rev-parse HEAD`))
        catch
        end
    end
    try
        ver = pkgversion(@__MODULE__)
        ver === nothing || push!(refs, "v$ver")
    catch
    end
    push!(refs, "master")
    push!(refs, "main")
    return unique!(refs)
end

# Per-process URL builder override, used only by the test suite to point at a
# local `file://` fixture instead of `media.githubusercontent.com`. In normal
# use this stays `nothing` and we hit GitHub.
const _LFS_URL_OVERRIDE = Ref{Union{Nothing,Function}}(nothing)

function _lfs_media_url(ref::AbstractString, relpath::AbstractString)
    override = _LFS_URL_OVERRIDE[]
    override === nothing || return override(ref, relpath)
    return "https://media.githubusercontent.com/media/$(_LFS_REPO)/$(ref)/$(relpath)"
end

# Can we create files in `dir`? Probes with an actual write so it catches
# read-only filesystems (EROFS, e.g. a container SIF) as well as plain
# permission denials (EACCES, e.g. a shared site install).
function _dir_writable(dir::AbstractString)
    isdir(dir) || return false
    probe = joinpath(dir, ".tt_write_probe_$(getpid())_$(rand(UInt32))")
    try
        touch(probe)
        rm(probe; force=true)
        return true
    catch
        return false
    end
end

# Per-package-version scratch space used to materialize models when the
# install itself is not writable. Keying by version invalidates the cache on
# upgrade; the oid check on every hit guards content within a version.
function _model_cache_dir()
    ver = try
        string(pkgversion(@__MODULE__))
    catch
        "dev"
    end
    return Scratch.get_scratch!(@__MODULE__, "models-v$ver")
end

# Cache destination for a model that cannot be materialized in place. Files
# under the package root mirror their relative layout (so e.g. an .onnx and
# its .data sidecar stay co-located); anything else gets a flat oid-keyed name.
function _cache_dest(path::AbstractString, oid::AbstractString)
    cache = _model_cache_dir()
    root = _pkg_root()
    p = normpath(abspath(path))
    if startswith(p, root)
        return joinpath(cache, relpath(p, root))
    else
        return joinpath(cache, "external", oid[1:16] * "-" * basename(p))
    end
end

# ONNX models may reference an external-data sidecar (<model>.onnx.data) that
# the runtime resolves relative to the .onnx file it opens. Whenever we hand
# back a (possibly cache-redirected) .onnx path, make sure its sidecar exists
# next to it too.
function _ensure_onnx_sidecar!(path::AbstractString, dest::AbstractString)
    endswith(lowercase(path), ".onnx") || return nothing
    sidecar = path * ".data"
    isfile(sidecar) || return nothing
    resolved = ensure_model_file!(sidecar)
    want = dest * ".data"
    if resolved != want && !isfile(want)
        mkpath(dirname(want))
        cp(resolved, want)
    end
    return nothing
end

"""
    ensure_model_file!(path::AbstractString) -> materialized_path

If `path` is a Git LFS pointer stub (Pkg installs don't run `git lfs pull`),
transparently materialize the real model bytes and return the path callers
should load. When the install is writable this happens in place, so the
returned path equals `path`; when it is read-only (e.g. a container SIF or a
shared site install) the bytes are materialized into a per-version scratch
cache instead and *that* path is returned — callers must always use the
returned path, not the argument.

The Git LFS pointer's `oid sha256:...` is the source of truth: every download
candidate (and every cache hit) is accepted only if its SHA-256 matches the
pointer's oid, so a ref that drifts ahead/behind the installed package cannot
silently swap in the wrong weights.
"""
function ensure_model_file!(path::AbstractString)
    isfile(path) || return path
    is_lfs_pointer(path) || return path

    dest = lock(_model_download_lock) do
        # Re-check under the lock: another task may have materialized it.
        is_lfs_pointer(path) || return path

        spec = _lfs_pointer_info(path)
        spec === nothing && error("Malformed Git LFS pointer at '$path'")

        root = _pkg_root()
        rel = relpath(path, root)
        in_place = _dir_writable(dirname(path))

        # 1) In a writable dev checkout, prefer `git lfs pull` — it uses the
        #    oid internally and bypasses any ref guessing.
        if in_place && isdir(joinpath(root, ".git"))
            try
                run(setenv(`git -C $root lfs pull --include $rel`; stderr=devnull, stdin=devnull))
                if !is_lfs_pointer(path) && _sha256_of_file(path) == spec.oid
                    return path
                end
            catch
            end
        end

        # Destination: in place when writable, else the scratch cache (where a
        # previous materialization may already satisfy the oid).
        target = path
        if !in_place
            target = _cache_dest(path, spec.oid)
            if isfile(target) && filesize(target) == spec.size && _sha256_of_file(target) == spec.oid
                return target
            end
            mkpath(dirname(target))
        end

        # 2) Fall back to `media.githubusercontent.com`, trying refs in order
        #    and accepting only the one whose bytes hash to the pointer's oid.
        attempted = String[]
        tmp = target * ".download"
        for ref in _candidate_refs()
            url = _lfs_media_url(ref, rel)
            push!(attempted, ref)
            try
                Downloads.download(url, tmp)
                if isfile(tmp) && filesize(tmp) == spec.size &&
                   !is_lfs_pointer(tmp) && _sha256_of_file(tmp) == spec.oid
                    mv(tmp, target; force=true)
                    return target
                end
            catch
            end
            rm(tmp; force=true)
        end

        error(
            "Could not materialize Git LFS model '$rel' " *
            "(oid sha256:$(spec.oid), size $(spec.size)). " *
            "Tried refs: $(join(attempted, ", ")). " *
            "Check network access to media.githubusercontent.com, or run " *
            "`git lfs pull` in a TurbulentTransport dev checkout."
        )
    end

    _ensure_onnx_sidecar!(path, dest)
    return dest
end

"""
    register_model_path!(path::String; prepend::Bool=true)

Register an additional model search directory. Called by package extensions.

NOTE: Registration must happen at extension load time (in `__init__`), before
any calls to `loadmodel()` or `loadmodelonce()`, to avoid cache misses.
"""
function register_model_path!(path::String; prepend::Bool=true)
    if prepend
        pushfirst!(_MODEL_SEARCH_PATHS, path)
    else
        push!(_MODEL_SEARCH_PATHS, path)
    end
end

"""
    resolve_model_path(spec::AbstractString; extensions=[".bson", ".onnx"])

Resolve a model name or path to an absolute file path.

Resolution order:
1. If `spec` is an existing file (absolute or relative), return it
2. Search registered provider paths (in order, providers first)
3. Search built-in models directory
4. Error with list of available models
"""
function resolve_model_path(spec::AbstractString; extensions::Vector{String}=[".bson", ".onnx"])
    # 1. If spec is an existing file (handles absolute AND relative paths)
    candidates = [spec; spec .* extensions]
    for candidate in candidates
        if isfile(candidate)
            return ensure_model_file!(candidate)
        end
    end

    # 2. Search provider paths (in registration order, newest first)
    for search_dir in _MODEL_SEARCH_PATHS
        for candidate in [spec; spec .* extensions]
            path = joinpath(search_dir, candidate)
            if isfile(path)
                return ensure_model_file!(path)
            end
        end
    end

    # 3. Built-in models (fallback)
    builtin_dir = joinpath(_pkg_root(), "models")
    for candidate in [spec; spec .* extensions]
        path = joinpath(builtin_dir, candidate)
        if isfile(path)
            return ensure_model_file!(path)
        end
    end

    error("Model '$spec' not found. Available: $(join(available_models(), ", "))")
end

"""
    available_models()

Return list of all available model names from registered providers and built-in.
Provider models appear first (and take precedence on name collision).
"""
function available_models()
    models = String[]
    all_paths = vcat(_MODEL_SEARCH_PATHS, [joinpath(dirname(@__DIR__), "models")])
    for dir in all_paths
        isdir(dir) || continue
        for f in readdir(dir)
            if endswith(f, ".bson") || endswith(f, ".onnx")
                push!(models, replace(f, r"\.(bson|onnx)$" => ""))
            end
        end
    end
    unique!(models)  # First-seen wins (providers before builtin)
end

"""
    available_qlnn_bundles()

Return list of available QLNN bundle directory names (one per subdirectory
that contains an `energy_regressor.bson`). Bundles from registered providers
appear first.
"""
function available_qlnn_bundles()
    bundles = String[]
    all_paths = vcat(_MODEL_SEARCH_PATHS, [joinpath(dirname(@__DIR__), "models")])
    for dir in all_paths
        isdir(dir) || continue
        for entry in readdir(dir; join=false)
            sub = joinpath(dir, entry)
            isdir(sub) || continue
            isfile(joinpath(sub, "energy_regressor.bson")) || continue
            push!(bundles, entry)
        end
    end
    unique!(bundles)
end
