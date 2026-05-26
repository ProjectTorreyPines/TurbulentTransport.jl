#= ======================== =#
#  Model path resolution API
#= ======================== =#

using Downloads

const _MODEL_SEARCH_PATHS = String[]  # Additional search paths from providers
const _LFS_POINTER_PREFIX = b"version https://git-lfs.github.com/spec/v1"
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
    prefix_len = min(sz, length(_LFS_POINTER_PREFIX))
    bytes = read(path, prefix_len)
    return startswith(bytes, _LFS_POINTER_PREFIX)
end

function _models_git_ref()
    if haskey(ENV, "TURBULENTTRANSPORT_MODELS_REF")
        ref = strip(ENV["TURBULENTTRANSPORT_MODELS_REF"])
        !isempty(ref) && return ref
    end
    ref_file = joinpath(_pkg_root(), "models", "MODELS_REF")
    if isfile(ref_file)
        ref = strip(read(ref_file, String))
        !isempty(ref) && return ref
    end
    if isdir(joinpath(_pkg_root(), ".git"))
        try
            return readchomp(`git -C $( _pkg_root()) rev-parse HEAD`)
        catch
        end
    end
    try
        return "v$(pkgversion(@__MODULE__))"
    catch
    end
    return "master"
end

function _lfs_media_url(relpath::AbstractString)
    ref = _models_git_ref()
    return "https://media.githubusercontent.com/media/$(_LFS_REPO)/$(ref)/$(relpath)"
end

"""
    ensure_model_file!(path::AbstractString) -> path

Materialize a model file when the on-disk copy is a Git LFS pointer. Pkg checkouts
do not run `git lfs pull`, so CI and fresh installs need this fallback.
"""
function ensure_model_file!(path::AbstractString)
    isfile(path) || return path
    is_lfs_pointer(path) || return path

    lock(_model_download_lock) do
        is_lfs_pointer(path) || return path

        root = _pkg_root()
        relpath = relpath(path, root)

        if isdir(joinpath(root, ".git"))
            try
                run(setenv(`git -C $root lfs pull --include $relpath`; stderr=devnull))
                is_lfs_pointer(path) || return path
            catch
            end
        end

        url = _lfs_media_url(relpath)
        tmp = path * ".download"
        try
            Downloads.download(url, tmp)
            if !isfile(tmp) || filesize(tmp) < 500 || is_lfs_pointer(tmp)
                rm(tmp; force=true)
                error(
                    "Model file '$(relpath)' is a Git LFS pointer and could not be downloaded " *
                    "from $(url). Set TURBULENTTRANSPORT_MODELS_REF to a valid commit/tag, or run " *
                    "`git lfs pull` in a full TurbulentTransport checkout."
                )
            end
            mv(tmp, path; force=true)
        catch err
            rm(tmp; force=true)
            if err isa ErrorException && occursin("Git LFS pointer", err.msg)
                rethrow()
            end
            error(
                "Model file '$(relpath)' is a Git LFS pointer and download failed from $(url): $(sprint(showerror, err))"
            )
        end
    end

    return path
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
