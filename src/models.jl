#= ======================== =#
#  Model path resolution API
#= ======================== =#

const _MODEL_SEARCH_PATHS = String[]  # Additional search paths from providers

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
        isfile(candidate) && return candidate
    end

    # 2. Search provider paths (in registration order, newest first)
    for search_dir in _MODEL_SEARCH_PATHS
        for candidate in [spec; spec .* extensions]
            path = joinpath(search_dir, candidate)
            isfile(path) && return path
        end
    end

    # 3. Built-in models (fallback)
    builtin_dir = joinpath(dirname(@__DIR__), "models")
    for candidate in [spec; spec .* extensions]
        path = joinpath(builtin_dir, candidate)
        isfile(path) && return path
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
