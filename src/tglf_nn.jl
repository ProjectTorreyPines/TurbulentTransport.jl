import Flux
import Flux: NNlib
import Dates
import Memoize
import StatsBase
import Measurements
import BSON

const log_suffix = "_log10"
const n_log_suffix = ncodeunits(log_suffix)

#= ====================================== =#
#  structs/constructors for the TGLFmodel
#= ====================================== =#
# TGLFmodel abstract type, since we could have different models
abstract type TGLFmodel end

# TGLFNNmodel
struct TGLFNNmodel <: TGLFmodel
    fluxmodel::Flux.Chain
    name::String
    date::Dates.DateTime
    xnames::Vector{String}
    ynames::Vector{String}
    xm::Vector{Float64}
    xσ::Vector{Float64}
    ym::Vector{Float64}
    yσ::Vector{Float64}
    xbounds::Array{Float64}
    ybounds::Array{Float64}
    nions::Int

    # Pre-initialized PooledChain for zero-allocation inference (not serialized)
    _pooled_chain::PooledChain
end

function Base.show(io::IO, mime::MIME"text/plain", model::TGLFNNmodel)
    println(io, "TGLFNNmodel")
    println(io, "name: $(length(model.name))")
    println(io, "date: $(model.date)")
    println(io, "nions: $(model.nions)")
    println(io, "xnames ($(length(model.xnames))): $(model.xnames)")
    return println(io, "ynames ($(length(model.ynames))): $(model.ynames)")
end

# TGLFNNensemble
struct TGLFNNensemble <: TGLFmodel
    models::Vector{TGLFNNmodel}
end

function Base.show(io::IO, mime::MIME"text/plain", ensemble::TGLFNNensemble)
    println(io, "TGLFNNensemble")
    println(io, "n models: $(length(ensemble.models))")
    return show(io, mime, ensemble.models[1])
end

function TGLFNNensemble(models::Vector{<:Any})
    return TGLFNNensemble(TGLFNNmodel[model for model in models])
end

function Base.getproperty(ensemble::TGLFNNensemble, field::Symbol)
    if field == :models
        return getfield(ensemble, field)
    elseif field == :fluxmodel
        error("Running TGLF ensemble like a model")
    else
        return getfield(ensemble.models[1], field)
    end
end

#= ====================================== =#
#  Zero-allocation field extraction cache
#= ====================================== =#

# Cache for field symbols Val types. Uses content hash (not objectid) so:
# - Same xnames content → same cache entry (memory efficient)
# - Cache size bounded by unique xnames patterns, not model instances
const _XNAMES_FIELD_SYMBOLS_CACHE = Dict{UInt64, Any}()

"""
    _get_xnames_without_log10_suffix(model::TGLFNNmodel)

Get cached `Val{Tuple{Symbol...}}` of InputTGLF field names from model's xnames.
Strips `_log10` suffix: `"BETAE_log10"` → `:BETAE`
"""
function _get_xnames_without_log10_suffix(model::TGLFNNmodel)
    key = hash(model.xnames)  # content-based, not objectid
    get!(_XNAMES_FIELD_SYMBOLS_CACHE, key) do
        symbols = Tuple(Symbol(endswith(x, log_suffix) ? x[1:end-n_log_suffix] : x) for x in model.xnames)
        Val(symbols)
    end
end

function _get_xnames_without_log10_suffix(ensemble::TGLFNNensemble)
    _get_xnames_without_log10_suffix(ensemble.models[1])
end

"""
    _extract_fields!(inputs::AbstractVector, obj, ::Val{symbols}, index::Int=0)

Zero-allocation field extraction using compile-time unrolled field access.
`symbols` is a tuple of field names known at compile time via `Val`.
`index` is the radial location index for error messages (default 0).
"""
@generated function _extract_fields!(inputs::AbstractVector, obj, ::Val{symbols}, index::Int=0) where {symbols}
    exprs = []
    for (i, s) in enumerate(symbols)
        push!(exprs, quote
            let value = getfield(obj, $(QuoteNode(s)))
                if ismissing(value)
                    _throw_missing_field_error($(QuoteNode(s)), index)
                end
                @inbounds inputs[$i] = value
            end
        end)
    end
    return Expr(:block, exprs..., :inputs)
end

# Function barrier: specializes on concrete Val type, avoiding dynamic dispatch in loop
function _extract_all_inputs!(inputs::AbstractMatrix, input_tglfs::Vector{InputTGLF{T}}, xnames_val::Val) where {T<:Real}
    for (i, input_tglf) in enumerate(input_tglfs)
        _extract_fields!(@view(inputs[:, i]), input_tglf, xnames_val, i)
    end
    return inputs
end

# Error function separated for hot path optimization (@noinline keeps it out of inlined code)
@noinline function _throw_missing_field_error(field::Symbol, index::Int)
    field_str = string(field)
    hint = ""
    if occursin("_5", field_str) || occursin("_6", field_str)
        hint = "\n\nHint: Missing species data (species 5 or 6). If using a TGLFNN model (e.g. 'stfpp' models), try setting:\n  act.ActorTGLF.lump_ions = false\nto ensure ion species are treated separately rather than lumped together."
    end
    error("TGLFNN input field '$field_str' is Missing at radial location $index. Check that all required equilibrium and profile data are properly initialized.$hint")
end

#= ====================================== =#
#  Pooled layer convenience methods
#= ====================================== =#

"""
    poolify(model::TGLFNNmodel)

Convert a TGLFNNmodel's fluxmodel to use pooled layers.
"""
poolify(model::TGLFNNmodel) = poolify(model.fluxmodel)

"""
    PooledChain(model::TGLFNNmodel)

Convenience constructor: creates a PooledChain from a TGLFNNmodel.
"""
PooledChain(model::TGLFNNmodel) = PooledChain(poolify(model.fluxmodel))

#= ============== =#
#  saving/loading  #
#= ============== =#
function mod2dict(model::TGLFNNmodel)
    savedict = Dict()
    for name in fieldnames(TGLFNNmodel)
        name === :_pooled_chain && continue  # Skip cache field (not serialized)
        value = getproperty(model, name)
        savedict[name] = value
    end
    return savedict
end

function mod2dict(ensemble::TGLFNNensemble)
    savedict = Dict()
    for (km, model) in enumerate(ensemble.models)
        savedict[km] = mod2dict(model)
    end
    return savedict
end

function savemodel(model::TGLFmodel, filename::AbstractString)
    if !endswith(filename, ".bson")
        filename = "$(filename).bson"
    end
    if startswith(filename, "/")
        fullpath = filename
    else
        fullpath = dirname(@__DIR__) * "/models/NN_ensembles/" * filename
    end
    BSON.bson(fullpath, mod2dict(model))
    return fullpath
end

Memoize.@memoize function loadmodelonce(filename::String)
    return loadmodel(filename)
end

function dict2mod(savedict::AbstractDict)
    args = []
    for name in fieldnames(TGLFNNmodel)
        if name == :fluxmodel
            savedict[name] = Flux.fmap(Flux.f64, savedict[name])
            push!(args, savedict[name])
        elseif name == :nions
            nions = maximum(map(m -> parse(Int, m[1]), filter(!isnothing, match.(r"_([0-9]+$)", savedict[:xnames])))) - 1
            push!(args, nions)
        elseif name === :_pooled_chain
            push!(args, PooledChain(poolify(savedict[:fluxmodel])))
        else
            push!(args, savedict[name])
        end
    end
    return TGLFNNmodel(args...)
end

function dict2ens(dict::Dict)
    return TGLFNNensemble([dict2mod(modict) for modict in values(dict)])
end

function loadmodel(filename::AbstractString)
    fullpath = resolve_model_path(filename; extensions=[".bson"])
    savedict = BSON.load(fullpath, @__MODULE__)
    if typeof(first(keys(savedict))) <: Integer
        return dict2ens(savedict)
    else
        return dict2mod(savedict)
    end
end

#= ==================================== =#
#  functions to get the fluxes solution
#= ==================================== =#

"""
    flux_array(fluxmodel::TGLFNNmodel, x::AbstractMatrix{T}; ...) where {T<:Real}

Batched inference: processes entire `[N_features, M_samples]` matrix in single forward pass.
"""
function flux_array(fluxmodel::TGLFNNmodel, x::AbstractMatrix{T}; warn_nn_train_bounds::Bool=true, fidelity::Symbol=:TGLFNN) where {T<:Real}
    nouts = length(fluxmodel.ynames)
    if fidelity == :GKNN
        nouts = div(nouts, 2)
    end
    yy = Matrix{T}(undef, nouts, size(x, 2))

    flux_array!(yy, fluxmodel, x; warn_nn_train_bounds, fidelity)
    return yy
end

"""
    flux_array(fluxmodel::TGLFNNmodel, x::AbstractVector{T}; ...) where {T<:Real}

Single-sample inference: processes one vector through the model.
"""
function flux_array(fluxmodel::TGLFNNmodel, x::AbstractVector{T}; warn_nn_train_bounds::Bool=true, fidelity::Symbol=:TGLFNN) where {T<:Real}
    nouts = length(fluxmodel.ynames)
    if fidelity == :GKNN
        nouts = div(nouts, 2)
    end
    yy = Vector{T}(undef, nouts)
    flux_array!(yy, fluxmodel, x; warn_nn_train_bounds, fidelity)
    return yy
end


"""
    flux_array!(out_y::AbstractMatrix{T}, fluxmodel::TGLFNNmodel, x::AbstractMatrix{T};
                warn_nn_train_bounds=true, fidelity=:TGLFNN) where {T<:Real}

In-place batched inference with **zero allocation** (requires AdaptiveArrayPools.jl v0.2.1+).

Processes `[N_features, M_samples]` matrix through the model and writes results to `out_y`.

# Arguments
- `out_y::AbstractMatrix{T}`: Pre-allocated output matrix `[N_outputs, M_samples]`
- `fluxmodel::TGLFNNmodel`: Neural network model
- `x::AbstractMatrix{T}`: Input features `[N_features, M_samples]`
- `warn_nn_train_bounds::Bool=true`: Warn if extrapolating beyond training bounds
- `fidelity::Symbol=:TGLFNN`: Output mode (`:TGLFNN` for denormalized, `:GKNN` for normalized)

"""
@with_pool pool function flux_array!(out_y::AbstractMatrix{T}, fluxmodel::TGLFNNmodel, x::AbstractMatrix{T}; warn_nn_train_bounds::Bool=true, fidelity::Symbol=:TGLFNN) where {T<:Real}
    N, M = size(x)  # N = input features, M = samples

    # acquire! returns Array (not ReshapedArray) to avoid boxing with non-concrete _pooled_chain
    xx = acquire!(pool, T, size(x))

    # Apply log10 transform where needed (determined by feature name)
    @inbounds for i in 1:N
        if contains(fluxmodel.xnames[i], log_suffix)
            for j in 1:M
                xx[i, j] = log10(x[i, j])
            end
        else
            for j in 1:M
                xx[i, j] = x[i, j]
            end
        end
    end

    # Validate bounds (check first sample only to avoid warning spam)
    if warn_nn_train_bounds
        for ix in 1:N
            val = xx[ix, 1]
            if isnan(val) || isinf(val)
                error("$(fluxmodel.xnames[ix]) = $(x[ix, 1]) is not allowed")
            elseif val < fluxmodel.xbounds[ix, 1]
                @warn("Extrapolation $(fluxmodel.xnames[ix])=$(val) is below training bound of $(fluxmodel.xbounds[ix, 1])")
            elseif val > fluxmodel.xbounds[ix, 2]
                @warn("Extrapolation $(fluxmodel.xnames[ix])=$(val) is above training bound of $(fluxmodel.xbounds[ix, 2])")
            end
        end
    end

    # Normalize inputs: (xx - mean) / std
    @inbounds for i in 1:N
        xm_i = fluxmodel.xm[i]
        xσ_i = fluxmodel.xσ[i]
        for j in 1:M
            xx[i, j] = (xx[i, j] - xm_i) / xσ_i
        end
    end

    # Forward pass through neural network (zero allocation via pooled layers)
    fluxmodel._pooled_chain(out_y, xx)

    if fidelity == :GKNN
        return out_y
    elseif fidelity == :TGLFNN
        # Denormalize outputs: out_y * yσ + ym
        nouts = size(out_y, 1)
        @inbounds for i in 1:nouts
            ym_i = fluxmodel.ym[i]
            yσ_i = fluxmodel.yσ[i]
            for j in 1:M
                out_y[i, j] = out_y[i, j] * yσ_i + ym_i
            end
        end
        return out_y
    else
        error("Unknown fidelity mode: $fidelity. Expected :GKNN or :TGLFNN")
    end
end


"""
    flux_array!(out_y::AbstractVector{T}, fluxmodel::TGLFNNmodel, x::AbstractVector{T};
                warn_nn_train_bounds=true, fidelity=:TGLFNN) where {T<:Real}

In-place single-sample inference with **zero allocation** (requires AdaptiveArrayPools.jl v0.2.1+).

Processes one input vector through the model and writes results to `out_y`.

# Arguments
- `out_y::AbstractVector{T}`: Pre-allocated output vector `[N_outputs]`
- `fluxmodel::TGLFNNmodel`: Neural network model
- `x::AbstractVector{T}`: Input features `[N_features]`
- `warn_nn_train_bounds::Bool=true`: Warn if extrapolating beyond training bounds
- `fidelity::Symbol=:TGLFNN`: Output mode (`:TGLFNN` for denormalized, `:GKNN` for normalized)

"""
@with_pool pool function flux_array!(out_y::AbstractVector{T}, fluxmodel::TGLFNNmodel, x::AbstractVector{T}; warn_nn_train_bounds::Bool=true, fidelity::Symbol=:TGLFNN) where {T<:Real}
    N = length(x)

    # acquire! returns Array (not ReshapedArray) to avoid boxing with non-concrete _pooled_chain
    xx = acquire!(pool, T, N)

    # Apply log10 transform where needed (determined by feature name)
    for (ix, name) in enumerate(fluxmodel.xnames)
        xx[ix] = contains(name, "_log10") ? log10(x[ix]) : x[ix]
    end

    # Validate bounds
    if warn_nn_train_bounds
        for ix in 1:N
            if isnan(xx[ix]) || isinf(xx[ix])
                error("$(fluxmodel.xnames[ix]) = $(x[ix]) is not allowed")
            elseif xx[ix] < fluxmodel.xbounds[ix, 1]
                @warn("Extrapolation $(fluxmodel.xnames[ix])=$(xx[ix]) is below training bound of $(fluxmodel.xbounds[ix, 1])")
            elseif xx[ix] > fluxmodel.xbounds[ix, 2]
                @warn("Extrapolation $(fluxmodel.xnames[ix])=$(xx[ix]) is above training bound of $(fluxmodel.xbounds[ix, 2])")
            end
        end
    end

    # Normalize inputs: (xx - mean) / std
    @. xx = (xx - fluxmodel.xm) / fluxmodel.xσ

    # Forward pass through neural network (zero allocation via pooled layers)
    fluxmodel._pooled_chain(out_y, xx)

    if fidelity == :GKNN
        return out_y
    elseif fidelity == :TGLFNN
        # Denormalize outputs: out_y * yσ + ym
        @. out_y = out_y * fluxmodel.yσ + fluxmodel.ym
        return out_y
    else
        error("Unknown fidelity mode: $fidelity. Expected :GKNN or :TGLFNN")
    end
end


#= ====================================== =#
#  Ensemble inference helpers (function barriers for type stability)
#= ====================================== =#

"""
    _flux_array_sequential!(all_yy, fluxensemble, x; ...) -> nothing

Sequential ensemble inference. Reuses single buffer across models (zero-alloc after warmup).
"""
@with_pool pool function _flux_array_sequential!(all_yy::AbstractArray{T,3}, fluxensemble::TGLFNNensemble, x::AbstractArray{T}; warn_nn_train_bounds::Bool=true, fidelity::Symbol=:TGLFNN) where {T<:Real}
    nouts, nsamples, nmodels = size(all_yy)
    each_y = acquire!(pool, T, nouts, nsamples)
    for k in 1:nmodels
        flux_array!(each_y, fluxensemble.models[k], x; warn_nn_train_bounds=(warn_nn_train_bounds && k == 1), fidelity)
        all_yy[:, :, k] = each_y
    end
end

# Vector input variant
@with_pool pool function _flux_array_sequential!(all_yy::AbstractMatrix{T}, fluxensemble::TGLFNNensemble, x::AbstractVector{T}; warn_nn_train_bounds::Bool=true, fidelity::Symbol=:TGLFNN) where {T<:Real}
    nouts, nmodels = size(all_yy)
    each_y = acquire!(pool, T, nouts)
    for k in 1:nmodels
        flux_array!(each_y, fluxensemble.models[k], x; warn_nn_train_bounds=(warn_nn_train_bounds && k == 1), fidelity)
        all_yy[:, k] = each_y
    end
end

"""
    _flux_array_threaded!(all_yy, fluxensemble, x; ...) -> nothing

Threaded ensemble inference. Each thread uses its own pool buffer (thread-safe).
"""
function _flux_array_threaded!(all_yy::AbstractArray{T,3}, fluxensemble::TGLFNNensemble, x::AbstractArray{T}; warn_nn_train_bounds::Bool=true, fidelity::Symbol=:TGLFNN) where {T<:Real}
    nouts, nsamples, nmodels = size(all_yy)
    Threads.@threads for k in 1:nmodels
        @with_pool thread_pool begin
            each_y = acquire!(thread_pool, T, nouts, nsamples)
            flux_array!(each_y, fluxensemble.models[k], x; warn_nn_train_bounds=(warn_nn_train_bounds && k == 1), fidelity)
            all_yy[:, :, k] = each_y
            nothing
        end
    end
end

# Vector input variant
function _flux_array_threaded!(all_yy::AbstractMatrix{T}, fluxensemble::TGLFNNensemble, x::AbstractVector{T}; warn_nn_train_bounds::Bool=true, fidelity::Symbol=:TGLFNN) where {T<:Real}
    nouts, nmodels = size(all_yy)
    Threads.@threads for k in 1:nmodels
        @with_pool thread_pool begin
            each_y = acquire!(thread_pool, T, nouts)
            flux_array!(each_y, fluxensemble.models[k], x; warn_nn_train_bounds=(warn_nn_train_bounds && k == 1), fidelity)
            all_yy[:, k] = each_y
            nothing
        end
    end
end

"""
    flux_array(fluxensemble::TGLFNNensemble, x::AbstractArray{T}; ...) where {T<:Real}

Ensemble batched inference: runs all models in parallel, returns mean (± std if `uncertain=true`).
"""
@with_pool pool function flux_array(fluxensemble::TGLFNNensemble, x::AbstractArray{T}; uncertain::Bool=false, warn_nn_train_bounds::Bool=true, fidelity::Symbol=:TGLFNN) where {T<:Real}
    nmodels = length(fluxensemble.models)
    nouts = length(fluxensemble.models[1].ynames)
    if fidelity == :GKNN
        nouts = div(nouts, 2)
    end
    nsamples = size(x, 2)

    # Store each model's output: (nouts, nsamples, nmodels) for efficient slice access
    all_yy = acquire!(pool, T, nouts, nsamples, nmodels)
    if Threads.nthreads() == 1
        _flux_array_sequential!(all_yy, fluxensemble, x; warn_nn_train_bounds, fidelity)
    else
        _flux_array_threaded!(all_yy, fluxensemble, x; warn_nn_train_bounds, fidelity)
    end

    # Compute mean using broadcasting
    mean_out = zeros(T, nouts, nsamples)
    @inbounds for k in 1:nmodels
        @. @views mean_out += all_yy[:, :, k]
    end
    mean_out ./= nmodels

    if uncertain && nmodels > 1
        if T <: Measurements.Measurement
            return mean_out
        else
            # Compute std using broadcasting
            std_out = zeros(T, nouts, nsamples)
            @inbounds for k in 1:nmodels
                @. @views std_out += (all_yy[:, :, k] - mean_out)^2
            end
            @. std_out = sqrt(std_out / (nmodels - 1))
            return Measurements.measurement.(mean_out, std_out)
        end
    else
        return mean_out
    end
end

"""
    flux_array(fluxensemble::TGLFNNensemble, x::AbstractVector{T}; ...) where {T<:Real}

Ensemble single-sample inference: runs all models on one vector, returns mean (± std if `uncertain=true`).
"""
@with_pool pool function flux_array(fluxensemble::TGLFNNensemble, x::AbstractVector{T}; uncertain::Bool=false, warn_nn_train_bounds::Bool=true, fidelity::Symbol=:TGLFNN) where {T<:Real}
    nmodels = length(fluxensemble.models)
    nouts = length(fluxensemble.models[1].ynames)
    if fidelity == :GKNN
        nouts = div(nouts, 2)
    end

    # Store each model's output: (nouts, nmodels) for efficient slice access
    all_yy = acquire!(pool, T, nouts, nmodels)
    if Threads.nthreads() == 1
        _flux_array_sequential!(all_yy, fluxensemble, x; warn_nn_train_bounds, fidelity)
    else
        _flux_array_threaded!(all_yy, fluxensemble, x; warn_nn_train_bounds, fidelity)
    end

    # Compute mean using broadcasting
    mean_out = zeros(T, nouts)
    @inbounds for k in 1:nmodels
        @. @views mean_out += all_yy[:, k]
    end
    mean_out ./= nmodels

    if uncertain && nmodels > 1
        # Compute std using broadcasting
        std_out = zeros(T, nouts)
        @inbounds for k in 1:nmodels
            @. @views std_out += (all_yy[:, k] - mean_out)^2
        end
        @. std_out = sqrt(std_out / (nmodels - 1))

        if T <: Measurements.Measurement
            return mean_out
        else
            return Measurements.measurement.(mean_out, std_out)
        end
    else
        return mean_out
    end
end

"""
    flux_array(fluxmodel::TGLFmodel, args...; ...)

Vararg convenience: reshapes scalar arguments into matrix and delegates to batched method.
"""
function flux_array(fluxmodel::TGLFmodel, args...; uncertain::Bool=false, warn_nn_train_bounds::Bool=true, fidelity::Symbol=:TGLFNN)
    args = reshape([k for k in args], (length(args), 1))
    return flux_array(fluxmodel, args; uncertain, warn_nn_train_bounds, fidelity)
end

function flux_solution(fluxmodel::TGLFmodel, args...; uncertain::Bool=false, warn_nn_train_bounds::Bool=true, fidelity::Symbol=:TGLFNN)
    return flux_solution(flux_array(fluxmodel, collect(args); uncertain, warn_nn_train_bounds, fidelity)...)
end

# functors for TGLFNNmodel
#= ======================= =#
function (fluxmodel::TGLFmodel)(x::AbstractArray; uncertain::Bool=false, warn_nn_train_bounds::Bool=true, fidelity::Symbol=:TGLFNN)
    return flux_array(fluxmodel, x; uncertain, warn_nn_train_bounds, fidelity)
end

function (fluxmodel::TGLFmodel)(args...; uncertain::Bool=false, warn_nn_train_bounds::Bool=true, fidelity::Symbol=:TGLFNN)
    return flux_solution(fluxmodel, args...; uncertain, warn_nn_train_bounds, fidelity)
end

#= ========== =#
#  run_tglfnn
#= ========== =#
"""
    run_tglfnn(input_tglf::InputTGLF; model_filename, warn_nn_train_bounds, uncertain=false, fidelity=:TGLFNN) -> GACODE.FluxSolution

Run TGLFNN starting from a InputTGLF, using a specific `model_filename`.

If the model is an ensemble of NNs, then the output can be uncertain (using the Measurements.jl package).

The warn_nn_train_bounds checks against the standard deviation of the inputs to warn if evaluation is likely outside of training bounds.

Returns a `flux_solution` structure

NOTE: Single-input convenience wrapper. Delegates to the vector version.

See [`run_tglfnn(::Vector{InputTGLF})`](@ref) for details.
"""
function run_tglfnn(input_tglf::InputTGLF{T}; model_filename::String, uncertain::Bool=false, warn_nn_train_bounds::Bool, fidelity::Symbol=:TGLFNN) where {T<:Real}
    return run_tglfnn([input_tglf]; model_filename, uncertain, warn_nn_train_bounds, fidelity)[1]
end

"""
    run_tglfnn(input_tglfs::Vector{InputTGLF{T}}; model_filename::String, uncertain::Bool=false, warn_nn_train_bounds::Bool, fidelity::Symbol=:TGLFNN) where {T<:Real}

Run TGLFNN for multiple InputTGLF, using a specific `model_filename`.

This is more efficient than running TGLFNN on each individual InputTGLFs.

If the model is an ensemble of NNs, then the output can be uncertain (using the Measurements.jl package).

The warn_nn_train_bounds checks against the standard deviation of the inputs to warn if evaluation is likely outside of training bounds.

Returns a vector of `flux_solution` structures
"""
@with_pool pool function run_tglfnn(input_tglfs::Vector{InputTGLF{T}}; model_filename::String, uncertain::Bool=false, warn_nn_train_bounds::Bool, fidelity::Symbol=:TGLFNN) where {T<:Real}
    if occursin("stfpp", model_filename) || occursin("tefpp", model_filename)
        for it in input_tglfs
            _apply_stfpp_transform!(it; dtf=0.5, device="")
        end
    end
    if model_filename == "sat3_em_d3d_azf-1" && fidelity == :GKNN
        tglfmod = loadmodelonce(model_filename * "_tglfnn24")
    else
        tglfmod = loadmodelonce(model_filename)
    end

    inputs = acquire_view!(pool, T, length(tglfmod.xnames), length(input_tglfs))

    # Extract input fields using @generated function for zero-allocation
    xnames_val = _get_xnames_without_log10_suffix(tglfmod)
    _extract_all_inputs!(inputs, input_tglfs, xnames_val)

    tmp = flux_array(tglfmod, inputs; uncertain, warn_nn_train_bounds, fidelity=:TGLFNN)

    # Handle models with radial-dependent variants
    k_rminloc = nothing
    if model_filename in ("sat0quench_em_d3d_azf+1_withnegD", "sat1_em_d3d_azf-1_withnegD", "sat2_em_d3d_azf-1_withnegD", "sat3_em_d3d_azf-1_withnegD")
        for (k, item) in enumerate(tglfmod.xnames)
            if item == "RMIN_LOC"
                k_rminloc = k
                break
            end
        end

        # For sat3_em_d3d_azf-1_withnegD + GKNN, blending is done in the GKNN block below
        if model_filename != "sat3_em_d3d_azf-1_withnegD" || fidelity != :GKNN
            if k_rminloc === nothing
                @warn "RMIN_LOC not found in xnames for radial-dependent model blending"
            else
                tglfmod2 = loadmodelonce(replace(model_filename, "d3d" => "d3dnearedge"))
                tglfmod3 = loadmodelonce(replace(model_filename, "d3d" => "d3dedge"))
                tmp2 = flux_array(tglfmod2, inputs; uncertain, warn_nn_train_bounds, fidelity=:TGLFNN)
                tmp3 = flux_array(tglfmod3, inputs; uncertain, warn_nn_train_bounds, fidelity=:TGLFNN)
                for i in eachindex(input_tglfs)
                    if inputs[k_rminloc, i] >= 0.881 && inputs[k_rminloc, i] < 0.975
                        tmp[:, i] .= tmp2[:, i]
                    elseif inputs[k_rminloc, i] >= 0.975
                        tmp[:, i] .= tmp3[:, i]
                    end
                end
            end
        end
    end
    if fidelity == :GKNN
        supported_gknn_models = ("sat3_em_d3d_azf-1", "sat3_em_d3d+mastu+nstx_azf-1", "sat3_em_d3d_azf-1_withnegD", "sat3_em_d3d_azf-1_gkdb", "sat2_em_d3d+mastu+nstx_azf-1", "sat3_em_d3d+mastu_azf-1")
        if !(model_filename in supported_gknn_models)
            error("GKNN fidelity is not supported for model '$model_filename'. Supported models are: $(join(supported_gknn_models, ", "))")
        end

        if model_filename == "sat3_em_d3d_azf-1"
            gk_inputs = acquire_view!(pool, T, size(inputs, 1) + 1, size(inputs, 2))
            gk_inputs[1:end-1, :] = inputs

            for (i, postfix) in enumerate(("_gknng24", "_gknnp24", "_gknne24", "_gknni24"))
                gk_inputs[end, :] = tmp[i, :]
                gknn_model = loadmodelonce(model_filename * postfix)
                err = flux_array(gknn_model, gk_inputs; uncertain, warn_nn_train_bounds, fidelity)[:]
                tmp[i, :] .*= err
            end
        elseif model_filename == "sat3_em_d3d_azf-1_withnegD"
            gk_inputs = acquire_view!(pool, T, size(inputs, 1) + 4, size(inputs, 2))
            gk_inputs[1:end-4, :] = inputs
            if k_rminloc === nothing
                @warn "RMIN_LOC not found in xnames for GKNN edge blending"
                gk_inputs[end-3:end, :] = tmp
                gknn31 = loadmodelonce(model_filename * "_gknn31")
                err = flux_array(gknn31, gk_inputs; uncertain, warn_nn_train_bounds, fidelity)
                tmp .*= err
            else
                # Load nearedge and edge base models
                tglfmod2 = loadmodelonce(replace(model_filename, "d3d" => "d3dnearedge"))
                tglfmod3 = loadmodelonce(replace(model_filename, "d3d" => "d3dedge"))
                tmp2 = flux_array(tglfmod2, inputs; uncertain, warn_nn_train_bounds, fidelity=:TGLFNN)
                tmp3 = flux_array(tglfmod3, inputs; uncertain, warn_nn_train_bounds, fidelity=:TGLFNN)

                # Core region (RMIN_LOC < 0.881): _gknn31 applied to d3d flux
                gknn31 = loadmodelonce(model_filename * "_gknn31")
                gk_inputs[end-3:end, :] = tmp
                err1 = flux_array(gknn31, gk_inputs; uncertain, warn_nn_train_bounds, fidelity)

                # Near-edge and edge regions: _gknn37 applied to nearedge/edge flux
                gknn37 = loadmodelonce(model_filename * "_gknn37")
                gk_inputs[end-3:end, :] = tmp2
                err2 = flux_array(gknn37, gk_inputs; uncertain, warn_nn_train_bounds, fidelity)
                gk_inputs[end-3:end, :] = tmp3
                err3 = flux_array(gknn37, gk_inputs; uncertain, warn_nn_train_bounds, fidelity)

                for i in eachindex(input_tglfs)
                    if inputs[k_rminloc, i] >= 0.881 && inputs[k_rminloc, i] < 0.975
                        tmp[:, i] .= tmp2[:, i] .* err2[:, i]
                    elseif inputs[k_rminloc, i] >= 0.975
                        tmp[:, i] .= tmp3[:, i] .* err3[:, i]
                    else
                        tmp[:, i] .*= err1[:, i]
                    end
                end
            end
        elseif model_filename in ("sat3_em_d3d+mastu+nstx_azf-1", "sat3_em_d3d_azf-1_gkdb", "sat2_em_d3d+mastu+nstx_azf-1")
            gk_inputs = acquire_view!(pool, T, size(inputs, 1) + 4, size(inputs, 2))
            gk_inputs[1:end-4, :] = inputs
            gk_inputs[end-3:end, :] = tmp

            gknn = loadmodelonce(model_filename * "_gknn31")
            err = flux_array(gknn, gk_inputs; uncertain, warn_nn_train_bounds, fidelity)
            tmp .*= err
            if model_filename == "sat3_em_d3d_azf-1_gkdb"
                gk_inputs[end-3:end, :] = tmp
                gkdb = loadmodelonce(model_filename * "_gknn31_cgyro")
                gkdb_err = flux_array(gkdb, gk_inputs; uncertain, warn_nn_train_bounds, fidelity)
                tmp .*= gkdb_err
            end
        elseif model_filename == "sat3_em_d3d+mastu_azf-1"
            gk_inputs = acquire_view!(pool, T, size(inputs, 1) + 4, size(inputs, 2))
            gk_inputs[1:end-4, :] = inputs
            gk_inputs[end-3:end, :] = tmp

            gknn = loadmodelonce(model_filename * "_gknn36")
            err = flux_array(gknn, gk_inputs; uncertain, warn_nn_train_bounds, fidelity)
            tmp .*= err
        end
    end

    sol = [flux_solution(@view(tmp[:, i])) for i in eachindex(input_tglfs)]
    return sol
end

"""
    run_tglfnn(data::Dict; model_filename::String, uncertain::Bool=false, warn_nn_train_bounds::Bool, fidelity::Symbol=:TGLFNN)

Run TGLFNN from a dictionary, using a specific `model_filename`.

If the model is an ensemble of NNs, then the output can be uncertain (using the Measurements.jl package).

The warn_nn_train_bounds checks against the standard deviation of the inputs to warn if evaluation is likely outside of training bounds.

Returns a dictionary with fluxes
"""
function run_tglfnn(data::Dict; model_filename::String, uncertain::Bool=false, warn_nn_train_bounds::Bool, fidelity::Symbol=:TGLFNN)
    if occursin("stfpp", model_filename) || occursin("tefpp", model_filename)
        _apply_stfpp_transform!(data; dtf=0.5, device="")
    end
    if model_filename == "sat3_em_d3d_azf-1" && fidelity == :GKNN
        tglfmod = loadmodelonce(model_filename * "_tglfnn24")
    else
        tglfmod = loadmodelonce(model_filename)
    end
    xnames = [replace(name, "_log10" => "") for name in tglfmod.xnames]
    x = collect(transpose(reduce(hcat, [Float64.(data[name]) for name in xnames])))
    y = tglfmod(x; uncertain, warn_nn_train_bounds, fidelity=:TGLFNN)
    if fidelity == :GKNN
        supported_gknn_models = ("sat3_em_d3d_azf-1", "sat3_em_d3d+mastu+nstx_azf-1", "sat3_em_d3d_azf-1_withnegD", "sat3_em_d3d_azf-1_gkdb", "sat2_em_d3d+mastu+nstx_azf-1", "sat3_em_d3d+mastu_azf-1")
        if !(model_filename in supported_gknn_models)
            error("GKNN fidelity is not supported for model '$model_filename'. Supported models are: $(join(supported_gknn_models, ", "))")
        end
        if model_filename == "sat3_em_d3d_azf-1"
            gknng = loadmodelonce(model_filename * "_gknng24")
            err_g = gknng(vcat(x, y[1])...; uncertain, warn_nn_train_bounds, fidelity)
            y[1] .*= err_g
            gknnp = loadmodelonce(model_filename * "_gknnp24")
            err_p = gknnp(vcat(x, y[2])...; uncertain, warn_nn_train_bounds, fidelity)
            y[2] .*= err_p
            gknne = loadmodelonce(model_filename * "_gknne24")
            err_e = gknne(vcat(x, y[3])...; uncertain, warn_nn_train_bounds, fidelity)
            y[3] .*= err_e
            gknni = loadmodelonce(model_filename * "_gknni24")
            err_i = gknni(vcat(x, y[4])...; uncertain, warn_nn_train_bounds, fidelity)
            y[4] .*= err_i
        elseif model_filename == "sat3_em_d3d_azf-1_withnegD"
            k_rminloc = findfirst(isequal("RMIN_LOC"), xnames)
            if k_rminloc === nothing
                @warn "RMIN_LOC not found in xnames for GKNN edge blending"
                gknn31 = loadmodelonce(model_filename * "_gknn31")
                err = gknn31(vcat(x, y)...; uncertain, warn_nn_train_bounds, fidelity)
                y .*= err
            else
                tglfmod2 = loadmodelonce(replace(model_filename, "d3d" => "d3dnearedge"))
                tglfmod3 = loadmodelonce(replace(model_filename, "d3d" => "d3dedge"))
                y2 = flux_array(tglfmod2, x; uncertain, warn_nn_train_bounds, fidelity=:TGLFNN)
                y3 = flux_array(tglfmod3, x; uncertain, warn_nn_train_bounds, fidelity=:TGLFNN)

                gknn31 = loadmodelonce(model_filename * "_gknn31")
                gknn37 = loadmodelonce(model_filename * "_gknn37")
                err1 = flux_array(gknn31, vcat(x, y); uncertain, warn_nn_train_bounds, fidelity)
                err2 = flux_array(gknn37, vcat(x, y2); uncertain, warn_nn_train_bounds, fidelity)
                err3 = flux_array(gknn37, vcat(x, y3); uncertain, warn_nn_train_bounds, fidelity)

                for i in axes(x, 2)
                    if x[k_rminloc, i] >= 0.881 && x[k_rminloc, i] < 0.975
                        y[:, i] .= y2[:, i] .* err2[:, i]
                    elseif x[k_rminloc, i] >= 0.975
                        y[:, i] .= y3[:, i] .* err3[:, i]
                    else
                        y[:, i] .*= err1[:, i]
                    end
                end
            end
        elseif model_filename in ("sat3_em_d3d+mastu+nstx_azf-1", "sat3_em_d3d_azf-1_gkdb", "sat2_em_d3d+mastu+nstx_azf-1")
            gknn = loadmodelonce(model_filename * "_gknn31")
            err = gknn(vcat(x, y)...; uncertain, warn_nn_train_bounds, fidelity)
            y .*= err
            if model_filename == "sat3_em_d3d_azf-1_gkdb"
                gkdb = loadmodelonce(model_filename * "_gknn31_cgyro")
                gkdb_err = gkdb(vcat(x, y)...; uncertain, warn_nn_train_bounds, fidelity)
                y .*= gkdb_err
            end
        elseif model_filename == "sat3_em_d3d+mastu_azf-1"
            gknn = loadmodelonce(model_filename * "_gknn36")
            err = gknn(vcat(x, y)...; uncertain, warn_nn_train_bounds, fidelity)
            y .*= err
        end
    end
    ynames = [replace(name, "OUT_" => "") for name in tglfmod.ynames]
    return Dict(name => y[k, :] for (k, name) in enumerate(ynames))
end

if !isdefined(@__MODULE__, :_ort_loaded)
    const _ort_loaded = Ref(false)
end
if !isdefined(@__MODULE__, :_sess_cache)
    const _sess_cache = Dict{Tuple{String,Int,Int}, Any}()
end

function _ensure_onnx_env!()
    if !haskey(ENV, "OMP_NUM_THREADS")
        nt = something(tryparse(Int, get(ENV, "SLURM_CPUS_PER_TASK", "")), 1)
        ENV["OMP_NUM_THREADS"] = string(max(nt, 1))
    end
    ENV["OMP_PROC_BIND"] = "false"
    ENV["KMP_AFFINITY"]  = "disabled"
    pop!(ENV, "GOMP_CPU_AFFINITY", nothing)
end

"Import ONNXRunTime only after env is finalized (no const alias!)."
function _load_ort!()
    _ort_loaded[] && return
    @eval import ONNXRunTime
    _ort_loaded[] = true
end

function _resolve_model_path(onnx_path::AbstractString)
    return resolve_model_path(onnx_path; extensions=[".onnx"])
end

function load_onnx_model(onnx_path::String; intra_threads::Int=1, inter_threads::Int=1)
    _ensure_onnx_env!(); _load_ort!()
    so = try
        s = ONNXRunTime.create_session_options()
        try ONNXRunTime.set_intra_op_num_threads!(s, intra_threads) catch end
        try ONNXRunTime.set_inter_op_num_threads!(s, inter_threads) catch end
        try ONNXRunTime.set_graph_optimization_level!(s, :ORT_ENABLE_ALL) catch end
        try ONNXRunTime.set_execution_mode_sequential!(s) catch end
        s
    catch
        nothing
    end
    onnx_path = _resolve_model_path(onnx_path)
    return isnothing(so) ?
        ONNXRunTime.load_inference(ONNXRunTime.testdatapath(onnx_path)) :
        ONNXRunTime.load_inference(ONNXRunTime.testdatapath(onnx_path); session_options=so)
end

function get_onnx_session(onnx_path::String; intra_threads::Int=1, inter_threads::Int=1)
    _ensure_onnx_env!(); _load_ort!()
    onnx_path = _resolve_model_path(onnx_path)
    key = (onnx_path, intra_threads, inter_threads)
    if haskey(_sess_cache, key)
        return _sess_cache[key]
    end
    so = try
        s = ONNXRunTime.create_session_options()
        try ONNXRunTime.set_intra_op_num_threads!(s, intra_threads) catch end
        try ONNXRunTime.set_inter_op_num_threads!(s, inter_threads) catch end
        try ONNXRunTime.set_graph_optimization_level!(s, :ORT_ENABLE_ALL) catch end
        try ONNXRunTime.set_execution_mode_sequential!(s) catch end
        s
    catch
        nothing
    end
    sess = isnothing(so) ?
        ONNXRunTime.load_inference(ONNXRunTime.testdatapath(onnx_path)) :
        ONNXRunTime.load_inference(ONNXRunTime.testdatapath(onnx_path); session_options=so)
    _sess_cache[key] = sess
    return sess
end

# Get (or build) a session once
function _session(onnx_path::String; intra_threads::Int=1, inter_threads::Int=1)
    try
        return get_onnx_session(onnx_path; intra_threads=intra_threads, inter_threads=inter_threads)
    catch
        return load_onnx_model(onnx_path; intra_threads=intra_threads, inter_threads=inter_threads)
    end
end

# Build X as [N, F] Float32 without intermediate allocations; supports InputTGLF{T}
function _build_X(input_tglfs::AbstractVector{TJLF.InputTGLF{T}}, xnames::Vector{String}) where {T}
    N = length(input_tglfs); F = length(xnames)
    X = Matrix{Float32}(undef, N, F)
    @inbounds for i in 1:N
        t = input_tglfs[i]
        for j in 1:F
            name = xnames[j]
            key  = replace(name, "_log10" => "")
            v    = getfield(t, Symbol(key))
            X[i, j] = occursin("_log10", name) ? log10(Float32(v)) : Float32(v)
        end
    end
    return X
end

# Extract output from ORT call (handles NamedTuple vs Dict and [M,N] vs [N,M])
@inline function _extract_Y(res, X_rows::Int, outdim::Int)
    out = hasproperty(res, :output) ? getfield(res, :output) : res["output"]
    Y = out
    # If returned as [M,N], flip to [N,M]
    if size(Y,1) == outdim && size(Y,2) == X_rows
        Y = permutedims(Y)
    end
    return Y
end

"""
    run_tglfnn_onnx(input_tglfs, onnx_path, xnames, ynames; intra_threads=1, inter_threads=1)

Run a TGLF-NN model exported to ONNX through ONNXRuntime, as an alternative to
the BSON/Flux path of [`run_tglfnn`](@ref).

`onnx_path` is resolved against the registered model search paths and the
built-in `models/` directory; `xnames`/`ynames` are the model's input/output
feature names (a trailing `_log10` on an `xname` triggers a `log10` transform of
that feature). `intra_threads`/`inter_threads` control the ONNXRuntime session
thread pools (sessions are cached and reused).

Three input forms are supported and dispatch to matching output shapes:
- a single `InputTGLF` -> a `Vector` of the (reordered) output fluxes,
- a `Vector{InputTGLF}` -> a `Vector{GACODE.FluxSolution}` (one per radius),
- a `Dict` of feature-name => vector -> a `Dict` of output-name => vector.
"""
function run_tglfnn_onnx(input_tglfs::AbstractVector{TJLF.InputTGLF{T}},
                         onnx_path::String,
                         xnames::Vector{String},
                         ynames::Vector{String};
                         intra_threads::Int=1, inter_threads::Int=1) where {T<:Real}

    sess = _session(onnx_path; intra_threads=intra_threads, inter_threads=inter_threads)
    X = _build_X(input_tglfs, xnames)                     # [N,F]
    res = sess((; input = X))                             # NamedTuple or Dict
    Y  = _extract_Y(res, size(X,1), length(ynames))       # [N,M]
    cols = [1, 4, 2, 3]
    Yv = @view Y[:, cols]

    N = size(Yv, 1)
    Tsol = typeof(flux_solution(Yv[1,1], Yv[1,2], Yv[1,3], Yv[1,4]))
    sol = Vector{Tsol}(undef, N)
    @inbounds for i in 1:N
        sol[i] = flux_solution(Yv[i,1], Yv[i,2], Yv[i,3], Yv[i,4])
    end
    return sol
end

function run_tglfnn_onnx(data::Dict,
                         onnx_path::String,
                         xnames::Vector{String},
                         ynames::Vector{String};
                         intra_threads::Int=1, inter_threads::Int=1)::Dict

    sess = _session(onnx_path; intra_threads=intra_threads, inter_threads=inter_threads)

    # Build X :: [N,F] from Dict data
    xclean = replace.(xnames, "_log10" => "")
    N = length(data[xclean[1]]); F = length(xnames)
    X = Matrix{Float32}(undef, N, F)
    @inbounds for j in 1:F
        col = data[xclean[j]]
        @assert length(col) == N
        if occursin("_log10", xnames[j])
            for i in 1:N; X[i,j] = log10(Float32(col[i])); end
        else
            for i in 1:N; X[i,j] = Float32(col[i]); end
        end
    end

    res = sess((; input = X))
    Y   = _extract_Y(res, size(X,1), length(ynames))      # [N,M]
    cols = [1, 4, 2, 3]
    Yv  = @view Y[:, cols]
    ynames_clean = replace.(ynames, "OUT_" => "")

    return Dict(name => @view(Yv[:,k]) for (k, name) in enumerate(ynames_clean))
end

function run_tglfnn_onnx(input_tglf::TJLF.InputTGLF{T},
                         onnx_path::String,
                         xnames::Vector{String},
                         ynames::Vector{String};
                         intra_threads::Int=1, inter_threads::Int=1) where {T<:Real}

    sess = _session(onnx_path; intra_threads=intra_threads, inter_threads=inter_threads)

    # X is 1×F
    F = length(xnames)
    X = Matrix{Float32}(undef, 1, F)
    @inbounds for j in 1:F
        name = xnames[j]
        key  = replace(name, "_log10" => "")
        v    = getfield(input_tglf, Symbol(key))
        X[1,j] = occursin("_log10", name) ? log10(Float32(v)) : Float32(v)
    end

    res = sess((; input = X))
    Y   = _extract_Y(res, 1, length(ynames))              # [1,M]
    cols = [1, 4, 2, 3]
    y    = vec(@view Y[:, cols])                          # length M (reordered)
    return y
end

"""
    flux_solution(xx::Vararg{T}) where {T<:Real}
    flux_solution(xx::AbstractVector{T}) where {T<:Real}

Construct a `FluxSolution` from scalar arguments or a vector.

Accepts either variadic arguments (scalars) or an `AbstractVector`.

    flux_solution(1.0, 2.0, 3.0, 4.0, 5.0, 6.0)
    flux_solution([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])

results in

    Qe = 1.0
    Qi = 2.0
    Γe = 3.0
    Γi = [4.0, 5.0]
    Πi = 6.0

NOTE: for backward compatibility with old TGLF-NN models, if number of arguments is 4 then

    flux_solution(1.0, 2.0, 3.0, 4.0)

results in

    Qe = 3.0
    Qi = 4.0
    Γe = 1.0
    Γi = []
    Πi = 2.0
"""
flux_solution(xx::Vararg{T}) where {T<:Real} = _flux_solution_impl(xx)
flux_solution(xx::AbstractVector{T}) where {T<:Real} = _flux_solution_impl(xx)

# Common implementation for both Vararg (Tuple) and AbstractVector
@inline function _flux_solution_impl(xx)
    T = eltype(xx)
    n_fields = length(xx)
    if n_fields == 4
        ENERGY_FLUX_e = 3
        ENERGY_FLUX_i = 4
        PARTICLE_FLUX_e = 1
        STRESS_TOR_i = 2
        return GACODE.FluxSolution{T}(xx[ENERGY_FLUX_e], xx[ENERGY_FLUX_i], xx[PARTICLE_FLUX_e], T[], xx[STRESS_TOR_i])
    else
        ENERGY_FLUX_e = n_fields - 1
        ENERGY_FLUX_i = n_fields
        PARTICLE_FLUX_e = 1
        PARTICLE_FLUX_i = 2:n_fields-3
        STRESS_TOR_i = n_fields - 2
        return GACODE.FluxSolution{T}(xx[ENERGY_FLUX_e], xx[ENERGY_FLUX_i], xx[PARTICLE_FLUX_e], T[xx[i] for i in PARTICLE_FLUX_i], xx[STRESS_TOR_i])
    end
end


export run_tglfnn, run_tglfnn_onnx, model_selector

# ------------------------------------------------------------
# Model Selector - Core Implementation
# ------------------------------------------------------------

# Helper to unwrap InputTGLFs wrapper object
_unwrap_input_tglfs(result) = hasproperty(result, :tglfs) ? getfield(result, :tglfs) : result

"""
    _model_selector_core(input_tglfs::Vector{InputTGLF{T}}, rho_values::Union{Vector,Nothing}; max_models::Int=3, filter_sat_rule::Union{Symbol,Nothing}=nothing, electromagnetic::Bool=true, ground_truth::Bool=true, show_fluxes::Bool=false, verbose::Bool=true) where T

Core implementation of model selector that works with InputTGLF vectors.

# Arguments
- `filter_sat_rule`: Optional saturation rule filter (e.g., :sat1). If provided, only models matching this sat rule are tested.
- `electromagnetic`: Filter to electromagnetic (true, default) or electrostatic (false) models
- `ground_truth`: Run TJLF and Fortran TGLF for ground truth comparison and accuracy-based ranking (default: true). Set to false to skip and rank by confidence only.
- `show_fluxes`: Print per-flux values and vectorized relative errors (default: false). When false, only the mean relative error is shown.

# Notes
- Automatically skips models with "gknn", "qlnn", or "edge" in their name (these are correction/specialized models, not standalone models)
"""
function _model_selector_core(input_tglfs::Vector{InputTGLF{T}}, rho_values::Union{Vector,Nothing}; max_models::Int=3, filter_sat_rule::Union{Symbol,Nothing}=nothing, electromagnetic::Bool=true, ground_truth::Bool=true, show_fluxes::Bool=false, verbose::Bool=true) where T
    # Make a copy to prevent in-place modifications by stfpp/tefpp models
    input_tglfs = deepcopy(input_tglfs)

    # Snapshot ALPHA_ZF before apply_presets! overwrites it (presets always set ALPHA_ZF=-1)
    original_alpha_zf = [it.ALPHA_ZF for it in input_tglfs]

    # Apply presets for TJLF consistency (same as TGLF Fortran USE_PRESETS=.true.)
    for it in input_tglfs
        apply_presets!(it)
    end

    # Run TJLF (Julia) for ground truth
    tjlf_sols = if !ground_truth
        nothing
    else
        try
            verbose && println("Running TJLF (Julia) for ground truth...")
            sols = [run_tjlf(it) for it in input_tglfs]
            verbose && println("  TJLF done\n")
            sols
        catch e
            verbose && @warn "TJLF unavailable, ranking by confidence only: $e"
            nothing
        end
    end

    # Try Fortran TGLF as additional ground truth (only if in PATH)
    tglf_sols = if !ground_truth
        nothing
    else
        try
            verbose && println("Running Fortran TGLF for ground truth...")
            sols = [run_tglf(it) for it in input_tglfs]
            verbose && println("  Fortran TGLF done\n")
            sols
        catch e
            verbose && println("  Fortran TGLF not available ($(typeof(e)))\n")
            nothing
        end
    end

    # Use Fortran TGLF if available, otherwise fall back to TJLF
    ground_truth_sols = tglf_sols !== nothing ? tglf_sols : tjlf_sols
    ground_truth_label = tglf_sols !== nothing ? "TGLF" : (tjlf_sols !== nothing ? "TJLF" : "none")

    # Print TJLF vs Fortran TGLF comparison to verify apply_presets! consistency
    if verbose && tjlf_sols !== nothing && tglf_sols !== nothing
        flux_names = ["PARTICLE_FLUX_e", "ENERGY_FLUX_e", "ENERGY_FLUX_i", "STRESS_TOR_i"]
        thresholds = [0.03, 0.3, 0.3, 0.3]
        rho_labels = rho_values !== nothing ? rho_values : collect(1:length(input_tglfs))
        all_agree = true
        flux_lines = String[]
        for (rho_idx, rho) in enumerate(rho_labels)
            tj = tjlf_sols[rho_idx]
            tf = tglf_sols[rho_idx]
            tjlf_fluxes = [tj.PARTICLE_FLUX_e, tj.ENERGY_FLUX_e, tj.ENERGY_FLUX_i, tj.STRESS_TOR_i]
            tglf_fluxes = [tf.PARTICLE_FLUX_e, tf.ENERGY_FLUX_e, tf.ENERGY_FLUX_i, tf.STRESS_TOR_i]
            push!(flux_lines, "  RMIN_LOC=$(round(rho, digits=3)):")
            for (name, tjlf_val, tglf_val, thresh) in zip(flux_names, tjlf_fluxes, tglf_fluxes, thresholds)
                rel_err = abs(tjlf_val - tglf_val) / (abs(tglf_val) + thresh)
                agree_1pct = rel_err < 0.01
                agree_abs = abs(tjlf_val - tglf_val) < thresh
                if !agree_1pct || !agree_abs
                    all_agree = false
                end
                push!(flux_lines, "    $name: TJLF=$(round(tjlf_val, digits=4))  TGLF=$(round(tglf_val, digits=4))  rel_err=$(round(rel_err, digits=4))  [1%:$(agree_1pct ? "✓" : "✗")  abs<$thresh:$(agree_abs ? "✓" : "✗")]")
            end
        end
        if all_agree
            println("Ground truth (TJLF Julia vs Fortran TGLF): ✓ all fluxes agree\n")
        else
            println("Ground truth (TJLF Julia vs Fortran TGLF): ✗ discrepancies found")
            foreach(println, flux_lines)
            println()
        end
    end

    # Get all available models and apply filters
    all_models = available_models()

    # Skip correction/specialized models
    all_models = filter(m -> !occursin("gknn", m) && !occursin("qlnn", m) && !occursin("edge", m), all_models)

    # Skip azf+1 models when ALPHA_ZF=-1, but not for sat0/sat0quench which rely on azf+1 models
    sat0_rules = (:sat0, :sat0quench)
    is_sat0 = filter_sat_rule !== nothing && filter_sat_rule in sat0_rules
    if !is_sat0 && all(==(-1.0), original_alpha_zf)
        all_models = filter(m -> !endswith(m, "azf+1"), all_models)
    end

    # Filter by electromagnetic vs electrostatic
    physics_str = electromagnetic ? "_em_" : "_es_"
    all_models = filter(m -> occursin(physics_str, m), all_models)

    # Filter by saturation rule if specified (unless :all)
    if filter_sat_rule !== nothing && filter_sat_rule != :all
        sat_str = string(filter_sat_rule)
        all_models = filter(m -> startswith(m, sat_str), all_models)
    end

    # Print filtering summary
    if verbose
        physics_type = electromagnetic ? "electromagnetic" : "electrostatic"
        println("Physics: $physics_type")
        if filter_sat_rule !== nothing && filter_sat_rule != :all
            println("Saturation rule: $(filter_sat_rule)")
        end
        println("Testing $(length(all_models)) models\n")
    end

    # Run all models and collect results
    all_results = Dict{String, Tuple{Bool, Any}}()

    for model_name in all_models
        try
            flux_sols = run_tglfnn(input_tglfs;
                                  model_filename=model_name,
                                  uncertain=true,
                                  warn_nn_train_bounds=false,
                                  fidelity=:TGLFNN)
            all_results[model_name] = (true, flux_sols)
        catch e
            all_results[model_name] = (false, e)
        end
    end

    successful_models = [name for (name, (success, _)) in all_results if success]

    # Rank models at each rho location
    rho_vec = rho_values === nothing ? collect(1:length(input_tglfs)) : rho_values
    rankings = []

    for (rho_idx, rho) in enumerate(rho_vec)
        tjlf_sol = ground_truth_sols === nothing ? nothing : ground_truth_sols[rho_idx]

        # Compute confidence and relative error to TJLF for each successful model
        model_scores = []

        for model_name in successful_models
            flux_sol = all_results[model_name][2][rho_idx]

            # Confidence: relative uncertainty = unc / (|value| + threshold)
            fluxes = [flux_sol.PARTICLE_FLUX_e, flux_sol.ENERGY_FLUX_e,
                      flux_sol.ENERGY_FLUX_i, flux_sol.STRESS_TOR_i]
            thresholds = [0.03, 0.3, 0.3, 0.3]

            confidence = 0.0
            for (flux, thresh) in zip(fluxes, thresholds)
                val = abs(Measurements.value(flux))
                unc = Measurements.uncertainty(flux)
                confidence += abs(unc / (val + thresh))
            end
            confidence /= length(fluxes)

            # Relative error to ground truth: |pred - gt| / (|gt| + 1e-6)
            rel_errors_vec = zeros(4)
            if tjlf_sol !== nothing
                tjlf_fluxes = [tjlf_sol.PARTICLE_FLUX_e, tjlf_sol.ENERGY_FLUX_e,
                              tjlf_sol.ENERGY_FLUX_i, tjlf_sol.STRESS_TOR_i]
                for (k, (flux, tjlf_val)) in enumerate(zip(fluxes, tjlf_fluxes))
                    pred_val = Measurements.value(flux)
                    rel_errors_vec[k] = abs(pred_val - tjlf_val) / (abs(tjlf_val) + 1e-6)
                end
            end
            rel_error = sum(rel_errors_vec) / length(rel_errors_vec)

            push!(model_scores, (model_name, confidence, rel_error, rel_errors_vec, flux_sol))
        end

        # Sort by rel_error when ground truth is available, otherwise by confidence
        sort!(model_scores, by=x -> tjlf_sol === nothing ? x[2] : x[3])

        # Deduplicate: group models with identical flux outputs
        flux_key(ms) = (
            Measurements.value(ms[5].PARTICLE_FLUX_e),
            Measurements.value(ms[5].ENERGY_FLUX_e),
            Measurements.value(ms[5].ENERGY_FLUX_i),
            Measurements.value(ms[5].STRESS_TOR_i)
        )
        seen_fluxes = Dict{NTuple{4,Float64}, Vector{String}}()
        unique_model_scores = []
        for ms in model_scores
            key = flux_key(ms)
            if !haskey(seen_fluxes, key)
                seen_fluxes[key] = [ms[1]]
                push!(unique_model_scores, ms)
            else
                push!(seen_fluxes[key], ms[1])
            end
        end

        n_top = min(max_models, length(unique_model_scores))

        push!(rankings, (
            rho=rho,
            top_models=[ms[1] for ms in unique_model_scores[1:n_top]],
            confidences=[ms[2] for ms in unique_model_scores[1:n_top]],
            rel_errors=[ms[3] for ms in unique_model_scores[1:n_top]],
            rel_errors_vec=[ms[4] for ms in unique_model_scores[1:n_top]],
            flux_outputs=[ms[5] for ms in unique_model_scores[1:n_top]]
        ))

        if verbose
            flux_str(Γ, Qe, Qi, Π) = "Γ=$(round(Γ, digits=3)) Qe=$(round(Qe, digits=3)) Qi=$(round(Qi, digits=3)) Π=$(round(Π, digits=3))"
            gt_flux_str = if tjlf_sol !== nothing && show_fluxes
                " [$(ground_truth_label): $(flux_str(tjlf_sol.PARTICLE_FLUX_e, tjlf_sol.ENERGY_FLUX_e, tjlf_sol.ENERGY_FLUX_i, tjlf_sol.STRESS_TOR_i))]"
            else
                ""
            end
            println("\nTop $n_top models at RMIN_LOC=$(round(rho, digits=3)):$gt_flux_str")
            for (i, ms) in enumerate(unique_model_scores[1:n_top])
                (name, conf, rel_err, rel_errs, flux_sol) = ms
                all_names = join(seen_fluxes[flux_key(ms)], ", ")
                err_str = if tjlf_sol !== nothing
                    if show_fluxes
                        errs = join([round(r, digits=3) for r in rel_errs], " ")
                        ", rel_err_$(ground_truth_label)=[$errs]"
                    else
                        ", rel_err_$(ground_truth_label)=$(round(rel_err, digits=4))"
                    end
                else
                    ""
                end
                nn_flux_suffix = if show_fluxes
                    " [NN: $(flux_str(
                        Measurements.value(flux_sol.PARTICLE_FLUX_e),
                        Measurements.value(flux_sol.ENERGY_FLUX_e),
                        Measurements.value(flux_sol.ENERGY_FLUX_i),
                        Measurements.value(flux_sol.STRESS_TOR_i)
                    ))]"
                else
                    ""
                end
                println("  $i. $all_names (conf=$(round(conf, digits=4))$err_str)$nn_flux_suffix")
            end
        end
    end

    return (
        rho_grid=rho_vec,
        rankings=rankings,
        input_tglfs=input_tglfs,
        all_results=all_results,
        tjlf_sols=tjlf_sols,
        tglf_sols=tglf_sols
    )
end

# ------------------------------------------------------------
# Model Selector - Public API with Multiple Dispatch
# ------------------------------------------------------------

"""
    model_selector(ods_path::String; rho_grid=range(0.1, 0.9, 9), electromagnetic=true, lump_ions=false, sat_rule=:sat3, filter_sat_rule=nothing, max_models=3, verbose=true)

Select the most confident TGLF-NN models for a given ODS file across a radial grid.

This function:
1. Loads an ODS (IMAS data structure) from the given path
2. Generates InputTGLF files for the specified rho grid
3. Runs available TGLF-NN models (with graceful failure handling)
4. Ranks models by confidence using the relative uncertainty metric
5. Returns the top `max_models` most confident models for each rho location

Note: Automatically skips "gknn", "qlnn", and "edge" models (correction/specialized models, not standalone)

# Arguments
- `ods_path::String`: Path to the ODS file (JSON or IMAS format)
- `rho_grid`: Radial grid points - can be a range, array, tuple, or single value (default: range(0.1, 0.9, 9))
- `electromagnetic::Bool`: Use electromagnetic (true, default) vs electrostatic (false). Controls both InputTGLF generation and model filtering.
- `lump_ions::Bool`: Lump ion species together (default: false)
- `sat_rule::Symbol`: Saturation rule for InputTGLF generation (default: :sat3)
- `filter_sat_rule::Union{Symbol,Nothing}`: Filter models by saturation rule. By default (nothing), automatically matches `sat_rule`. Set to a specific symbol (e.g., :sat1) to override, or set to :all to test all models. (default: nothing, auto-matches sat_rule)
- `max_models::Int`: Number of top models to return per rho (default: 3)
- `ground_truth::Bool`: Run TJLF/TGLF for accuracy-based ranking (default: true). Set to false for faster confidence-only ranking.
- `show_fluxes::Bool`: Print per-flux values and vectorized relative errors (default: false).
- `verbose::Bool`: Print progress information (default: true)

# Returns
A NamedTuple with:
- `rho_grid`: The radial grid used
- `rankings`: Vector of NamedTuples (one per rho) containing:
  - `rho`: The rho value
  - `top_models`: Vector of model names (most confident first)
  - `confidences`: Vector of confidence scores (lower is better)
  - `flux_outputs`: Vector of flux solutions for top models
- `input_tglfs`: Vector of InputTGLF structures used
- `all_results`: Dict mapping model_name => (success, result/error) for all models

# Example
```julia
using TurbulentTransport, IMAS

# Load and analyze an ODS file - defaults to sat3, auto-matches filter to sat3
results = model_selector("/path/to/ods.json"; rho_grid=range(0.1, 0.9, 9))

# Single radial location (scalar, array, or tuple all work)
results = model_selector("/path/to/ods.json"; rho_grid=0.5)
results = model_selector("/path/to/ods.json"; rho_grid=[0.5])
results = model_selector("/path/to/ods.json"; rho_grid=(0.5,))

# Use sat1 - filter automatically matches
results = model_selector("/path/to/ods.json"; sat_rule=:sat1)

# Test all models regardless of sat rule
results = model_selector("/path/to/ods.json"; filter_sat_rule=:all)

# Check top models at first rho location
println("Top models at ρ=\$(results.rankings[1].rho):")
for (i, (model, conf)) in enumerate(zip(results.rankings[1].top_models, results.rankings[1].confidences))
    println("  \$i. \$model (confidence: \$(round(conf, digits=4)))")
end
```
"""
function model_selector(ods_path::String;
                       rho_grid=range(0.1, 0.9, 9),
                       electromagnetic::Bool=true,
                       lump_ions::Bool=false,
                       sat_rule::Symbol=:sat3,
                       filter_sat_rule::Union{Symbol,Nothing}=nothing,
                       max_models::Int=3,
                       ground_truth::Bool=true,
                       show_fluxes::Bool=false,
                       verbose::Bool=true)

    verbose && println("Loading ODS from: $ods_path")
    dd = IMAS.json2imas(ods_path; error_on_missing_coordinates=false, show_warnings=false)

    # Convert rho_grid to vector (handles scalar, tuple, range, etc.)
    rho_vec = rho_grid isa Number ? [Float64(rho_grid)] : collect(rho_grid)

    verbose && println("Generating InputTGLF for $(length(rho_vec)) radial locations\n")
    input_tglfs = _unwrap_input_tglfs(InputTGLF(dd, rho_vec, sat_rule, electromagnetic, lump_ions))

    # Auto-match filter_sat_rule to sat_rule if not specified
    effective_filter = filter_sat_rule === nothing ? sat_rule : filter_sat_rule

    return _model_selector_core(input_tglfs, rho_vec; max_models, filter_sat_rule=effective_filter, electromagnetic, ground_truth, show_fluxes, verbose)
end

"""
    model_selector(dd::IMAS.dd; rho_grid=range(0.1, 0.9, 9), electromagnetic=true, lump_ions=false, sat_rule=:sat3, filter_sat_rule=nothing, max_models=3, verbose=true)

Select the most confident TGLF-NN models for a given IMAS data structure across a radial grid.

Similar to the String path method, but takes an already-loaded IMAS.dd object (e.g., from FUSE initialization).

Note: Automatically skips "gknn", "qlnn", and "edge" models (correction/specialized models, not standalone)

# Arguments
- `dd::IMAS.dd`: IMAS data structure (already loaded/initialized)
- `rho_grid`: Radial grid points - can be a range, array, tuple, or single value (default: range(0.1, 0.9, 9))
- `electromagnetic::Bool`: Use electromagnetic (true, default) vs electrostatic (false). Controls both InputTGLF generation and model filtering.
- `lump_ions::Bool`: Lump ion species together (default: false)
- `sat_rule::Symbol`: Saturation rule for InputTGLF generation (default: :sat3)
- `filter_sat_rule::Union{Symbol,Nothing}`: Filter models by saturation rule. By default (nothing), automatically matches `sat_rule`. Set to a specific symbol (e.g., :sat1) to override, or set to :all to test all models. (default: nothing, auto-matches sat_rule)
- `max_models::Int`: Number of top models to return per rho (default: 3)
- `ground_truth::Bool`: Run TJLF/TGLF for accuracy-based ranking (default: true). Set to false for faster confidence-only ranking.
- `show_fluxes::Bool`: Print per-flux values and vectorized relative errors (default: false).
- `verbose::Bool`: Print progress information (default: true)

# Returns
Same structure as the String-based method (see [`model_selector(::String)`](@ref))

# Example
```julia
using TurbulentTransport, IMAS, FUSE

# Initialize with FUSE
ini, act = FUSE.case_parameters(:D3D, :L_mode)
dd = IMAS.dd()
FUSE.init(dd, ini, act)

# Run model selector on the initialized dd
results = model_selector(dd; rho_grid=[0.3, 0.5, 0.7])

# Or use different sat rule
results = model_selector(dd; rho_grid=0.5, sat_rule=:sat1)
```
"""
function model_selector(dd::IMAS.dd;
                       rho_grid=range(0.1, 0.9, 9),
                       electromagnetic::Bool=true,
                       lump_ions::Bool=false,
                       sat_rule::Symbol=:sat3,
                       filter_sat_rule::Union{Symbol,Nothing}=nothing,
                       max_models::Int=3,
                       ground_truth::Bool=true,
                       show_fluxes::Bool=false,
                       verbose::Bool=true)

    rho_vec = rho_grid isa Number ? [Float64(rho_grid)] : collect(rho_grid)

    verbose && println("Generating InputTGLF for $(length(rho_vec)) radial locations\n")
    input_tglfs = _unwrap_input_tglfs(InputTGLF(dd, rho_vec, sat_rule, electromagnetic, lump_ions))

    effective_filter = filter_sat_rule === nothing ? sat_rule : filter_sat_rule

    return _model_selector_core(input_tglfs, rho_vec; max_models, filter_sat_rule=effective_filter, electromagnetic, ground_truth, show_fluxes, verbose)
end

"""
    model_selector(input_tglfs::Vector{InputTGLF{T}}; filter_sat_rule=:sat3, electromagnetic=true, max_models=3, verbose=true) where T

Select the most confident TGLF-NN models for a vector of InputTGLF objects.

This function:
1. Takes pre-generated InputTGLF objects
2. Runs available TGLF-NN models (with graceful failure handling)
3. Ranks models by confidence using the relative uncertainty metric
4. Returns the top `max_models` most confident models for each input

# Arguments
- `input_tglfs::Vector{InputTGLF{T}}`: Vector of InputTGLF structures (already contain all radial location info)
- `filter_sat_rule::Union{Symbol,Nothing}`: Filter models by saturation rule. Default is :sat3. Set to :all to test all models. (default: :sat3)
- `electromagnetic::Bool`: Filter to electromagnetic (true, default) or electrostatic (false) models
- `max_models::Int`: Number of top models to return per input (default: 3)
- `ground_truth::Bool`: Run TJLF/TGLF for accuracy-based ranking (default: true). Set to false for faster confidence-only ranking.
- `show_fluxes::Bool`: Print per-flux values and vectorized relative errors (default: false).
- `verbose::Bool`: Print progress information (default: true)

# Returns
Same structure as the String-based method (see [`model_selector(::String)`](@ref)).

# Example
#=julia
using TurbulentTransport, IMAS

# Pre-generate InputTGLF structures with sat3
dd = IMAS.json2imas("/path/to/ods.json")
rho_grid = [0.3, 0.5, 0.7]
input_tglfs = InputTGLF(dd, rho_grid, :sat3, true, false)

# Run model selector
results = model_selector(input_tglfs)

# Or filter to sat1 models
input_tglfs_sat1 = InputTGLF(dd, rho_grid, :sat1, true, false)
results = model_selector(input_tglfs_sat1; filter_sat_rule=:sat1)
=#
"""
function model_selector(input_tglfs::Vector{InputTGLF{T}};
                       filter_sat_rule::Union{Symbol,Nothing}=:sat3,
                       electromagnetic::Bool=true,
                       max_models::Int=3,
                       ground_truth::Bool=true,
                       show_fluxes::Bool=false,
                       verbose::Bool=true) where T

    verbose && println("Analyzing $(length(input_tglfs)) InputTGLF objects\n")
    rho_values = [inp.RMIN_LOC for inp in input_tglfs]
    return _model_selector_core(input_tglfs, rho_values; max_models, filter_sat_rule, electromagnetic, ground_truth, show_fluxes, verbose)
end

"""
    model_selector(input_tglf::InputTGLF{T}; filter_sat_rule=:sat3, electromagnetic=true, max_models=3, verbose=true) where T

Select the most confident TGLF-NN models for a single InputTGLF object.

Convenience method that wraps a single InputTGLF in a vector and extracts the single ranking result.

# Arguments
- `input_tglf::InputTGLF{T}`: Single InputTGLF structure (already contains all radial location info)
- `filter_sat_rule::Union{Symbol,Nothing}`: Filter models by saturation rule. Default is :sat3. Set to :all to test all models. (default: :sat3)
- `electromagnetic::Bool`: Filter to electromagnetic (true, default) or electrostatic (false) models
- `max_models::Int`: Number of top models to return (default: 3)
- `ground_truth::Bool`: Run TJLF/TGLF for accuracy-based ranking (default: true). Set to false for faster confidence-only ranking.
- `show_fluxes::Bool`: Print per-flux values and vectorized relative errors (default: false).
- `verbose::Bool`: Print progress information (default: true)

# Returns
A NamedTuple with:
- `rho`: The rho value (or 1 if not provided)
- `top_models`: Vector of model names (most confident first)
- `confidences`: Vector of confidence scores (lower is better)
- `flux_outputs`: Vector of flux solutions for top models
- `input_tglf`: The InputTGLF structure used
- `all_results`: Dict mapping model_name => (success, result/error) for all models

# Example
#=julia
using TurbulentTransport, IMAS

# Create a single InputTGLF with sat3
dd = IMAS.json2imas("/path/to/ods.json")
input_tglf = InputTGLF(dd, [0.5], :sat3, true, false)[1]

# Find best models
result = model_selector(input_tglf)

println("Top model: \$(result.top_models[1])")
println("Confidence: \$(result.confidences[1])")
=#
"""
function model_selector(input_tglf::InputTGLF{T};
                       filter_sat_rule::Union{Symbol,Nothing}=:sat3,
                       electromagnetic::Bool=true,
                       max_models::Int=3,
                       ground_truth::Bool=true,
                       show_fluxes::Bool=false,
                       verbose::Bool=true) where T

    verbose && println("Analyzing single InputTGLF\n")
    rho_values = [input_tglf.RMIN_LOC]
    full_results = _model_selector_core([input_tglf], rho_values; max_models, filter_sat_rule, electromagnetic, ground_truth, show_fluxes, verbose)

    ranking = full_results.rankings[1]
    return (rho=ranking.rho, top_models=ranking.top_models, confidences=ranking.confidences,
            rel_errors=ranking.rel_errors, rel_errors_vec=ranking.rel_errors_vec,
            flux_outputs=ranking.flux_outputs,
            input_tglf=input_tglf, all_results=full_results.all_results,
            tjlf_sols=full_results.tjlf_sols, tglf_sols=full_results.tglf_sols)
end

# ------------------------------------------------------------
# Helper: species splitting transform for stfpp models
# Replicates training-time dictionary manipulation for inference.
# dtf: deuterium-tritium fraction assigned to new AS_2 (remainder to AS_3)
# device: optional device string ("ukstep" -> NS=4 else NS=5)
function _apply_stfpp_transform!(t::InputTGLF; dtf::Float64=0.5, device::AbstractString="")
    # This mirrors the training-time dictionary manipulation exactly:
    # 1. Rename *_4 -> *_5, *_3 -> *_4, drop original *_5
    # 2. Duplicate *_2 -> *_3
    # 3. Split AS_2 into AS_2 (dtf) and AS_3 (1-dtf)
    # 4. Set MASS_2, MASS_3, NS
    # Skip if already unbundled (MASS_3 already near target ~1.49760 within 1%)
    try
        m3 = getfield(t, :MASS_3)
        if !(m3 === missing) && isfinite(m3) && abs(m3 - 1.49760)/1.49760 < 0.01
            return t
        end
    catch
        # ignore if field access fails
    end
    Tt = typeof(t)
    orig = Dict{String,Any}(String(f) => getfield(t,f) for f in fieldnames(Tt))
    temp = Dict{String,Any}()
    for (key,val) in orig
        if endswith(key, "_4")
            temp[replace(key, "_4" => "_5")] = val
        elseif endswith(key, "_3")
            temp[replace(key, "_3" => "_4")] = val
        elseif endswith(key, "_5")
            continue  # drop original _5
        else
            temp[key] = val
        end
    end
    # Duplicate *_2 -> *_3
    for (key,val) in collect(temp)
        if occursin("_2", key)
            temp[replace(key, "_2" => "_3")] = val
        end
    end
    # Remember original AS_2 before splitting (from temp now)
    if haskey(temp, "AS_2") && temp["AS_2"] !== missing
        original_as_2 = temp["AS_2"]
        temp["AS_2"] = original_as_2 * dtf
        temp["AS_3"] = original_as_2 * (1 - dtf)
    end
    # Set masses per spec
    temp["MASS_2"] = 1.0
    temp["MASS_3"] = 1.49760170089
    temp["NS"] = device == "ukstep" ? 4 : 5
    # Write back to struct fields (missing if dropped)
    for f in fieldnames(Tt)
        fname = String(f)
        if haskey(temp, fname)
            val = temp[fname]
            try
                setfield!(t, f, val)
            catch
                # Ignore type mismatch silently
            end
        else
            # If this was a dropped *_5 (original) ensure it's missing
            if endswith(fname, "_5")
                try; setfield!(t, f, missing); catch; end
            end
        end
    end
    return t
end

function _apply_stfpp_transform!(data::Dict; dtf::Float64=0.5, device::AbstractString="")
    # Assume data values are vectors (as in run_tglfnn(dict))
    # Skip if already unbundled: check MASS_3 first element (or scalar) ~ 1.49760 within 1%
    if haskey(data, "MASS_3")
        m3val = data["MASS_3"]
        m3 = try
            isa(m3val, AbstractArray) ? m3val[begin] : m3val
        catch
            nothing
        end
        if m3 isa Number && isfinite(m3) && abs(m3 - 1.49760)/1.49760 < 0.01
            return data
        end
    end
    tempdict = Dict{String,Any}()
    for (key,val) in data
        if endswith(key, "_4")
            tempdict[replace(key, "_4" => "_5")] = val
        elseif endswith(key, "_3")
            tempdict[replace(key, "_3" => "_4")] = val
        elseif endswith(key, "_5")
            # skip
        else
            tempdict[key] = val
        end
    end
    # Duplicate _2 -> _3
    for (key,val) in collect(tempdict)
        if occursin("_2", key)
            tempdict[replace(key, "_2" => "_3")] = val
        end
    end
    if haskey(tempdict, "AS_2")
        original_as_2 = tempdict["AS_2"]
        tempdict["AS_2"] = original_as_2 .* dtf
        tempdict["AS_3"] = original_as_2 .* (1 - dtf)
    end
    N = haskey(tempdict, "AS_2") ? length(tempdict["AS_2"]) : (haskey(tempdict, "AS_1") ? length(tempdict["AS_1"]) : 0)
    if N > 0
        tempdict["MASS_2"] = fill(1.0, N)
        tempdict["MASS_3"] = fill(1.49760170089, N)
        tempdict["NS"] = fill(device == "ukstep" ? 4 : 5, N)
    else
        tempdict["MASS_2"] = 1.0
        tempdict["MASS_3"] = 1.49760170089
        tempdict["NS"] = device == "ukstep" ? 4 : 5
    end
    # Write back
    empty!(data)
    for (k,v) in tempdict
        data[k] = v
    end
    return data
end