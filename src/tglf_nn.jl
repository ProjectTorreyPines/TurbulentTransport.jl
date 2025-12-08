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

# External cache for field symbols Val types (avoids Any field in struct)
const _XNAMES_FIELD_SYMBOLS_CACHE = Dict{UInt64, Any}()

"""
    _get_xnames_without_log10_suffix(model::TGLFNNmodel)

Get cached `Val{Tuple{Symbol...}}` of InputTGLF field names from model's xnames.
Strips `_log10` suffix: `"BETAE_log10"` → `:BETAE`

Uses external cache to avoid type instability from storing in struct.
Used by `_extract_fields!` for zero-allocation field extraction.
"""
function _get_xnames_without_log10_suffix(model::TGLFNNmodel)
    key = objectid(model.xnames)
    get!(_XNAMES_FIELD_SYMBOLS_CACHE, key) do
        symbols = Tuple(Symbol(endswith(x, log_suffix) ? x[1:end-n_log_suffix] : x) for x in model.xnames)
        Val(symbols)
    end
end

function _get_xnames_without_log10_suffix(ensemble::TGLFNNensemble)
    _get_xnames_without_log10_suffix(ensemble.models[1])
end

"""
    _extract_fields!(inputs::AbstractVector, obj, ::Val{symbols}, index::Int=0) where {symbols}

Zero-allocation field extraction with `ismissing` validation.
`symbols` is a tuple of field names (Symbols) known at compile time.
`index` is the radial location index for error messages (default 0 for single InputTGLF).

This function is `@generated` to unroll field access at compile time,
avoiding the boxing overhead of dynamic `getfield(obj, runtime_symbol)`.
The `ismissing` check adds negligible overhead due to branch prediction.
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
    if !endswith(filename, ".bson")
        filename = "$(filename).bson"
    end
    if startswith(filename, "/")
        fullpath = filename
    else
        fullpath = dirname(@__DIR__) * "/models/" * filename
        if !isfile(fullpath)
            error("TGLFNN model $filename does not exist. Possible nn models are:\n    $(join(available_models(),"\n    ",))")
        end
    end
    savedict = BSON.load(fullpath, @__MODULE__)
    if typeof(first(keys(savedict))) <: Integer
        return dict2ens(savedict)
    else
        return dict2mod(savedict)
    end
end

function available_models()
    models_dir = joinpath(dirname(@__DIR__), "models")
    return [replace(model, r"\.(bson|onnx)$" => "") for model in readdir(models_dir) if endswith(model, ".bson") || endswith(model, ".onnx")]
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
function flux_array(fluxmodel::TGLFNNmodel, x::AbstractVector{T}; warn_nn_train_bounds::Bool=true, fidelity::Symbol=:TGLFNN, xx::AbstractVector{T}=similar(x)) where {T<:Real}
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

    # unsafe_acquire! returns Array (not ReshapedArray) to avoid boxing with non-concrete _pooled_chain
    xx = unsafe_acquire!(pool, T, size(x))
    
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

    # unsafe_acquire! returns Array (not ReshapedArray) to avoid boxing with non-concrete _pooled_chain
    xx = unsafe_acquire!(pool, T, N)

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
    each_y = unsafe_acquire!(pool, T, nouts, nsamples)
    all_yy = unsafe_acquire!(pool, T, nouts, nsamples, nmodels)

    # Threads.@threads for k in 1:nmodels
    for k in 1:nmodels
        # tmp[:, :, k] = flux_array(fluxensemble.models[k], x; warn_nn_train_bounds=(warn_nn_train_bounds && k == 1), fidelity)
        flux_array!(each_y, fluxensemble.models[k], x; warn_nn_train_bounds=(warn_nn_train_bounds && k == 1), fidelity)
        all_yy[:, :, k] = each_y
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
    each_y = unsafe_acquire!(pool, T, nouts)
    all_yy = unsafe_acquire!(pool, T, nouts, nmodels)
    for k in 1:nmodels
        flux_array!(each_y, fluxensemble.models[k], x; warn_nn_train_bounds=(warn_nn_train_bounds && k == 1), fidelity)
        all_yy[:, k] = each_y
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
    if occursin("stfpp", model_filename)
        for it in input_tglfs
            _apply_stfpp_transform!(it; dtf=0.5, device="")
        end
    end
    if model_filename in ("sat3_em_d3d_azf-1",) && fidelity == :GKNN
        tglfmod = loadmodelonce(model_filename * "_tglfnn24")
    else
        tglfmod = loadmodelonce(model_filename)
    end

    inputs = acquire!(pool, T, length(tglfmod.xnames), length(input_tglfs))

    # Extract input fields using @generated function for zero-allocation
    xnames_val = _get_xnames_without_log10_suffix(tglfmod)
    for (i, input_tglf) in enumerate(input_tglfs)
        _extract_fields!(@view(inputs[:, i]), input_tglf, xnames_val)
    end

    tmp = flux_array(tglfmod, inputs; uncertain, warn_nn_train_bounds, fidelity=:TGLFNN)
    if fidelity == :GKNN
        supported_gknn_models = ("sat3_em_d3d_azf-1", "sat3_em_d3d+mastu+nstx_azf-1", "sat3_em_d3d_azf-1_withnegD", "sat3_em_d3d_azf-1_gkdb", "sat2_em_d3d+mastu+nstx_azf-1", "sat3_em_d3d+mastu_azf-1")
        if !(model_filename in supported_gknn_models)
            error("GKNN fidelity is not supported for model '$model_filename'. Supported models are: $(join(supported_gknn_models, ", "))")
        end

        if model_filename == "sat3_em_d3d_azf-1"
            gk_inputs = acquire!(pool, T, size(inputs, 1) + 1, size(inputs, 2))
            gk_inputs[1:end-1, :] = inputs

            for (i, postfix) in enumerate(("_gknng24", "_gknnp24", "_gknne24", "_gknni24"))
                gk_inputs[end, :] = tmp[i, :]
                gknn_model = loadmodelonce(model_filename * postfix)
                err = flux_array(gknn_model, gk_inputs; uncertain, warn_nn_train_bounds, fidelity)[:]
                tmp[i, :] .*= err
            end
        elseif model_filename in ("sat3_em_d3d+mastu+nstx_azf-1", "sat3_em_d3d_azf-1_withnegD", "sat3_em_d3d_azf-1_gkdb", "sat2_em_d3d+mastu+nstx_azf-1")
            gk_inputs = acquire!(pool, T, size(inputs, 1) + 4, size(inputs, 2))
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
            gk_inputs = acquire!(pool, T, size(inputs, 1) + 4, size(inputs, 2))
            gk_inputs[1:end-4, :] = inputs
            gk_inputs[end-3:end, :] = tmp

            gknn = loadmodelonce(model_filename * "_gknn36")
            err = flux_array(gknn, gk_inputs; uncertain, warn_nn_train_bounds, fidelity)
            tmp .*= err
        end
    end

    sol = [flux_solution(tmp[:, i]...) for i in eachindex(input_tglfs)]
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
    if occursin("stfpp", model_filename)
        _apply_stfpp_transform!(data; dtf=0.5, device="")
    end
    if model_filename in ("sat3_em_d3d_azf-1",) && fidelity == :GKNN
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
        if model_filename in ("sat3_em_d3d_azf-1",)
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
        elseif model_filename in ("sat3_em_d3d+mastu+nstx_azf-1", "sat3_em_d3d_azf-1_withnegD", "sat3_em_d3d_azf-1_gkdb", "sat2_em_d3d+mastu+nstx_azf-1")
            gknn = loadmodelonce(model_filename * "_gknn31")
            err = gknn(vcat(x, y)...; uncertain, warn_nn_train_bounds, fidelity)
            y .*= err
            if model_filename in ("sat3_em_d3d_azf-1_gkdb",)
                gkdb = loadmodelonce(model_filename * "_gknn31_cgyro")
                gkdb_err = gkdb(vcat(x, y)...; uncertain, warn_nn_train_bounds, fidelity)
                y .*= gkdb_err
            end
        elseif model_filename in ["sat3_em_d3d+mastu_azf-1"]
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
    if !occursin("/models/", onnx_path)
        onnx_path = joinpath(dirname(@__DIR__), "models",
                             endswith(onnx_path, ".onnx") ? onnx_path : onnx_path * ".onnx")
    end
    isfile(onnx_path) || error("TGLFNN model does not exist in $onnx_path")
    return onnx_path
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

Constructor used to handle PARTICLE_FLUX_i entered as a set of scalars instead of an array

    flux_solution(1.0, 2.0, 3.0, 4.0, 5.0, 6.0)

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
function flux_solution(xx::Vararg{T}) where {T<:Real}
    n_fields = length(xx)
    if n_fields == 4
        ENERGY_FLUX_e = 3
        ENERGY_FLUX_i = 4
        PARTICLE_FLUX_e = 1
        STRESS_TOR_i = 2
        sol = GACODE.FluxSolution{T}(xx[ENERGY_FLUX_e], xx[ENERGY_FLUX_i], xx[PARTICLE_FLUX_e], T[], xx[STRESS_TOR_i])
    else
        ENERGY_FLUX_e = n_fields - 1
        ENERGY_FLUX_i = n_fields
        PARTICLE_FLUX_e = 1
        PARTICLE_FLUX_i = 2:n_fields-3
        STRESS_TOR_i = n_fields - 2
        sol = GACODE.FluxSolution{T}(xx[ENERGY_FLUX_e], xx[ENERGY_FLUX_i], xx[PARTICLE_FLUX_e], T[xx[i] for i in PARTICLE_FLUX_i], xx[STRESS_TOR_i])
    end
    return sol
end

export run_tglfnn, run_tglfnn_onnx

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