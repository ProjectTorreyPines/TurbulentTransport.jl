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

#= ============== =#
#  saving/loading  #
#= ============== =#
function mod2dict(model::TGLFNNmodel)
    savedict = Dict()
    for name in fieldnames(TGLFNNmodel)
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
function flux_array(fluxmodel::TGLFNNmodel, x::AbstractMatrix{T}; warn_nn_train_bounds::Bool=true, fidelity::Symbol=:TGLFNN) where {T}
    N, _ = size(x)
    xx = Vector{T}(undef, N)
    return hcat(collect(map(x0 -> flux_array(fluxmodel, x0; warn_nn_train_bounds, fidelity, xx), eachslice(x; dims=2)))...)
end

function flux_array(fluxmodel::TGLFNNmodel, x::AbstractVector{T}; warn_nn_train_bounds::Bool=true, fidelity::Symbol=:TGLFNN, xx::Vector{T}=similar(x)) where {T}
    for (ix, name) in enumerate(fluxmodel.xnames)
        xx[ix] = contains(name, "_log10") ? log10(x[ix]) : x[ix]
    end
    if warn_nn_train_bounds # training bounds are on the original data but after log10
        for ix in eachindex(xx)
            if isnan(xx[ix]) || isinf(xx[ix])
                error("$(fluxmodel.xnames[ix]) = $(x[ix]) is not allowed")
            elseif xx[ix] < fluxmodel.xbounds[ix, 1]
                @warn("Extrapolation $(fluxmodel.xnames[ix])=$(minimum(xx[ix,:])) is below training bound of $(fluxmodel.xbounds[ix,:])")
            elseif xx[ix] > fluxmodel.xbounds[ix, 2]
                @warn("Extrapolation $(fluxmodel.xnames[ix])=$(maximum(xx[ix,:])) is above training bound of $(fluxmodel.xbounds[ix,:])")
            end
        end
    end
    @. xx = (xx - fluxmodel.xm) / fluxmodel.xσ
    yy = fluxmodel.fluxmodel(xx)
    if fidelity == :GKNN
        return yy
    elseif fidelity == :TGLFNN
        @. yy = yy * fluxmodel.yσ + fluxmodel.ym
        return yy
    end
end

function flux_array(fluxensemble::TGLFNNensemble, x::AbstractArray{T}; uncertain::Bool=false, warn_nn_train_bounds::Bool=true, fidelity::Symbol=:TGLFNN) where {T<:Real}
    nmodels = length(fluxensemble.models)
    nouts = length(fluxensemble.models[1].ynames)
    if fidelity == :GKNN
        nouts = div(nouts, 2)
    end
    nsamples = size(x)[2]

    tmp = zeros(T, nmodels, nouts, nsamples)
    Threads.@threads for k in 1:length(fluxensemble.models)
        tmp[k, :, :] = flux_array(fluxensemble.models[k], x; warn_nn_train_bounds=(warn_nn_train_bounds && k == 1), fidelity)
    end

    mean, std = StatsBase.mean_and_std(tmp, 1; corrected=true)
    if uncertain && nmodels > 1
        if T <: Measurements.Measurement
            return mean[1, :, :]
        else
            return Measurements.measurement.(mean[1, :, :], std[1, :, :])
        end
    else
        return mean[1, :, :]
    end
end

function flux_array(fluxensemble::TGLFNNensemble, x::AbstractVector{T}; uncertain::Bool=false, warn_nn_train_bounds::Bool=true, fidelity::Symbol=:TGLFNN) where {T<:Real}
    nmodels = length(fluxensemble.models)
    nouts = length(fluxensemble.models[1].ynames)
    if fidelity == :GKNN
        nouts = div(nouts, 2)
    end

    tmp = zeros(T, nmodels, nouts)
    Threads.@threads for k in 1:length(fluxensemble.models)
        tmp[k, :] = flux_array(fluxensemble.models[k], x; warn_nn_train_bounds=(warn_nn_train_bounds && k == 1), fidelity)
    end

    mean, std = StatsBase.mean_and_std(tmp, 1; corrected=true)
    if uncertain
        if T <: Measurements.Measurement
            return mean[1, :]
        else
            return Measurements.measurement.(mean[1, :], std[1, :])
        end
    else
        return mean[1, :]
    end
end

function flux_array(fluxmodel::TGLFmodel, args...; uncertain::Bool=false, warn_nn_train_bounds::Bool=true, fidelity::Symbol=:TGLFNN)
    args = reshape([k for k in args], (length(args), 1))
    return flux_array(fluxmodel, args; uncertain, warn_nn_train_bounds, fidelity)
end

function flux_solution(fluxmodel::TGLFmodel, args...; uncertain::Bool=false, warn_nn_train_bounds::Bool=true, fidelity::Symbol=:TGLFNN)
    return flux_solution(flux_array(fluxmodel, collect(args); uncertain, warn_nn_train_bounds, fidelity)...)
end

#= ======================= =#
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
    run_tglfnn(input_tglf::InputTGLF; model_filename::String, uncertain::Bool=false, warn_nn_train_bounds::Bool, fidelity::Symbol=:TGLFNN)

Run TGLFNN starting from a InputTGLF, using a specific `model_filename`.

If the model is an ensemble of NNs, then the output can be uncertain (using the Measurements.jl package).

The warn_nn_train_bounds checks against the standard deviation of the inputs to warn if evaluation is likely outside of training bounds.

Returns a `flux_solution` structure
"""
function run_tglfnn(input_tglf::InputTGLF{T}; model_filename::String, uncertain::Bool=false, warn_nn_train_bounds::Bool, fidelity::Symbol=:TGLFNN) where {T<:Real}
    # In-place transform for any model containing "stfpp"
    if occursin("stfpp", model_filename)
        _apply_stfpp_transform!(input_tglf; dtf=0.5, device="")
    end
    if model_filename in ["sat3_em_d3d_azf-1"] && fidelity == :GKNN
        tglfmod = loadmodelonce(model_filename * "_tglfnn24")
    else
        tglfmod = loadmodelonce(model_filename)
    end
    inputs = zeros(length(tglfmod.xnames))
    for (k, item) in enumerate(tglfmod.xnames)
        item = replace(item, "_log10" => "")
        inputs[k] = getfield(input_tglf, Symbol(item))
    end
    sol = tglfmod(inputs...; uncertain, warn_nn_train_bounds, fidelity=:TGLFNN)
    if fidelity == :GKNN
        supported_gknn_models = ["sat3_em_d3d_azf-1", "sat3_em_d3d+mastu+nstx_azf-1", "sat3_em_d3d_azf-1_withnegD", "sat3_em_d3d_azf-1_gkdb", "sat2_em_d3d+mastu+nstx_azf-1", "sat3_em_d3d+mastu_azf-1"]
        if !(model_filename in supported_gknn_models)
            error("GKNN fidelity is not supported for model '$model_filename'. Supported models are: $(join(supported_gknn_models, ", "))")
        end
        base_fluxes = [sol.ENERGY_FLUX_e, sol.ENERGY_FLUX_i, sol.PARTICLE_FLUX_e, sol.STRESS_TOR_i]
        if model_filename in ["sat3_em_d3d_azf-1"]
            gknne = loadmodelonce(model_filename * "_gknne24")
            err_e = flux_array(gknne, vcat(inputs, base_fluxes[1]); uncertain, warn_nn_train_bounds, fidelity)
            gknni = loadmodelonce(model_filename * "_gknni24")
            err_i = flux_array(gknni, vcat(inputs, base_fluxes[2]); uncertain, warn_nn_train_bounds, fidelity)
            gknng = loadmodelonce(model_filename * "_gknng24")
            err_g = flux_array(gknng, vcat(inputs, base_fluxes[3]); uncertain, warn_nn_train_bounds, fidelity)
            gknnp = loadmodelonce(model_filename * "_gknnp24")
            err_p = flux_array(gknnp, vcat(inputs, base_fluxes[4]); uncertain, warn_nn_train_bounds, fidelity)
            sol = GACODE.FluxSolution{T}(
                base_fluxes[1] * err_e[1],  # ENERGY_FLUX_e
                base_fluxes[2] * err_i[1],  # ENERGY_FLUX_i
                base_fluxes[3] * err_g[1],  # PARTICLE_FLUX_e
                Float64[],                   # PARTICLE_FLUX_i (empty for this model)
                base_fluxes[4] * err_p[1]   # STRESS_TOR_i
            )
        elseif model_filename in ["sat3_em_d3d+mastu+nstx_azf-1", "sat3_em_d3d_azf-1_withnegD", "sat3_em_d3d_azf-1_gkdb", "sat2_em_d3d+mastu+nstx_azf-1"]
            gknn = loadmodelonce(model_filename * "_gknn31")
            err = flux_array(gknn, vcat(inputs, [base_fluxes[3], base_fluxes[4], base_fluxes[1], base_fluxes[2]]); uncertain, warn_nn_train_bounds, fidelity)
            sol = GACODE.FluxSolution{T}(
                sol.ENERGY_FLUX_e * err[3],
                sol.ENERGY_FLUX_i * err[4],
                sol.PARTICLE_FLUX_e * err[1],
                Float64[],                   # PARTICLE_FLUX_i (empty for this model)
                sol.STRESS_TOR_i * err[2]
            )
            if model_filename in ["sat3_em_d3d_azf-1_gkdb"]
                gkdb = loadmodelonce(model_filename * "_gknn31_cgyro")
                gkdb_err = flux_array(gkdb, vcat(inputs, [sol.PARTICLE_FLUX_e, sol.STRESS_TOR_i, sol.ENERGY_FLUX_e, sol.ENERGY_FLUX_i]); uncertain, warn_nn_train_bounds, fidelity)
                sol = GACODE.FluxSolution{T}(
                    sol.ENERGY_FLUX_e * gkdb_err[3],
                    sol.ENERGY_FLUX_i * gkdb_err[4],
                    sol.PARTICLE_FLUX_e * gkdb_err[1],
                    Float64[],                   # PARTICLE_FLUX_i (empty for this model)
                    sol.STRESS_TOR_i * gkdb_err[2]
                )
            end
        elseif model_filename in ["sat3_em_d3d+mastu_azf-1"]
            gknn = loadmodelonce(model_filename * "_gknn36")
            err = flux_array(gknn, vcat(inputs, [base_fluxes[3], base_fluxes[4], base_fluxes[1], base_fluxes[2]]); uncertain, warn_nn_train_bounds, fidelity)
            sol = GACODE.FluxSolution{T}(
                sol.ENERGY_FLUX_e * err[3],
                sol.ENERGY_FLUX_i * err[4],
                sol.PARTICLE_FLUX_e * err[1],
                Float64[],                   # PARTICLE_FLUX_i (empty for this model)
                sol.STRESS_TOR_i * err[2]
            )
        end
    end
    return sol
end

"""
    run_tglfnn(input_tglfs::Vector{InputTGLF{T}}; model_filename::String, uncertain::Bool=false, warn_nn_train_bounds::Bool, fidelity::Symbol=:TGLFNN) where {T<:Real}

Run TGLFNN for multiple InputTGLF, using a specific `model_filename`.

This is more efficient than running TGLFNN on each individual InputTGLFs.

If the model is an ensemble of NNs, then the output can be uncertain (using the Measurements.jl package).

The warn_nn_train_bounds checks against the standard deviation of the inputs to warn if evaluation is likely outside of training bounds.

Returns a vector of `flux_solution` structures
"""
function run_tglfnn(input_tglfs::Vector{InputTGLF{T}}; model_filename::String, uncertain::Bool=false, warn_nn_train_bounds::Bool, fidelity::Symbol=:TGLFNN) where {T<:Real}
    if occursin("stfpp", model_filename)
        for it in input_tglfs
            _apply_stfpp_transform!(it; dtf=0.5, device="")
        end
    end
    if model_filename in ["sat3_em_d3d_azf-1"] && fidelity == :GKNN
        tglfmod = loadmodelonce(model_filename * "_tglfnn24")
    else
        tglfmod = loadmodelonce(model_filename)
    end
    inputs = zeros(T, length(tglfmod.xnames), length(input_tglfs))
    for (i, input_tglf) in enumerate(input_tglfs)
        for (k, item) in enumerate(tglfmod.xnames)
            if endswith(item, log_suffix)
                subitem = SubString(item, firstindex(item), prevind(item, lastindex(item), n_log_suffix))
                value = getfield(input_tglf, Symbol(subitem))
            else
                value = getfield(input_tglf, Symbol(item))
            end
            if ismissing(value)
                hint = ""
                if occursin("_5", item) || occursin("_6", item)
                    hint = "\n\nHint: Missing species data (species 5). If using a TGLFNN model (e.g. 'stfpp' models), try setting:\n  act.ActorTGLF.lump_ions = false\nto ensure ion species are treated separately rather than lumped together."
                end
                error("TGLFNN input field '$(item)' is Missing at radial location $(i). Check that all required equilibrium and profile data are properly initialized.$(hint)")
            end
            inputs[k, i] = value
        end
    end
    tmp = flux_array(tglfmod, inputs; uncertain, warn_nn_train_bounds, fidelity=:TGLFNN)
    if fidelity == :GKNN
        supported_gknn_models = ["sat3_em_d3d_azf-1", "sat3_em_d3d+mastu+nstx_azf-1", "sat3_em_d3d_azf-1_withnegD", "sat3_em_d3d_azf-1_gkdb", "sat2_em_d3d+mastu+nstx_azf-1", "sat3_em_d3d+mastu_azf-1"]
        if !(model_filename in supported_gknn_models)
            error("GKNN fidelity is not supported for model '$model_filename'. Supported models are: $(join(supported_gknn_models, ", "))")
        end
        if model_filename in ["sat3_em_d3d_azf-1"]
            gknng = loadmodelonce(model_filename * "_gknng24")
            err_g = flux_array(gknng, vcat(inputs, reshape(tmp[1, :], 1, :)); uncertain, warn_nn_train_bounds, fidelity)
            tmp[1, :] .*= err_g[1, :]
            gknnp = loadmodelonce(model_filename * "_gknnp24")
            err_p = flux_array(gknnp, vcat(inputs, reshape(tmp[2, :], 1, :)); uncertain, warn_nn_train_bounds, fidelity)
            tmp[2, :] .*= err_p[1, :]
            gknne = loadmodelonce(model_filename * "_gknne24")
            err_e = flux_array(gknne, vcat(inputs, reshape(tmp[3, :], 1, :)); uncertain, warn_nn_train_bounds, fidelity)
            tmp[3, :] .*= err_e[1, :]
            gknni = loadmodelonce(model_filename * "_gknni24")
            err_i = flux_array(gknni, vcat(inputs, reshape(tmp[4, :], 1, :)); uncertain, warn_nn_train_bounds, fidelity)
            tmp[4, :] .*= err_i[1, :]
        elseif model_filename in ["sat3_em_d3d+mastu+nstx_azf-1", "sat3_em_d3d_azf-1_withnegD", "sat3_em_d3d_azf-1_gkdb", "sat2_em_d3d+mastu+nstx_azf-1"]
            gknn = loadmodelonce(model_filename * "_gknn31")
            err = flux_array(gknn, vcat(inputs, tmp); uncertain, warn_nn_train_bounds, fidelity)
            tmp .*= err
            if model_filename in ["sat3_em_d3d_azf-1_gkdb"]
                gkdb = loadmodelonce(model_filename * "_gknn31_cgyro")
                gkdb_err = flux_array(gkdb, vcat(inputs, tmp); uncertain, warn_nn_train_bounds, fidelity)
                tmp .*= gkdb_err
            end
        elseif model_filename in ["sat3_em_d3d+mastu_azf-1"]
            gknn = loadmodelonce(model_filename * "_gknn36")
            err = flux_array(gknn, vcat(inputs, tmp); uncertain, warn_nn_train_bounds, fidelity)
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
    if model_filename in ["sat3_em_d3d_azf-1"] && fidelity == :GKNN
        tglfmod = loadmodelonce(model_filename * "_tglfnn24")
    else
        tglfmod = loadmodelonce(model_filename)
    end
    xnames = [replace(name, "_log10" => "") for name in tglfmod.xnames]
    x = collect(transpose(reduce(hcat, [Float64.(data[name]) for name in xnames])))
    y = tglfmod(x; uncertain, warn_nn_train_bounds, fidelity=:TGLFNN)
    if fidelity == :GKNN
        supported_gknn_models = ["sat3_em_d3d_azf-1", "sat3_em_d3d+mastu+nstx_azf-1", "sat3_em_d3d_azf-1_withnegD", "sat3_em_d3d_azf-1_gkdb", "sat2_em_d3d+mastu+nstx_azf-1", "sat3_em_d3d+mastu_azf-1"]
        if !(model_filename in supported_gknn_models)
            error("GKNN fidelity is not supported for model '$model_filename'. Supported models are: $(join(supported_gknn_models, ", "))")
        end
        if model_filename in ["sat3_em_d3d_azf-1"]
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
        elseif model_filename in ["sat3_em_d3d+mastu+nstx_azf-1", "sat3_em_d3d_azf-1_withnegD", "sat3_em_d3d_azf-1_gkdb", "sat2_em_d3d+mastu+nstx_azf-1"]
            gknn = loadmodelonce(model_filename * "_gknn31")
            err = gknn(vcat(x, y)...; uncertain, warn_nn_train_bounds, fidelity)
            y .*= err
            if model_filename in ["sat3_em_d3d_azf-1_gkdb"]
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