import Flux
import BSON
import Dates
import Memoize

#= ====================================== =#
#  FINN Model Types and Loading
#= ====================================== =#

struct FINNmodel
    fluxmodel::Flux.Chain
    name::String
    date::Dates.DateTime
    xnames::Vector{String}
    ynames::Vector{String}
    xm::Vector{Float64}
    xσ::Vector{Float64}
    xbounds::Array{Float64}
    ym::Vector{Float64}
    yσ::Vector{Float64}
    ybounds::Array{Float64}
end

struct FINNensemble
    models::Vector{FINNmodel}
end

function Base.show(io::IO, ::MIME"text/plain", model::FINNmodel)
    println(io, "FINNmodel")
    println(io, "name: $(model.name)")
    println(io, "date: $(model.date)")
    println(io, "xnames ($(length(model.xnames))): $(model.xnames)")
    return println(io, "ynames ($(length(model.ynames))): $(model.ynames)")
end

function Base.show(io::IO, ::MIME"text/plain", ens::FINNensemble)
    println(io, "FINNensemble ($(length(ens.models)) models)")
    return show(io, MIME"text/plain"(), ens.models[1])
end

function Base.getproperty(ensemble::FINNensemble, field::Symbol)
    if field == :models
        return getfield(ensemble, field)
    else
        return getfield(ensemble.models[1], field)
    end
end

function _dict2finn(dict::AbstractDict)
    fluxmodel = Flux.fmap(Flux.f64, dict[:fluxmodel])
    return FINNmodel(
        fluxmodel,
        String(dict[:name]),
        dict[:date],
        String.(dict[:xnames]),
        String.(dict[:ynames]),
        Float64.(vec(dict[:xm])),
        Float64.(vec(dict[:xσ])),
        Float64.(dict[:xbounds]),
        Float64.(vec(dict[:ym])),
        Float64.(vec(dict[:yσ])),
        Float64.(dict[:ybounds])
    )
end

function _dict2finn_ensemble(dict::Dict)
    return FINNensemble([_dict2finn(modict) for modict in values(dict)])
end

"""
    load_finn_model(filename::AbstractString) -> FINNmodel | FINNensemble

Load a FINN model from a `.bson` file. `filename` is resolved against the
registered model search paths and the built-in `models/` directory (see
[`run_finn`](@ref)). Returns a `FINNensemble` when the saved dictionary uses
integer keys (one entry per ensemble member), otherwise a single `FINNmodel`.
"""
function load_finn_model(filename::AbstractString)
    fullpath = resolve_model_path(filename; extensions=[".bson"])
    savedict = BSON.load(fullpath, @__MODULE__)
    if typeof(first(keys(savedict))) <: Integer
        return _dict2finn_ensemble(savedict)
    else
        return _dict2finn(savedict)
    end
end

Memoize.@memoize function load_finn_model_once(filename::String)
    return load_finn_model(filename)
end

#= ====================================== =#
#  FINN Prediction
#= ====================================== =#

function predict_finn(model::FINNmodel, x::AbstractVector{T}) where {T<:Real}
    xn = (Float64.(x) .- model.xm) ./ model.xσ
    yn = model.fluxmodel(xn)
    return yn .* model.yσ .+ model.ym
end

function predict_finn(model::FINNmodel, x::AbstractMatrix{T}) where {T<:Real}
    xn = (Float64.(x) .- model.xm) ./ model.xσ
    yn = model.fluxmodel(xn)
    return yn .* model.yσ .+ model.ym
end

function predict_finn(ensemble::FINNensemble, x::AbstractVecOrMat{T}) where {T<:Real}
    preds = [predict_finn(m, x) for m in ensemble.models]
    return sum(preds) ./ length(preds)
end

#= ====================================== =#
#  FINN inference from IMAS data
#= ====================================== =#


"""
    run_finn(dd::IMAS.dd, rho_transport::AbstractVector{<:Real};
             model_filename::String, warn_nn_train_bounds::Bool=false, MXH_modes::Int=1)

Run FINN (Flux-matcher Inversion Neural Network) to directly predict
flux-matched TGLF gradients from geometry and sources.

Instead of iteratively solving for profiles (as in the flux matcher),
FINN directly predicts the converged gradients in a single forward pass (< 1 ms).

Returns a NamedTuple with:
  - `RLTS_1`: electron temperature gradient a/L_Te at each rho
  - `RLTS_2`: ion temperature gradient a/L_Ti at each rho
  - `RLNS_1`: electron density gradient a/L_ne at each rho
  - `VEXB_SHEAR`: ExB shearing rate at each rho
  - `rho`: rho_transport grid used
"""
function run_finn(dd::IMAS.dd, rho_transport::AbstractVector{<:Real};
                  model_filename::String,
                  warn_nn_train_bounds::Bool=false,
                  MXH_modes::Int=1)

    finn_model = load_finn_model_once(model_filename)

    input_tglfs = InputTGLF(dd, rho_transport, :sat3, true, true; MXH_modes)
    if hasproperty(input_tglfs, :tglfs)
        input_tglfs = input_tglfs.tglfs
    end

    sources_gB = GACODE.sources_to_gyrobohm(dd, rho_transport)

    N = length(rho_transport)
    n_features = length(finn_model.xnames)
    inputs = Matrix{Float64}(undef, n_features, N)

    for (i, rho) in enumerate(rho_transport)
        it = input_tglfs[i]
        for (j, xname) in enumerate(finn_model.xnames)
            if xname == "rho"
                inputs[j, i] = rho
            elseif xname == "Qe"
                inputs[j, i] = sources_gB.Qe[i]
            elseif xname == "Qi"
                inputs[j, i] = sources_gB.Qi[i]
            elseif xname == "Ge"
                inputs[j, i] = sources_gB.Ge[i]
            elseif xname == "Pi"
                inputs[j, i] = sources_gB.Pi[i]
            else
                val = getfield(it, Symbol(xname))
                if ismissing(val)
                    error("FINN input field '$xname' is Missing at rho=$rho")
                end
                inputs[j, i] = Float64(val)
            end
        end
    end

    if warn_nn_train_bounds
        for j in 1:n_features
            for i in 1:N
                if inputs[j, i] < finn_model.xbounds[j, 1]
                    @warn "FINN: $(finn_model.xnames[j])=$(inputs[j,i]) below training bound $(finn_model.xbounds[j,1]) at rho=$(rho_transport[i])"
                elseif inputs[j, i] > finn_model.xbounds[j, 2]
                    @warn "FINN: $(finn_model.xnames[j])=$(inputs[j,i]) above training bound $(finn_model.xbounds[j,2]) at rho=$(rho_transport[i])"
                end
            end
        end
    end

    output = predict_finn(finn_model, inputs)

    ynames_clean = [replace(yn, "OUT_" => "") for yn in finn_model.ynames]
    result = Dict{String,Vector{Float64}}()
    for (k, name) in enumerate(ynames_clean)
        result[name] = output[k, :]
    end

    return (
        RLTS_1=get(result, "RLTS_1", zeros(N)),
        RLTS_2=get(result, "RLTS_2", zeros(N)),
        RLNS_1=get(result, "RLNS_1", zeros(N)),
        VEXB_SHEAR=get(result, "VEXB_SHEAR", zeros(N)),
        rho=collect(rho_transport)
    )
end

export run_finn, load_finn_model
