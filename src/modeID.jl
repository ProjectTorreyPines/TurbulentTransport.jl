import Flux
import BSON
import Dates
import Memoize

#= ====================================== =#
#  ModeID Model Types and Loading
#= ====================================== =#

struct ModeIDmodel
    fluxmodel::Flux.Chain
    name::String
    date::Dates.DateTime
    xnames::Vector{String}
    ynames::Vector{String}  # Mode names: ["MTM", "KBM", "TEM", "ITG", "ETG"]
    xm::Vector{Float32}
    xσ::Vector{Float32}
    xbounds::Array{Float32}
end

struct ModeIDensemble
    models::Vector{ModeIDmodel}
end

function Base.show(io::IO, ::MIME"text/plain", model::ModeIDmodel)
    println(io, "ModeIDmodel")
    println(io, "name: $(model.name)")
    println(io, "date: $(model.date)")
    println(io, "xnames ($(length(model.xnames))): $(model.xnames)")
    return println(io, "ynames ($(length(model.ynames))): $(model.ynames)")
end

function Base.show(io::IO, ::MIME"text/plain", ens::ModeIDensemble)
    println(io, "ModeIDensemble ($(length(ens.models)) models)")
    return show(io, MIME"text/plain"(), ens.models[1])
end

function Base.getproperty(ensemble::ModeIDensemble, field::Symbol)
    if field == :models
        return getfield(ensemble, field)
    else
        return getfield(ensemble.models[1], field)
    end
end

function _dict2modeid(dict::AbstractDict)
    fluxmodel = Flux.fmap(Flux.f32, dict[:fluxmodel])
    return ModeIDmodel(
        fluxmodel,
        String(dict[:name]),
        dict[:date],
        String.(dict[:xnames]),
        String.(dict[:ynames]),
        Float32.(vec(dict[:xm])),
        Float32.(vec(dict[:xσ])),
        Float32.(dict[:xbounds])
    )
end

function _dict2modeid_ensemble(dict::Dict)
    return ModeIDensemble([_dict2modeid(modict) for modict in values(dict)])
end

function load_modeid_model(filename::AbstractString)
    fullpath = resolve_model_path(filename; extensions=[".bson"])
    savedict = BSON.load(fullpath, @__MODULE__)
    if typeof(first(keys(savedict))) <: Integer
        return _dict2modeid_ensemble(savedict)
    else
        return _dict2modeid(savedict)
    end
end

Memoize.@memoize function load_modeid_model_once(filename::String)
    return load_modeid_model(filename)
end

#= ====================================== =#
#  ModeID Prediction
#= ====================================== =#

function predict_modeid(model::ModeIDmodel, x::AbstractVector{T}) where {T<:Real}
    xn = (Float32.(x) .- model.xm) ./ model.xσ
    logits = model.fluxmodel(xn)
    return Flux.softmax(logits)
end

function predict_modeid(model::ModeIDmodel, x::AbstractMatrix{T}) where {T<:Real}
    xn = (Float32.(x) .- model.xm) ./ model.xσ
    logits = model.fluxmodel(xn)
    return Flux.softmax(logits)
end

function predict_modeid(ensemble::ModeIDensemble, x::AbstractVecOrMat{T}) where {T<:Real}
    preds = [predict_modeid(m, x) for m in ensemble.models]
    return sum(preds) ./ length(preds)
end

#= ====================================== =#
#  Log-transformed feature helpers
#= ====================================== =#

const _LOG10_FEATURES = Dict(
    "BETAE_log10" => "BETAE",
    "DEBYE_log10" => "DEBYE",
    "XNUE_log10" => "XNUE"
)

#= ====================================== =#
#  ModeID NN Results
#= ====================================== =#

const _YNAME_TO_MODE = Dict(
    "MTM" => MTM, "OUT_MODE_MTM" => MTM,
    "KBM" => KBM, "OUT_MODE_KBM" => KBM,
    "TEM" => TEM, "OUT_MODE_TEM" => TEM,
    "ITG" => ITG, "OUT_MODE_ITG" => ITG,
    "ETG" => ETG, "OUT_MODE_ETG" => ETG
)

"""
    NNModeIdentification{T<:Real}

Results of neural-network turbulence mode identification at a single radial location.

# Fields
- `probabilities`: Dict mapping each `TurbulenceMode` to its predicted probability
- `dominant_mode`: mode with the highest predicted probability
- `dominant_mode_fraction`: probability of the dominant mode (for consistency with `TJLFModeIdentification`)
"""
struct NNModeIdentification{T<:Real} <: AbstractModeIdentification{T}
    probabilities::Dict{TurbulenceMode,T}
    dominant_mode::TurbulenceMode
    dominant_mode_fraction::T
end

function Base.show(io::IO, ::MIME"text/plain", mid::NNModeIdentification)
    println(io, "NNModeIdentification:")
    println(io, "  Dominant mode: $(mid.dominant_mode) ($(round(mid.dominant_mode_fraction * 100; digits=1))%)")
    println(io, "  Probabilities:")
    for mode in instances(TurbulenceMode)
        prob = mid.probabilities[mode]
        prob < 0.001 && continue
        println(io, "    $(MODE_LABELS[mode]): $(round(prob * 100; digits=1))%")
    end
end

#= ====================================== =#
#  ModeID NN inference from IMAS data
#= ====================================== =#

"""
    run_modeid_nn(dd::IMAS.dd, rho_transport::AbstractVector{<:Real};
                  model_filename::String, warn_nn_train_bounds::Bool=false, MXH_modes::Int=1)

Run ModeID neural network to predict the dominant turbulence mode at each radial location.

Instead of running TJLF quasilinear analysis (which takes seconds per radial point),
the ModeID NN directly classifies the dominant mode from TGLF input parameters in a
single forward pass (< 1 ms for all radial points).

The model takes 34 TGLF input parameters (geometry, gradients, collisionality, etc.)
and outputs a 5-class softmax probability vector: [MTM, KBM, TEM, ITG, ETG].

Returns a `Vector{NNModeIdentification}` with one result per rho point.
"""
function run_modeid_nn(dd::IMAS.dd, rho_transport::AbstractVector{<:Real};
                       model_filename::String,
                       warn_nn_train_bounds::Bool=false,
                       MXH_modes::Int=1)

    modeid_model = load_modeid_model_once(model_filename)

    input_tglfs = InputTGLF(dd, rho_transport, :sat3, true, true; MXH_modes)
    if hasproperty(input_tglfs, :tglfs)
        input_tglfs = input_tglfs.tglfs
    end

    N = length(rho_transport)
    n_features = length(modeid_model.xnames)
    inputs = Matrix{Float32}(undef, n_features, N)

    for (i, rho) in enumerate(rho_transport)
        it = input_tglfs[i]
        for (j, xname) in enumerate(modeid_model.xnames)
            if haskey(_LOG10_FEATURES, xname)
                # Log10-transformed feature: get base parameter and apply log10
                base_name = _LOG10_FEATURES[xname]
                val = getfield(it, Symbol(base_name))
                if ismissing(val)
                    error("ModeID input field '$base_name' is Missing at rho=$rho")
                end
                inputs[j, i] = Float32(log10(max(Float64(val), 1e-10)))
            else
                val = getfield(it, Symbol(xname))
                if ismissing(val)
                    error("ModeID input field '$xname' is Missing at rho=$rho")
                end
                inputs[j, i] = Float32(val)
            end
        end
    end

    if warn_nn_train_bounds
        for j in 1:n_features
            for i in 1:N
                if inputs[j, i] < modeid_model.xbounds[j, 1]
                    @warn "ModeID: $(modeid_model.xnames[j])=$(inputs[j,i]) below training bound $(modeid_model.xbounds[j,1]) at rho=$(rho_transport[i])"
                elseif inputs[j, i] > modeid_model.xbounds[j, 2]
                    @warn "ModeID: $(modeid_model.xnames[j])=$(inputs[j,i]) above training bound $(modeid_model.xbounds[j,2]) at rho=$(rho_transport[i])"
                end
            end
        end
    end

    output = predict_modeid(modeid_model, inputs)

    # Convert output probabilities to NNModeIdentification structs
    results = Vector{NNModeIdentification{Float64}}(undef, N)
    for i in 1:N
        probs = Dict{TurbulenceMode,Float64}()
        max_prob = -Inf
        dominant = ITG
        for (k, yname) in enumerate(modeid_model.ynames)
            mode = _YNAME_TO_MODE[yname]
            p = Float64(output[k, i])
            probs[mode] = p
            if p > max_prob
                max_prob = p
                dominant = mode
            end
        end
        results[i] = NNModeIdentification{Float64}(probs, dominant, max_prob)
    end

    return results
end

export NNModeIdentification, run_modeid_nn, load_modeid_model
