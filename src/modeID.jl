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

"""
    load_modeid_model(filename::AbstractString) -> ModeIDmodel | ModeIDensemble

Load a ModeID turbulence-mode classifier from a `.bson` file. `filename` is
resolved against the registered model search paths and the built-in `models/`
directory (see [`run_modeid_nn`](@ref)). Returns a `ModeIDensemble` when the
saved dictionary uses integer keys (one entry per ensemble member), otherwise a
single `ModeIDmodel`.
"""
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

    output_modes = [_YNAME_TO_MODE[yn] for yn in modeid_model.ynames]

    # Convert output probabilities to NNModeIdentification structs
    results = Vector{NNModeIdentification{Float64}}(undef, N)
    for i in 1:N
        probs = Dict{TurbulenceMode,Float64}()
        max_prob = -Inf
        dominant = ITG
        for (k, mode) in enumerate(output_modes)
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

#= ====================================== =#
#  QLNN-based Mode Identification
#= ====================================== =#

"""
    run_modeid_qlnn(input_tjlfs::Vector{InputTJLF{T}}; kw...) -> Vector{TJLFModeIdentification{T}}

Run the QLNN bundle to predict QL weights and eigenvalues, then classify turbulence
modes using the same QL-weight-ratio and frequency-sign logic as `identify_modes`
(the TJLF path).

Returns `TJLFModeIdentification` structs with every field the TJLF path produces:
`mode_per_ky`, `energy_flux_per_mode`, `dominant_mode`, `dominant_mode_fraction`,
`ky_spectrum`, and `flux_solution`.

# Keywords
- `bundle_name::AbstractString="QLNN"`: QLNN bundle directory name
- `warn_nn_train_bounds::Bool=false`: warn when NN inputs exceed training bounds
- `stability_threshold::Real=0.5`: P(unstable) threshold for hard gate
- `em_threshold`, `ion_electron_threshold`, `ky_etg`: mode classification parameters
  (same semantics as `identify_modes`)
"""
function run_modeid_qlnn(input_tjlfs::Vector{InputTJLF{T}};
                         bundle_name::AbstractString="QLNN",
                         warn_nn_train_bounds::Bool=false,
                         stability_threshold::Real=0.5,
                         em_threshold::Real=0.5,
                         ion_electron_threshold::Real=0.5,
                         ky_etg::Real=2.0) where {T<:Real}
    nr = length(input_tjlfs)
    nr == 0 && return TJLFModeIdentification{T}[]

    bundle = loadqlnnbundleonce(String(bundle_name))
    pred = _run_qlnn_predict(input_tjlfs, bundle; warn_nn_train_bounds)
    nf = pred.info_e.nf
    ns = pred.info_e.ns

    results = Vector{TJLFModeIdentification{T}}(undef, nr)
    Threads.@threads for r in 1:nr
        c0 = pred.chunk_starts[r]
        nk = pred.nky_r[r]
        cols = c0:c0+nk-1
        ks = pred.ky_spectrums[r]
        mask = nothing
        if pred.P_unstable !== nothing
            mask = Bool[pred.P_unstable[c0+j-1] >= stability_threshold for j in 1:nk]
        end
        results[r] = _qlnn_classify_modeid_radial(
            input_tjlfs[r], pred.sat_params_v[r], ks,
            view(pred.Y_energy,   :, cols),
            view(pred.Y_particle, :, cols),
            view(pred.Y_momentum, :, cols),
            view(pred.Y_eig,      :, cols),
            mask, pred.info_e, pred.info_p, pred.info_m,
            pred.eig_norm_by_ky, nf, ns;
            em_threshold, ion_electron_threshold, ky_etg
        )
    end
    return results
end

"""
    run_modeid_qlnn(dd::IMAS.dd, rho_transport; kw...) -> Vector{TJLFModeIdentification}

Convenience: build `InputTJLF` from IMAS `dd` and rho grid, then classify modes via QLNN.
Accepts the same `sat_rule`, `electromagnetic`, `lump_ions`, `MXH_modes` keywords as
`ActorModeID` to control the TGLF input construction.
"""
function run_modeid_qlnn(dd::IMAS.dd, rho_transport::AbstractVector{<:Real};
                         bundle_name::AbstractString="QLNN",
                         sat_rule::Symbol=:sat3,
                         electromagnetic::Bool=true,
                         lump_ions::Bool=true,
                         MXH_modes::Int=1,
                         warn_nn_train_bounds::Bool=false,
                         stability_threshold::Real=0.5,
                         em_threshold::Real=0.5,
                         ion_electron_threshold::Real=0.5,
                         ky_etg::Real=2.0)
    input_tglfs = InputTGLF(dd, rho_transport, sat_rule, electromagnetic, lump_ions; MXH_modes)
    if hasproperty(input_tglfs, :tglfs)
        input_tglfs = input_tglfs.tglfs
    end
    input_tjlfs = InputTJLF{Float64}[InputTJLF{Float64}(input_tglfs[k]) for k in eachindex(rho_transport)]
    return run_modeid_qlnn(input_tjlfs;
        bundle_name, warn_nn_train_bounds, stability_threshold,
        em_threshold, ion_electron_threshold, ky_etg)
end

# Per-radial-point helper: pack QLNN predictions into the QL tensor, classify
# modes from QL-weight ratios and eigenvalue frequency, run the TJLF saturation
# rule, and return a full TJLFModeIdentification.
function _qlnn_classify_modeid_radial(
    input_tjlf::TJLF.InputTJLF{T},
    sat_params,
    ky_spectrum::AbstractVector,
    y_energy::AbstractMatrix,
    y_particle::AbstractMatrix,
    y_momentum::AbstractMatrix,
    y_eig::AbstractMatrix,
    mask::Union{Nothing,AbstractVector{Bool}},
    info_e, info_p, info_m,
    eig_norm_by_ky::Bool,
    nf::Int, ns::Int;
    em_threshold::Real=0.5,
    ion_electron_threshold::Real=0.5,
    ky_etg::Real=2.0
) where {T<:Real}
    nky = length(ky_spectrum)
    @assert nf >= 2 "QLNN mode identification requires phi and Apar fields (nf=$nf)"
    @assert ns >= 2 "QLNN mode identification requires ≥2 species (ns=$ns)"

    # Pack QL tensor
    QL = zeros(T, nf, ns, 1, nky, 5)
    _qlnn_pack_qlweight!(QL, y_energy,   info_e, _QLNN_TARGET_TYPE_IDX[:energy],   ky_spectrum, mask)
    _qlnn_pack_qlweight!(QL, y_particle, info_p, _QLNN_TARGET_TYPE_IDX[:particle], ky_spectrum, mask)
    _qlnn_pack_qlweight!(QL, y_momentum, info_m, _QLNN_TARGET_TYPE_IDX[:momentum], ky_spectrum, mask)

    # Build Γ matrix and recover per-ky frequency (ω)
    Γ = zeros(T, 1, nky)
    frequencies = Vector{T}(undef, nky)
    Tz = zero(T)
    has_omega = size(y_eig, 1) >= 2
    @inbounds for k in 1:nky
        γ = y_eig[1, k]
        if eig_norm_by_ky
            γ *= T(ky_spectrum[k])
        end
        if mask !== nothing && !mask[k]
            γ = Tz
        end
        Γ[1, k] = γ
        if has_omega
            ω = y_eig[2, k]
            if eig_norm_by_ky
                ω *= T(ky_spectrum[k])
            end
            frequencies[k] = ω
        else
            frequencies[k] = Tz
        end
    end

    # --- Mode classification (same logic as identify_modes in tjlf.jl) ---
    ky_vec = collect(T, ky_spectrum)
    etg_indices = findall(>(ky_etg), ky_vec)
    signetg = zero(T)
    ignore_sign = true
    if !isempty(etg_indices)
        etg_avg = sum(frequencies[etg_indices]) / length(etg_indices)
        if isfinite(etg_avg) && etg_avg != 0
            signetg = sign(etg_avg)
            ignore_sign = false
        end
    end

    mode_per_ky = Vector{TurbulenceMode}(undef, nky)
    for j in 1:nky
        freq = frequencies[j]
        ky_val = ky_vec[j]
        qlw_es_e = abs(QL[1, 1, 1, j, 2])
        qlw_em_e = abs(QL[2, 1, 1, j, 2])
        qlw_es_i = abs(QL[1, 2, 1, j, 2])
        qlw_em_i = abs(QL[2, 2, 1, j, 2])

        em_ratio = qlw_es_e > eps(T) ? qlw_em_e / qlw_es_e : (qlw_em_e > eps(T) ? T(Inf) : zero(T))
        ie_ratio = qlw_es_e > eps(T) ? qlw_es_i / qlw_es_e : zero(T)

        if !ignore_sign
            mode_per_ky[j] = _classify_mode_with_sign(em_ratio, freq, ky_val, ie_ratio, signetg, em_threshold, ion_electron_threshold, ky_etg)
        else
            emi_ratio = qlw_em_e > eps(T) ? qlw_em_i / qlw_em_e : zero(T)
            mode_per_ky[j] = _classify_mode_no_sign(em_ratio, ky_val, ie_ratio, emi_ratio, em_threshold, ion_electron_threshold, ky_etg)
        end
    end

    # --- Flux integration via TJLF saturation rule ---
    QL_flux_out, flux_spectrum = TJLF.sum_ky_spectrum(input_tjlf, sat_params, Γ, QL)

    energy_flux_ky = vec(sum(@view(flux_spectrum[:, :, :, :, 2]); dims=(1, 2, 3)))

    energy_flux_per_mode = Dict{TurbulenceMode,T}(mode => zero(T) for mode in instances(TurbulenceMode))
    for j in 1:nky
        energy_flux_per_mode[mode_per_ky[j]] += energy_flux_ky[j]
    end

    dominant = first(instances(TurbulenceMode))
    max_flux = typemin(T)
    for mode in instances(TurbulenceMode)
        if energy_flux_per_mode[mode] > max_flux
            max_flux = energy_flux_per_mode[mode]
            dominant = mode
        end
    end

    total_abs = sum(abs, values(energy_flux_per_mode))
    dominant_fraction = total_abs > 0 ? abs(energy_flux_per_mode[dominant]) / total_abs : one(T)

    flux_sol = GACODE.FluxSolution{T}(
        TJLF.Qe(QL_flux_out), TJLF.Qi(QL_flux_out),
        TJLF.Γe(QL_flux_out), TJLF.Γi(QL_flux_out), TJLF.Πi(QL_flux_out))

    return TJLFModeIdentification{T}(mode_per_ky, energy_flux_per_mode, dominant, dominant_fraction, ky_vec, flux_sol)
end

export NNModeIdentification, run_modeid_nn, load_modeid_model, run_modeid_qlnn
