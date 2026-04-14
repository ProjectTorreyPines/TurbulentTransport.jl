"""
    TurbulenceMode

Plasma turbulence mode types identified from quasilinear transport analysis.

Colors follow the convention: ITG=green, TEM=orange, KBM=violet, ETG=blue, MTM=red.
"""
@enum TurbulenceMode ITG TEM KBM ETG MTM

const MODE_COLORS = Dict(ITG => :green, TEM => :orange, KBM => :violet, ETG => :blue, MTM => :red)
const MODE_LABELS = Dict(ITG => "ITG", TEM => "TEM", KBM => "KBM", ETG => "ETG", MTM => "MTM")

"""
    TJLFModeIdentification{T<:Real}

Results of turbulence mode identification from a TJLF run at a single radial location.

# Fields
- `mode_per_ky`: classified mode at each ky (based on the most unstable linear mode)
- `energy_flux_per_mode`: total energy flux contribution by mode type (summed over ky points classified as that mode)
- `dominant_mode`: mode type with the largest energy flux contribution
- `dominant_mode_fraction`: fraction of total |energy flux| from the dominant mode
- `ky_spectrum`: ky values used in the TJLF run
- `flux_solution`: total quasilinear fluxes from the TJLF run
"""
struct TJLFModeIdentification{T<:Real}
    mode_per_ky::Vector{TurbulenceMode}
    energy_flux_per_mode::Dict{TurbulenceMode,T}
    dominant_mode::TurbulenceMode
    dominant_mode_fraction::T
    ky_spectrum::Vector{T}
    flux_solution::GACODE.FluxSolution{T}
end

function Base.show(io::IO, ::MIME"text/plain", mid::TJLFModeIdentification)
    println(io, "TJLFModeIdentification:")
    println(io, "  Dominant mode: $(mid.dominant_mode) ($(round(mid.dominant_mode_fraction * 100; digits=1))% of total |flux|)")
    println(io, "  Energy flux by mode:")
    for mode in instances(TurbulenceMode)
        flux = mid.energy_flux_per_mode[mode]
        flux == 0 && continue
        println(io, "    $(MODE_LABELS[mode]): $(round(flux; sigdigits=4))")
    end
end

function _classify_mode_with_sign(em_ratio, freq, ky_val, ie_ratio, signetg, em_threshold, ion_electron_threshold, ky_etg)
    if em_ratio > em_threshold
        return freq * signetg >= 0 ? MTM : KBM
    elseif freq * signetg > 0
        if ie_ratio >= ion_electron_threshold
            return TEM
        else
            return ky_val > ky_etg ? ETG : TEM
        end
    else
        return ITG
    end
end

function _classify_mode_no_sign(em_ratio, ky_val, ie_ratio, emi_ratio, em_threshold, ion_electron_threshold, ky_etg)
    if em_ratio > em_threshold
        return emi_ratio > 0.25 ? KBM : MTM
    else
        if ie_ratio >= 1.5
            return ITG
        elseif ie_ratio >= ion_electron_threshold
            return TEM
        else
            return ky_val > ky_etg ? ETG : TEM
        end
    end
end

"""
    identify_modes(tjlf_result::NamedTuple, input_tjlf::InputTJLF; kw...)

Classify turbulence modes from pre-computed TJLF results and determine the dominant
mode driving heat flux.

Mode classification at each ky uses the most unstable linear mode's quasilinear weights
and frequency to assign one of: ITG, TEM, KBM, ETG, MTM. The saturated flux spectrum
is then partitioned by mode type to determine which mode drives the most transport.

# Classification rules (when ETG-scale frequency sign is available)
- **MTM**: EM-dominated (`|QL_Apar_e/QL_phi_e| > em_threshold`) and electron-direction frequency
- **KBM**: EM-dominated and ion-direction frequency
- **TEM**: ES-dominated, electron-direction, and either significant ion QL weight or low ky
- **ETG**: ES-dominated, electron-direction, weak ion QL weight, and high ky
- **ITG**: ES-dominated and ion-direction frequency

# Arguments
- `tjlf_result`: output of `TJLF.run(input_tjlf)` (NamedTuple with QL_weights, eigenvalue, flux_spectrum, QL_flux_out)
- `input_tjlf`: the InputTJLF used for the run

# Keywords
- `em_threshold::Real=0.5`: EM/ES QL weight ratio threshold separating electromagnetic (MTM/KBM) from electrostatic (ITG/TEM/ETG) modes
- `ion_electron_threshold::Real=0.5`: ion/electron ES QL weight ratio threshold for TEM vs ETG classification
- `ky_etg::Real=2.0`: ky threshold above which ES electron-direction modes are classified as ETG
"""
function identify_modes(
    tjlf_result::NamedTuple,
    input_tjlf::InputTJLF{T};
    em_threshold::Real=0.5,
    ion_electron_threshold::Real=0.5,
    ky_etg::Real=2.0
) where {T<:Real}
    QL_weights = tjlf_result.QL_weights          # (nf, ns, nm, nky, ntype)
    eigenvalue = tjlf_result.eigenvalue           # (nm, nky, 2): [:,:,1]=gamma, [:,:,2]=freq
    flux_spectrum = tjlf_result.flux_spectrum     # (nf, ns, nm, nky, ntype)
    QL_flux_out = tjlf_result.QL_flux_out         # (nf, ns, ntype) integrated

    @assert !isempty(QL_flux_out) "Mode identification requires USE_TRANSPORT_MODEL=true"

    ky_spectrum = collect(T, input_tjlf.KY_SPECTRUM)
    nky = length(ky_spectrum)
    ns = size(QL_weights, 2)
    @assert ns >= 2 "Need at least 2 species (electrons + ion) for mode identification"

    # Determine frequency sign convention: average frequency at ETG-scale ky
    frequencies = eigenvalue[1, :, 2]
    etg_indices = findall(>(ky_etg), ky_spectrum)
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
        freq = eigenvalue[1, j, 2]
        ky_val = ky_spectrum[j]

        # QL weights: field 1=phi, 2=A∥; species 1=electron, 2=first ion; mode 1; type 2=energy
        qlw_es_e = abs(QL_weights[1, 1, 1, j, 2])
        qlw_em_e = abs(QL_weights[2, 1, 1, j, 2])
        qlw_es_i = abs(QL_weights[1, 2, 1, j, 2])
        qlw_em_i = abs(QL_weights[2, 2, 1, j, 2])

        em_ratio = qlw_es_e > eps(T) ? qlw_em_e / qlw_es_e : (qlw_em_e > eps(T) ? T(Inf) : zero(T))
        ie_ratio = qlw_es_e > eps(T) ? qlw_es_i / qlw_es_e : zero(T)

        if !ignore_sign
            mode_per_ky[j] = _classify_mode_with_sign(em_ratio, freq, ky_val, ie_ratio, signetg, em_threshold, ion_electron_threshold, ky_etg)
        else
            emi_ratio = qlw_em_e > eps(T) ? qlw_em_i / qlw_em_e : zero(T)
            mode_per_ky[j] = _classify_mode_no_sign(em_ratio, ky_val, ie_ratio, emi_ratio, em_threshold, ion_electron_threshold, ky_etg)
        end
    end

    # Energy flux at each ky: sum flux_spectrum over fields (1), species (2), modes (3) for energy channel (type=2)
    energy_flux_ky = vec(sum(@view(flux_spectrum[:, :, :, :, 2]); dims=(1, 2, 3)))

    energy_flux_per_mode = Dict{TurbulenceMode,T}(mode => zero(T) for mode in instances(TurbulenceMode))
    for j in 1:nky
        energy_flux_per_mode[mode_per_ky[j]] += energy_flux_ky[j]
    end

    # Dominant mode: largest signed energy flux (most outward transport), matching Python convention
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

    return TJLFModeIdentification{T}(mode_per_ky, energy_flux_per_mode, dominant, dominant_fraction, ky_spectrum, flux_sol)
end

"""
    identify_modes(input_tjlf::InputTJLF; kw...)

Run TJLF and classify turbulence modes. Returns a `TJLFModeIdentification` containing
the mode classification at each ky, energy flux breakdown by mode, and the dominant mode.
"""
function identify_modes(input_tjlf::InputTJLF{T}; kw...) where {T<:Real}
    tjlf_result = TJLF.run(input_tjlf)
    return identify_modes(tjlf_result, input_tjlf; kw...)
end

"""
    identify_modes(input_tglf::InputTGLF; kw...)

Convert InputTGLF to InputTJLF, run TJLF, and classify turbulence modes.
"""
function identify_modes(input_tglf::InputTGLF; kw...)
    input_tjlf = InputTJLF{Float64}(input_tglf)
    return identify_modes(input_tjlf; kw...)
end

"""
    identify_modes(input_tjlfs::Vector{InputTJLF{T}}; kw...) -> Vector{TJLFModeIdentification{T}}

Run TJLF and classify modes at multiple radial locations (threaded).
"""
function identify_modes(input_tjlfs::Vector{InputTJLF{T}}; kw...) where {T<:Real}
    results = Vector{TJLFModeIdentification{T}}(undef, length(input_tjlfs))
    Threads.@threads for idx in eachindex(input_tjlfs)
        results[idx] = identify_modes(input_tjlfs[idx]; kw...)
    end
    return results
end

export TurbulenceMode, ITG, TEM, KBM, ETG, MTM
export TJLFModeIdentification, identify_modes
export MODE_COLORS, MODE_LABELS

# ========== Original run_tjlf functions ==========

function run_tjlf(input_tjlf::InputTJLF{T}) where {T<:Real}
    QL_flux_out = TJLF.run_tjlf(input_tjlf)
    return GACODE.FluxSolution{T}(TJLF.Qe(QL_flux_out), TJLF.Qi(QL_flux_out), TJLF.Γe(QL_flux_out), TJLF.Γi(QL_flux_out), TJLF.Πi(QL_flux_out))
end

function run_tjlf(input_tglf::InputTGLF)
    input_tjlf = InputTJLF{Float64}(input_tglf)
    return run_tjlf(input_tjlf)
end
