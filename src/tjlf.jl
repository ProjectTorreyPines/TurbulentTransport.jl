"""
    TurbulenceMode

Plasma turbulence mode types identified from quasilinear transport analysis.

Colors follow the convention: ITG=green, TEM=orange, KBM=violet, ETG=blue, MTM=red.
"""
@enum TurbulenceMode ITG TEM KBM ETG MTM

@doc "Ion Temperature Gradient (ITG) mode: ion-scale, electrostatic, ion-direction frequency." ITG
@doc "Trapped Electron Mode (TEM): ion-scale, electrostatic, electron-direction frequency." TEM
@doc "Kinetic Ballooning Mode (KBM): electromagnetic, ion-direction frequency." KBM
@doc "Electron Temperature Gradient (ETG) mode: electron-scale (high ky), electrostatic." ETG
@doc "Micro-Tearing Mode (MTM): electromagnetic, electron-direction frequency." MTM

"""
    MODE_COLORS::Dict{TurbulenceMode,Symbol}

Plot color associated with each [`TurbulenceMode`](@ref):
ITG=green, TEM=orange, KBM=violet, ETG=blue, MTM=red.
"""
const MODE_COLORS = Dict(ITG => :green, TEM => :orange, KBM => :violet, ETG => :blue, MTM => :red)

"""
    MODE_LABELS::Dict{TurbulenceMode,String}

Short string label (`"ITG"`, `"TEM"`, ...) for each [`TurbulenceMode`](@ref),
used in plots and printed summaries.
"""
const MODE_LABELS = Dict(ITG => "ITG", TEM => "TEM", KBM => "KBM", ETG => "ETG", MTM => "MTM")

"""
    AbstractModeIdentification{T}

Abstract supertype for all mode identification results.
Subtypes must have fields `dominant_mode::TurbulenceMode` and `dominant_mode_fraction::T`.
"""
abstract type AbstractModeIdentification{T<:Real} end

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
struct TJLFModeIdentification{T<:Real} <: AbstractModeIdentification{T}
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
    classify_modes_core(ky_spectrum, freq_by_ky, qlw_es_e, qlw_em_e, qlw_es_i, qlw_em_i,
                        energy_flux_by_ky; em_threshold=0.5, ion_electron_threshold=0.5,
                        ky_etg=2.0, signetg_ref=nothing)
        -> NamedTuple (mode_per_ky, energy_flux_per_mode, dominant_mode,
                       dominant_fraction, signetg, ignore_sign)

Pure per-ky turbulence-mode classifier extracted from [`identify_modes`](@ref).
Operates directly on raw per-ky arrays so callers that already have the spectra in
memory (e.g. a TGLF-database parser) can classify without an `InputTJLF`/TJLF run.

All array arguments are length `nky`:
- `freq_by_ky`: eigenmode frequency at each ky (for the eigenmode being classified)
- `qlw_es_e/qlw_em_e/qlw_es_i/qlw_em_i`: QL energy weights (type=energy) for the same
  eigenmode — electrostatic/electromagnetic (phi/Apar) x electron/ion
- `energy_flux_by_ky`: (saturated) energy flux at each ky, partitioned by the mode
  classified there

The ETG-scale frequency sign is normally derived from `freq_by_ky` at `ky>ky_etg`.
Pass `signetg_ref` to reuse a sign computed from another eigenmode (e.g. use the
most-unstable eigenmode's sign when classifying the subdominant one); `signetg_ref=0`
forces the no-sign classification path. The returned `signetg`/`ignore_sign` can be
fed back as `signetg_ref` for a subsequent call.
"""
function classify_modes_core(
    ky_spectrum::AbstractVector,
    freq_by_ky::AbstractVector,
    qlw_es_e::AbstractVector,
    qlw_em_e::AbstractVector,
    qlw_es_i::AbstractVector,
    qlw_em_i::AbstractVector,
    energy_flux_by_ky::AbstractVector;
    em_threshold::Real=0.5,
    ion_electron_threshold::Real=0.5,
    ky_etg::Real=2.0,
    signetg_ref::Union{Nothing,Real}=nothing
)
    T = float(promote_type(eltype(ky_spectrum), eltype(freq_by_ky)))
    nky = length(ky_spectrum)

    # Determine frequency sign convention: average frequency at ETG-scale ky
    local signetg::T, ignore_sign::Bool
    if signetg_ref === nothing
        etg_indices = findall(>(ky_etg), ky_spectrum)
        signetg = zero(T)
        ignore_sign = true
        if !isempty(etg_indices)
            etg_avg = sum(@view freq_by_ky[etg_indices]) / length(etg_indices)
            if isfinite(etg_avg) && etg_avg != 0
                signetg = T(sign(etg_avg))
                ignore_sign = false
            end
        end
    else
        signetg = T(signetg_ref)
        ignore_sign = signetg == 0
    end

    mode_per_ky = Vector{TurbulenceMode}(undef, nky)
    for j in 1:nky
        freq = freq_by_ky[j]
        ky_val = ky_spectrum[j]

        es_e = abs(qlw_es_e[j])
        em_e = abs(qlw_em_e[j])
        es_i = abs(qlw_es_i[j])
        em_i = abs(qlw_em_i[j])

        em_ratio = es_e > eps(T) ? em_e / es_e : (em_e > eps(T) ? T(Inf) : zero(T))
        ie_ratio = es_e > eps(T) ? es_i / es_e : zero(T)

        if !ignore_sign
            mode_per_ky[j] = _classify_mode_with_sign(em_ratio, freq, ky_val, ie_ratio, signetg, em_threshold, ion_electron_threshold, ky_etg)
        else
            emi_ratio = em_e > eps(T) ? em_i / em_e : zero(T)
            mode_per_ky[j] = _classify_mode_no_sign(em_ratio, ky_val, ie_ratio, emi_ratio, em_threshold, ion_electron_threshold, ky_etg)
        end
    end

    energy_flux_per_mode = Dict{TurbulenceMode,T}(mode => zero(T) for mode in instances(TurbulenceMode))
    for j in 1:nky
        energy_flux_per_mode[mode_per_ky[j]] += T(energy_flux_by_ky[j])
    end

    # Dominant mode: largest signed energy flux (most outward transport)
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

    return (mode_per_ky=mode_per_ky, energy_flux_per_mode=energy_flux_per_mode,
        dominant_mode=dominant, dominant_fraction=dominant_fraction,
        signetg=signetg, ignore_sign=ignore_sign)
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

    # Per-ky quantities for the most-unstable eigenmode (mode index 1).
    # QL weights: field 1=phi, 2=A∥; species 1=electron, 2=first ion; type 2=energy
    freq_by_ky = T[eigenvalue[1, j, 2] for j in 1:nky]
    qlw_es_e = T[abs(QL_weights[1, 1, 1, j, 2]) for j in 1:nky]
    qlw_em_e = T[abs(QL_weights[2, 1, 1, j, 2]) for j in 1:nky]
    qlw_es_i = T[abs(QL_weights[1, 2, 1, j, 2]) for j in 1:nky]
    qlw_em_i = T[abs(QL_weights[2, 2, 1, j, 2]) for j in 1:nky]

    # Energy flux at each ky: sum flux_spectrum over fields (1), species (2), modes (3) for energy channel (type=2)
    energy_flux_ky = vec(sum(@view(flux_spectrum[:, :, :, :, 2]); dims=(1, 2, 3)))

    res = classify_modes_core(
        ky_spectrum, freq_by_ky, qlw_es_e, qlw_em_e, qlw_es_i, qlw_em_i, energy_flux_ky;
        em_threshold, ion_electron_threshold, ky_etg)

    flux_sol = GACODE.FluxSolution{T}(
        TJLF.Qe(QL_flux_out), TJLF.Qi(QL_flux_out),
        TJLF.Γe(QL_flux_out), TJLF.Γi(QL_flux_out), TJLF.Πi(QL_flux_out))

    return TJLFModeIdentification{T}(res.mode_per_ky, res.energy_flux_per_mode,
        res.dominant_mode, res.dominant_fraction, ky_spectrum, flux_sol)
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
export AbstractModeIdentification, TJLFModeIdentification, identify_modes, classify_modes_core
export MODE_COLORS, MODE_LABELS

# ========== Fluctuation spectra reconstruction ==========

"""
    TJLFFluctuationSpectra{T<:Real}

2D fluctuation spectrum reconstruction from a TJLF SAT saturation run.

The potential spectrum is reconstructed from the Staebler spectral-shift Lorentzian
kx distribution ([Staebler et al. PoP 2016, NF 2017, NF 2021]):

```
|φ|²(kx, ky, m) = phinorm(ky, m) · L(u(kx)) / L(u(0))
        L(u)  = 1 / [ (1 + ay · u²)² · (1 + |ax · u|^exp_ax)² ]
        u(kx) = (kx - kx0_e(ky)) · ky / kx_width(ky)
```

Normalization convention: `phi2[kx=0, ky, m] == phinorm(ky, m)` matches the TGLF
`phinorm` (the field intensity at lab-frame `kx = 0`). The Lorentzian peaks at
`kx = kx0_e`, where `|φ|² = phinorm / L(u(0))`.

Density fluctuation spectra for species `s` are given by
`|δnₛ|²(kx, ky, m) = N_weight[s, m, ky] · |φ|²(kx, ky, m)`, populated only when
species-resolved density QL weights are supplied.

# Fields
- `kx::Vector{T}`          — kx grid (ρ_s units)
- `ky::Vector{T}`          — ky grid (ρ_s units; copy of `input_tjlf.KY_SPECTRUM`)
- `phi2::Array{T,3}`       — `|φ|²(kx, ky, mode)`, shape `(nkx, nky, nmodes)`
- `phi2_ky::Matrix{T}`     — peak intensity `phinorm(ky, mode)`, shape `(nky, nmodes)`
- `kx_width::Vector{T}`    — SAT kx-Lorentzian width per ky, shape `(nky,)`
- `kx0_e::Vector{T}`       — spectral shift per ky, shape `(nky,)`
- `ax::T`, `ay::T`, `exp_ax::Int` — SAT Lorentzian coefficients
- `sat_rule::Int`          — SAT rule used
- `density2::Union{Nothing, Array{T,4}}` — `|δnₛ|²(kx, ky, mode, species)`,
   shape `(nkx, nky, nmodes, ns)`, or `nothing` if density weights were not supplied.
"""
struct TJLFFluctuationSpectra{T<:Real}
    kx::Vector{T}
    ky::Vector{T}
    phi2::Array{T,3}
    phi2_ky::Matrix{T}
    kx_width::Vector{T}
    kx0_e::Vector{T}
    ax::T
    ay::T
    exp_ax::Int
    sat_rule::Int
    density2::Union{Nothing,Array{T,4}}
end

function Base.show(io::IO, ::MIME"text/plain", fs::TJLFFluctuationSpectra{T}) where {T}
    nkx, nky, nmodes = size(fs.phi2)
    println(io, "TJLFFluctuationSpectra{$T}:")
    println(io, "  SAT_RULE=$(fs.sat_rule), ax=$(fs.ax), ay=$(fs.ay), exp_ax=$(fs.exp_ax)")
    println(io, "  kx: $nkx points in [$(round(minimum(fs.kx); sigdigits=3)), $(round(maximum(fs.kx); sigdigits=3))]")
    println(io, "  ky: $nky points in [$(round(minimum(fs.ky); sigdigits=3)), $(round(maximum(fs.ky); sigdigits=3))]")
    println(io, "  nmodes=$nmodes")
    println(io, "  phi2: Array{$T,3} of size $(size(fs.phi2))")
    if fs.density2 === nothing
        println(io, "  density2: (not computed — pass `density_weights` to include)")
    else
        println(io, "  density2: Array{$T,4} of size $(size(fs.density2)) (nkx, nky, nmodes, nspecies)")
    end
end

@inline function _sat_lorentz_shape(u, ax, ay, exp_ax)
    ay_term = (1 + ay * u^2)
    ax_term = (1 + abs(ax * u)^exp_ax)
    return 1 / (ay_term^2 * ax_term^2)
end

# Lorentzian-shaped reconstruction of `|φ|²(kx, ky, mode)` from the saturation-rule
# output. Shared by `fluctuation_spectra` (TJLF path) and `qlnn_fluctuation_spectra`
# (QLNN path) so both call sites apply the exact same `_sat_lorentz_shape` math.
#
# `phinorm`, `kx_width`, `kx0_e`, `ax`, `ay`, `exp_ax` come from
# `TJLF.intensity_sat(... ; return_phi_params=true)`.
function _build_phi2_from_sat_params(phinorm::AbstractMatrix{T},
                                     kx_width::AbstractVector{T},
                                     kx0_e::AbstractVector{T},
                                     ax::T, ay::T, exp_ax::Int,
                                     ky::AbstractVector{T};
                                     kx::Union{Nothing,AbstractVector}=nothing,
                                     n_kx::Int=128,
                                     kx_max_sigma::Real=6.0) where {T<:Real}
    nky = length(ky)
    nmodes = size(phinorm, 2)
    @assert size(phinorm, 1) == nky "_build_phi2_from_sat_params: phinorm rows ($(size(phinorm,1))) ≠ nky ($nky)"
    @assert length(kx_width) == nky "_build_phi2_from_sat_params: kx_width length ≠ nky"
    @assert length(kx0_e)    == nky "_build_phi2_from_sat_params: kx0_e length ≠ nky"

    kx_vec::Vector{T} = if kx === nothing
        lorentz_hw = ay > 0 ? kx_width ./ (ky .* sqrt(ay)) : kx_width ./ ky
        kx_half = maximum(lorentz_hw) * T(kx_max_sigma) + maximum(abs.(kx0_e))
        collect(T, range(-kx_half, kx_half; length=n_kx))
    else
        collect(T, kx)
    end
    nkx = length(kx_vec)

    phi2 = zeros(T, nkx, nky, nmodes)
    @inbounds for j in 1:nky
        inv_w = ky[j] / kx_width[j]
        u0 = -kx0_e[j] * inv_w
        L0 = _sat_lorentz_shape(u0, ax, ay, exp_ax)
        invL0 = L0 > 0 ? one(T) / L0 : one(T)
        for m in 1:nmodes
            pn = phinorm[j, m]
            pn == 0 && continue
            for i in 1:nkx
                u = (kx_vec[i] - kx0_e[j]) * inv_w
                phi2[i, j, m] = pn * _sat_lorentz_shape(u, ax, ay, exp_ax) * invL0
            end
        end
    end
    return kx_vec, phi2
end

"""
    fluctuation_spectra(tjlf_result::NamedTuple, input_tjlf::InputTJLF; kw...) -> TJLFFluctuationSpectra
    fluctuation_spectra(input_tjlf::InputTJLF; kw...)
    fluctuation_spectra(input_tglf::InputTGLF; kw...)

Reconstruct the 2D `(kx, ky, mode)` potential fluctuation spectrum `|φ|²` from a
TJLF saturation run using the Staebler SAT1/SAT2 Lorentzian kx distribution.

Only the potential spectrum is computable from standard `TJLF.run` outputs. To
additionally obtain per-species density fluctuation spectra, pass the density
QL weights via `density_weights` (shape `(nspecies, nmodes, nky)`), as produced
internally by TJLF's `get_QL_weights` (field `N_weight`).

# Arguments
- `tjlf_result`: output of `TJLF.run(input_tjlf)` (NamedTuple with `QL_weights`, `eigenvalue`, ...)
- `input_tjlf`: the `InputTJLF` used for that run (its `KY_SPECTRUM` must be populated)

# Keywords
- `kx::Union{Nothing, AbstractVector}=nothing` — custom kx grid (ρ_s units). If `nothing`, an auto-sized symmetric grid is generated covering `±kx_max_sigma` Lorentzian half-widths around the largest `|kx0_e|`.
- `n_kx::Int=128` — number of kx points when auto-sizing
- `kx_max_sigma::Real=6.0` — kx half-range in Lorentzian half-widths when auto-sizing
- `density_weights::Union{Nothing, AbstractArray}=nothing` — density QL weights `N_weight[species, mode, ky]`. If supplied, `|δnₛ|²` is stored in `density2`.

# Notes
- Requires `SAT_RULE ∈ (1, 2, 3)` and `ALPHA_QUENCH == 0.0` (spectral-shift model active). Under the quench rule (`ALPHA_QUENCH != 0`) TGLF sets `ax = ay = 0` so no Lorentzian kx shape is defined.
- `phi2_ky[j, m]` is the TJLF `phinorm` at `(ky[j], mode m)`, which equals `|φ|²` at lab-frame `kx = 0`. The spectrum peaks at `kx = kx0_e(ky[j])` with value `phi2_ky[j, m] / L(-kx0_e·ky/kx_width)`.
- Integrate `phi2` numerically in kx if a kx-integrated intensity is needed — the prefactor `L(u(0))⁻¹ · (π/(2√ay)) · (kx_width/ky)` is not absorbed into `phi2`.
"""
function fluctuation_spectra(
    tjlf_result::NamedTuple,
    input_tjlf::InputTJLF{T};
    kx::Union{Nothing,AbstractVector}=nothing,
    n_kx::Int=128,
    kx_max_sigma::Real=6.0,
    density_weights::Union{Nothing,AbstractArray}=nothing
) where {T<:Real}
    sat_rule = input_tjlf.SAT_RULE
    @assert sat_rule in (1, 2, 3) "fluctuation_spectra requires SAT_RULE ∈ {1,2,3} (got $sat_rule)"
    @assert input_tjlf.ALPHA_QUENCH == 0.0 "fluctuation_spectra requires ALPHA_QUENCH=0 (spectral-shift model)"
    @assert !isempty(input_tjlf.KY_SPECTRUM) && !any(isnan, input_tjlf.KY_SPECTRUM) "input_tjlf.KY_SPECTRUM must be populated (run TJLF first)"

    QL_weights = tjlf_result.QL_weights
    gamma_matrix = tjlf_result.eigenvalue[:, :, 1]  # (nmodes, nky)

    satParams = TJLF.get_sat_params(input_tjlf)

    # SAT2/SAT3 need zonal-mixing params from the first-pass gammas
    zonal_kwargs = if sat_rule in (2, 3)
        most_unstable_gamma_fp = gamma_matrix[1, :]
        vzf, kymax, jmax = TJLF.get_zonal_mixing(input_tjlf, satParams, most_unstable_gamma_fp)
        (; vzf_out_param=vzf, kymax_out_param=kymax, jmax_out_param=jmax)
    else
        NamedTuple()
    end

    params = TJLF.intensity_sat(
        input_tjlf, satParams, gamma_matrix, QL_weights, T(2.0), true;
        zonal_kwargs...
    )

    phinorm::Matrix{T} = params.phinorm
    kx_width::Vector{T} = params.kx_width
    kx0_e::Vector{T} = params.kx0_e
    ax::T = T(params.ax)
    ay::T = T(params.ay)
    exp_ax::Int = Int(params.exp_ax)

    ky = collect(T, input_tjlf.KY_SPECTRUM)
    nky = length(ky)
    nmodes = size(phinorm, 2)

    kx_vec, phi2 = _build_phi2_from_sat_params(phinorm, kx_width, kx0_e,
                                               ax, ay, exp_ax, ky;
                                               kx=kx, n_kx=n_kx,
                                               kx_max_sigma=kx_max_sigma)

    density2 = _density_spectrum_from_weights(density_weights, phi2, nmodes, nky, T)

    return TJLFFluctuationSpectra{T}(
        kx_vec, ky, phi2, copy(phinorm), kx_width, kx0_e, ax, ay, exp_ax, sat_rule, density2
    )
end

function fluctuation_spectra(input_tjlf::InputTJLF{T}; kw...) where {T<:Real}
    tjlf_result = TJLF.run(input_tjlf)
    return fluctuation_spectra(tjlf_result, input_tjlf; kw...)
end

function fluctuation_spectra(input_tglf::InputTGLF; kw...)
    input_tjlf = InputTJLF{Float64}(input_tglf)
    return fluctuation_spectra(input_tjlf; kw...)
end

function _density_spectrum_from_weights(::Nothing, ::Array, ::Int, ::Int, ::Type)
    return nothing
end

function _density_spectrum_from_weights(dw::AbstractArray, phi2::Array{T,3}, nmodes::Int, nky::Int, ::Type{T}) where {T<:Real}
    @assert ndims(dw) == 3 "density_weights must be a 3D array with shape (nspecies, nmodes, nky)"
    ns = size(dw, 1)
    @assert size(dw, 2) == nmodes "density_weights size(2)=$(size(dw,2)) must equal nmodes=$nmodes"
    @assert size(dw, 3) == nky    "density_weights size(3)=$(size(dw,3)) must equal nky=$nky"
    nkx = size(phi2, 1)
    d2 = Array{T,4}(undef, nkx, nky, nmodes, ns)
    @inbounds for s in 1:ns, m in 1:nmodes, j in 1:nky
        w = T(dw[s, m, j])
        for i in 1:nkx
            d2[i, j, m, s] = w * phi2[i, j, m]
        end
    end
    return d2
end

"""
    radial_correlation_length(fs; method=:hwhm, pad=4) -> (Lr, L_avg)

Radial correlation length (ρ_s units) derived from the kx power spectrum in
`fs::TJLFFluctuationSpectra`.

At each `ky[j]`, the 2-point radial autocorrelation is the inverse Fourier
transform of the (modes-summed) kx power:

    C(Δr; ky) = ∫ |φ̂|²(kx, ky) e^{i kx Δr} dkx  /  ∫ |φ̂|²(kx, ky) dkx .

The carrier from the spectral shift `kx0_e(ky)` is stripped by taking the
magnitude `|C|`; the envelope's width defines `Lr[j]`.

# Arguments
- `fs`: spectrum returned by `fluctuation_spectra`.

# Keywords
- `method = :hwhm` : half-width at half-maximum (Δr where |C| drops to 0.5).
  Alternatives: `:efold` (1/e point), `:integral` (∫|C|dΔr / |C(0)|).
- `n_dr = 512`     : number of Δr samples for the direct DFT.
- `dr_max`         : extent of the Δr grid (default 8× the longest radial
  wavelength supported by `fs.kx`).
- `modes = nothing`: optional mode slice; default sums over all modes.

# Returns
A `NamedTuple` with
- `Lr::Vector{T}` — `Lr[j]` in ρ_s units, for each `fs.ky[j]`.
- `L_avg::T`     — intensity-weighted mean, weight = `sum(phi2_ky; dims=mode)`.

Both are in ρ_s units; multiply by `ρ_s/a` (from `GACODE.rho_s` and
`eqt.boundary.minor_radius`) for units of the minor radius.
"""
function radial_correlation_length(fs::TJLFFluctuationSpectra{T};
    method::Symbol = :hwhm,
    n_dr::Int = 512,
    dr_max::Union{Nothing,Real} = nothing,
    modes::Union{Nothing,AbstractVector{<:Integer}} = nothing,
) where {T<:Real}
    kx = fs.kx
    ky = fs.ky
    length(kx) >= 4 || throw(ArgumentError("need ≥4 kx points"))
    Δkx = kx[2] - kx[1]
    @assert maximum(abs, Base.diff(kx) .- Δkx) < 1e-6 * Δkx "fs.kx must be uniformly spaced"

    # Δr grid: default span = 8× the longest supported wavelength of the kx grid.
    Δr_max = T(dr_max === nothing ? 8 * (2π / (Δkx * length(kx))) : dr_max)
    Δr_max <= 0 && throw(ArgumentError("dr_max must be > 0"))
    Δr_grid = collect(range(zero(T), Δr_max; length=n_dr))

    # modes-summed kx power at each ky
    P = if modes === nothing
        dropdims(sum(fs.phi2; dims=3); dims=3)    # (nkx, nky)
    else
        dropdims(sum(view(fs.phi2, :, :, modes); dims=3); dims=3)
    end

    Lr = similar(ky, T)
    # Precompute kx-by-Δr table of cos/sin to vectorise the direct DFT:
    #   C(Δr; ky) = Σ_kx P(kx, ky) e^{i kx Δr} Δkx
    #   |C|       = √((Σ P cos)² + (Σ P sin)²)
    cosK = cos.(kx .* Δr_grid')       # (nkx, n_dr)
    sinK = sin.(kx .* Δr_grid')

    for j in eachindex(ky)
        Pj  = @view P[:, j]
        C_re = vec(Pj' * cosK)
        C_im = vec(Pj' * sinK)
        Cmag = @. sqrt(C_re^2 + C_im^2)
        C0   = Cmag[1]
        if C0 == 0
            Lr[j] = zero(T)
            continue
        end
        Lr[j] = if method === :integral
            T(sum((Cmag[2:end] .+ Cmag[1:end-1]) .* (Δr_grid[2] / 2)) / C0)
        else
            target = method === :efold ? C0 / ℯ : T(0.5) * C0
            idx = findfirst(<(target), Cmag)
            idx === nothing ? Δr_grid[end] : Δr_grid[idx]
        end
    end

    # intensity-weighted average over ky (phi2_ky is phinorm per (ky, mode))
    w = if modes === nothing
        vec(sum(fs.phi2_ky; dims=2))
    else
        vec(sum(view(fs.phi2_ky, :, modes); dims=2))
    end
    L_avg = sum(Lr .* w) / max(sum(w), eps(T))
    return (Lr = Lr, L_avg = L_avg)
end

export TJLFFluctuationSpectra, fluctuation_spectra, radial_correlation_length

# ========== Original run_tjlf functions ==========

function run_tjlf(input_tjlf::InputTJLF{T}) where {T<:Real}
    QL_flux_out = TJLF.run_tjlf(input_tjlf)
    return GACODE.FluxSolution{T}(TJLF.Qe(QL_flux_out), TJLF.Qi(QL_flux_out), TJLF.Γe(QL_flux_out), TJLF.Γi(QL_flux_out), TJLF.Πi(QL_flux_out))
end

function run_tjlf(input_tglf::InputTGLF)
    input_tjlf = InputTJLF{Float64}(input_tglf)
    return run_tjlf(input_tjlf)
end
