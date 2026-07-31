module QLNN_CUDA_Ext

# ==========================================================================
#  GPU-accelerated QLNN inference + GPU-resident |φ|²(kx, ky) reconstruction.
#
#  This extension is loaded automatically when both `TurbulentTransport` and
#  `CUDA` are present (Julia 1.9+ Pkg-extension mechanism, registered via
#  TurbulentTransport/Project.toml `[extensions]` block).
#
#  The base package keeps its CPU-only contract (`predict`, `run_qlnn`,
#  `qlnn_fluctuation_spectra` — all pure-Julia + Flux-CPU). When a user wants
#  a GPU forward pass they call `qlnn_to_gpu(...)` first and then
#  `qlnn_fluctuation_spectra_gpu(...)`. Both names are loaded into
#  `TurbulentTransport` via stub overrides defined here.
# ==========================================================================

using TurbulentTransport
using CUDA
import Flux
import TJLF

# Re-export internal helpers from the base package.
const _build_phi2_from_sat_params = TurbulentTransport._build_phi2_from_sat_params
const log_suffix                  = TurbulentTransport.log_suffix
const _qlnn_fill_xs!              = TurbulentTransport._qlnn_fill_xs!
const _qlnn_fill_xs_with_eig!     = TurbulentTransport._qlnn_fill_xs_with_eig!
const _qlnn_pack_qlweight!        = TurbulentTransport._qlnn_pack_qlweight!
const _QLNN_TARGET_TYPE_IDX       = TurbulentTransport._QLNN_TARGET_TYPE_IDX
const _qlnn_parse_qlweight_ynames = TurbulentTransport._qlnn_parse_qlweight_ynames
const _qlnn_recover_eig           = TurbulentTransport._qlnn_recover_eig
const _qlnn_width_uses_chained_eig = TurbulentTransport._qlnn_width_uses_chained_eig
const _run_qlnn_predict            = TurbulentTransport._run_qlnn_predict

# ==========================================================================
#  GPU-resident model wrapper
# ==========================================================================

"""
    QLNNmodelGPU

GPU sibling of `QLNNmodel`: holds a `Flux.Chain` already on the GPU plus
`CuArray{Float32}` normalization buffers. Construct via `qlnn_to_gpu`.
"""
struct QLNNmodelGPU
    fluxmodel::Any                   # Flux.Chain on GPU (Flux.gpu(...))
    target::Symbol
    drop_bpar::Bool
    normalize_by_ky::Bool
    xnames::Vector{String}
    ynames::Vector{String}
    xm_d::CuVector{Float32}
    xσ_d::CuVector{Float32}
    ym_d::CuVector{Float32}
    yσ_d::CuVector{Float32}
    log_mask_d::CuVector{Float32}    # 1.0 where the feature is `_log10`-coded, else 0.0
end

"""
    QLNNensembleGPU

GPU sibling of `QLNNensemble`. `predict_gpu` averages across members.
"""
struct QLNNensembleGPU
    models::Vector{QLNNmodelGPU}
end

"""
    QLNNbundleGPU

GPU-side bundle (energy, particle, momentum, eigenvalue + optional stability
classifier + optional width regressor). Width is folded into the bundle on
the GPU side so all NN heads can be fanned out together.
"""
struct QLNNbundleGPU
    energy::Union{QLNNmodelGPU,QLNNensembleGPU}
    particle::Union{QLNNmodelGPU,QLNNensembleGPU}
    momentum::Union{QLNNmodelGPU,QLNNensembleGPU}
    eigenvalue::Union{QLNNmodelGPU,QLNNensembleGPU}
    stability::Union{Nothing,QLNNmodelGPU,QLNNensembleGPU}
    width::Union{Nothing,QLNNmodelGPU,QLNNensembleGPU}
    momentum_sign::Float64
    dir::String
end

# Forwarding accessors so callers can read xnames / target / etc. on either
# wrapper without caring whether they have a single model or an ensemble.
function Base.getproperty(ens::QLNNensembleGPU, field::Symbol)
    if field === :models
        return getfield(ens, field)
    elseif field === :fluxmodel
        error("QLNNensembleGPU: cannot access fluxmodel directly; iterate `ens.models`")
    else
        return getfield(ens.models[1], field)
    end
end

# ==========================================================================
#  Conversion: CPU model -> GPU model
# ==========================================================================

# Build the per-feature `_log10` mask once per model so the GPU normalization
# kernel can do a fused log10 + z-score pass without per-row branching.
function _make_log_mask_d(xnames::Vector{String})
    n = length(xnames)
    mask = zeros(Float32, n)
    @inbounds for i in 1:n
        mask[i] = endswith(xnames[i], log_suffix) ? 1.0f0 : 0.0f0
    end
    return CuArray(mask)
end

"""
    qlnn_to_gpu(model::QLNNmodel) -> QLNNmodelGPU
    qlnn_to_gpu(model::QLNNensemble) -> QLNNensembleGPU
    qlnn_to_gpu(bundle::QLNNbundle; width=nothing) -> QLNNbundleGPU

Move the QLNN model(s) to the active CUDA device. The `Flux.Chain` is moved
via `Flux.gpu(...)`; the per-feature normalization vectors are uploaded as
`CuArray{Float32}`. The GPU forward path stays Float32 (matching the trained
weights) and the result is cast to `Float64` on the way out, so downstream
TJLF saturation-rule code keeps its FP64 contract.

`qlnn_to_gpu(bundle; width=width_model)` lifts an optional `width_model`
into the bundle so GPU-side `qlnn_fluctuation_spectra_gpu` can fan width
inference out alongside the rest of the heads.
"""
TurbulentTransport.qlnn_to_gpu(model::TurbulentTransport.QLNNmodel) =
    _qlnn_model_to_gpu(model)

function _qlnn_model_to_gpu(model::TurbulentTransport.QLNNmodel)
    # Bypass Flux's MLDataDevices machinery (which silently falls back to CPU
    # without the cuDNN.jl trigger package and would leave the chain as
    # `Matrix{Float64}`). Functors.fmap recurses through the Chain and pushes
    # every leaf array to GPU as Float32 via `CUDA.cu`, matching the trained
    # weight precision and giving us a real `CuArray`-backed chain.
    chain_gpu = Flux.fmap(CUDA.cu, model.fluxmodel)
    return QLNNmodelGPU(
        chain_gpu,
        model.target, model.drop_bpar, model.normalize_by_ky,
        copy(model.xnames), copy(model.ynames),
        CuArray(Float32.(model.xm)),
        CuArray(Float32.(model.xσ)),
        CuArray(Float32.(model.ym)),
        CuArray(Float32.(model.yσ)),
        _make_log_mask_d(model.xnames),
    )
end

TurbulentTransport.qlnn_to_gpu(ens::TurbulentTransport.QLNNensemble) =
    QLNNensembleGPU(QLNNmodelGPU[_qlnn_model_to_gpu(m) for m in ens.models])

function TurbulentTransport.qlnn_to_gpu(bundle::TurbulentTransport.QLNNbundle;
                                        width::Union{Nothing,TurbulentTransport.AbstractQLNNmodel}=nothing)
    energy_g     = TurbulentTransport.qlnn_to_gpu(bundle.energy)
    particle_g   = TurbulentTransport.qlnn_to_gpu(bundle.particle)
    momentum_g   = TurbulentTransport.qlnn_to_gpu(bundle.momentum)
    eigenvalue_g = TurbulentTransport.qlnn_to_gpu(bundle.eigenvalue)
    stability_g  = bundle.stability === nothing ? nothing :
                   TurbulentTransport.qlnn_to_gpu(bundle.stability)
    width_g      = width === nothing ? nothing : TurbulentTransport.qlnn_to_gpu(width)
    return QLNNbundleGPU(energy_g, particle_g, momentum_g, eigenvalue_g, stability_g, width_g,
                         bundle.momentum_sign, bundle.dir)
end

# ==========================================================================
#  GPU forward pass: fused log10 + z-score + chain + denormalize
# ==========================================================================

# Build a normalized GPU input matrix in one fused broadcast.
# `x_d` (`Float32`, `(nfeat, n)`) is the raw feature matrix on device; on
# return `xx_d` is `(x_d - xm) / xσ` with optional log10 applied per-row
# according to `log_mask_d` (1.0 = log10, 0.0 = identity).
function _gpu_normalize_inputs(x_d::CuMatrix{Float32}, model::QLNNmodelGPU)
    xm = reshape(model.xm_d, :, 1)        # (nfeat, 1) for broadcast
    xσ = reshape(model.xσ_d, :, 1)
    msk = reshape(model.log_mask_d, :, 1)
    # `log10(max(v, eps)` with `eps = 1f-30` so a stray zero in the input
    # doesn't trigger -Inf and propagate NaNs through the chain.
    return @. ((1f0 - msk) * x_d + msk * log10(max(x_d, 1f-30)) - xm) / xσ
end

@inline _gpu_denormalize(yn_d::CuMatrix{Float32}, model::QLNNmodelGPU) = begin
    ym = reshape(model.ym_d, :, 1)
    yσ = reshape(model.yσ_d, :, 1)
    @. yn_d * yσ + ym
end

"""
    predict_gpu(model::QLNNmodelGPU, x_d::CuMatrix{Float32}) -> CuMatrix{Float32}
    predict_gpu(model::QLNNensembleGPU, x_d::CuMatrix{Float32}) -> CuMatrix{Float32}

GPU equivalent of `predict`. Runs log10 + z-score + chain + denormalize on
device and returns a `CuMatrix{Float32}` of shape `(length(ynames), n_samp)`.

Ensembles average member predictions; the Flux chain forward pass on GPU
already releases its scratch back to the CUDA pool.
"""
function predict_gpu(model::QLNNmodelGPU, x_d::CuMatrix{Float32})
    xx_d = _gpu_normalize_inputs(x_d, model)
    yn_d = model.fluxmodel(xx_d)
    return _gpu_denormalize(yn_d, model)
end

function predict_gpu(ens::QLNNensembleGPU, x_d::CuMatrix{Float32})
    M = length(ens.models)
    @assert M > 0 "predict_gpu: empty ensemble"
    acc = predict_gpu(ens.models[1], x_d)
    if M == 1
        return acc
    end
    for k in 2:M
        acc = acc .+ predict_gpu(ens.models[k], x_d)
    end
    return acc ./ Float32(M)
end

# ==========================================================================
#  GPU helper: build phi2 from intensity_sat parameters
# ==========================================================================

# Inline GPU shape function (mirrors `_sat_lorentz_shape` from
# TurbulentTransport/src/tjlf.jl). The CUDA broadcast / kernel below uses
# this same expression element-wise.
@inline _gpu_lorentz(u::Float64, ax::Float64, ay::Float64, exp_ax::Int) =
    1.0 / ((1.0 + ay * u^2)^2 * (1.0 + abs(ax * u)^exp_ax)^2)

# Build `phi2(kx, ky, mode)` on the GPU. `phinorm`, `kx_width`, `kx0_e`,
# `ky` all live on host (saturation-rule outputs are tiny; the heavy data
# is the Lorentzian-shaped tensor itself, which we keep on device).
function _build_phi2_gpu(phinorm::AbstractMatrix{Float64},
                          kx_width::AbstractVector{Float64},
                          kx0_e::AbstractVector{Float64},
                          ax::Float64, ay::Float64, exp_ax::Int,
                          ky::AbstractVector{Float64};
                          n_kx::Int = 128, kx_max_sigma::Real = 6.0,
                          kx::Union{Nothing,AbstractVector} = nothing)
    nky = length(ky)
    nmodes = size(phinorm, 2)
    @assert size(phinorm, 1) == nky "_build_phi2_gpu: phinorm rows ≠ nky"

    kx_vec = if kx === nothing
        lorentz_hw = ay > 0 ? kx_width ./ (ky .* sqrt(ay)) : kx_width ./ ky
        kx_half = maximum(lorentz_hw) * kx_max_sigma + maximum(abs.(kx0_e))
        collect(Float64, range(-kx_half, kx_half; length=n_kx))
    else
        collect(Float64, kx)
    end
    nkx = length(kx_vec)

    # Precompute per-ky scalars on host: inv_w[j], invL0[j], pn[j, m].
    inv_w  = ky ./ kx_width                   # (nky,)
    u0     = -kx0_e .* inv_w                  # (nky,)
    L0     = [ _gpu_lorentz(u0[j], ax, ay, exp_ax) for j in 1:nky ]
    invL0  = [ L > 0 ? 1.0 / L : 1.0 for L in L0 ]

    # Upload to GPU and broadcast the Lorentzian shape over (nkx, nky, nmodes).
    kx_d     = CuArray(kx_vec)                # (nkx,)
    inv_w_d  = CuArray(inv_w)                 # (nky,)
    kx0_e_d  = CuArray(kx0_e)                 # (nky,)
    invL0_d  = CuArray(invL0)                 # (nky,)
    phinorm_d = CuArray(phinorm)              # (nky, nmodes)
    kxr  = reshape(kx_d, nkx, 1, 1)
    kw   = reshape(inv_w_d, 1, nky, 1)
    kx0r = reshape(kx0_e_d, 1, nky, 1)
    iL0  = reshape(invL0_d, 1, nky, 1)
    pnr  = reshape(phinorm_d, 1, nky, nmodes)

    # The broadcast below keeps everything on device (no host scratch).
    # `exp_ax` is hoisted out of the broadcast as a `Float64` literal so
    # CUDA.jl's integer-exponent fastpath fuses cleanly into the kernel.
    exp_ax_f = Float64(exp_ax)
    phi2_d = @. pnr * iL0 / ((1.0 + ay * ((kxr - kx0r) * kw)^2)^2 *
                              (1.0 + abs(ax * ((kxr - kx0r) * kw))^exp_ax_f)^2)

    return kx_d, phi2_d
end

# ==========================================================================
#  qlnn_fluctuation_spectra_gpu: full Phase A on GPU + GPU phi_amp injection
# ==========================================================================

"""
    qlnn_fluctuation_spectra_gpu(input_tjlfs::Vector{InputTJLF},
                                 bundle_gpu::QLNNbundleGPU;
                                 stability_threshold=0.5,
                                 n_kx=128, kx_max_sigma=6.0,
                                 kx=nothing) -> Vector{NamedTuple}

GPU sibling of `qlnn_fluctuation_spectra`. NN forward passes run on the
active CUDA device; the saturation-rule call (`TJLF.intensity_sat`) is on
the CPU (small + branchy); the Lorentzian `phi2` and `phi_amp = sqrt(sum
phi2; dims=3)` build is back on GPU so downstream
`build_frame_rasters_gpu` can consume `phi_amp_d` without an H2D copy.

Per-radial NamedTuple fields:
  - `fs::TJLFFluctuationSpectra{Float64}`     (CPU mirror; same shape as the CPU path)
  - `eigenvalue::Array{Float64,3}`            `(1, nky, 2)` (γ, ω) on CPU
  - `phi_amp_d::CuMatrix{Float64}`            `(nkx, nky)` on device
  - `kx_d::CuVector{Float64}`                 on device
  - `ky_d::CuVector{Float64}`                 on device
  - `ω_ky_d::CuVector{Float64}`               physical ω per ky, on device
"""
function TurbulentTransport.qlnn_fluctuation_spectra_gpu(
        input_tjlfs::Vector{TJLF.InputTJLF{T}},
        bundle_gpu::QLNNbundleGPU;
        stability_threshold::Real = 0.5,
        n_kx::Int = 128,
        kx_max_sigma::Real = 6.0,
        kx::Union{Nothing,AbstractVector} = nothing) where {T<:Real}

    nr = length(input_tjlfs)
    if nr == 0
        return NamedTuple[]
    end

    # SAT_RULE / ALPHA_QUENCH precondition (mirrors the CPU path).
    for (r, it) in enumerate(input_tjlfs)
        @assert it.SAT_RULE in (1, 2, 3) "qlnn_fluctuation_spectra_gpu requires SAT_RULE ∈ {1,2,3} (got $(it.SAT_RULE) at r=$r)"
        @assert it.ALPHA_QUENCH == 0.0 "qlnn_fluctuation_spectra_gpu requires ALPHA_QUENCH=0 (r=$r)"
    end

    # ----- Phase 1: CPU per-radial setup (sat_params + KY_SPECTRUM + xs_all) -----
    # We don't have a GPU-resident equivalent, but this is tiny relative to
    # the NN forward + Lorentzian build.
    info_e = _qlnn_parse_qlweight_ynames(bundle_gpu.energy.ynames)
    info_p = _qlnn_parse_qlweight_ynames(bundle_gpu.particle.ynames)
    info_m = _qlnn_parse_qlweight_ynames(bundle_gpu.momentum.ynames)
    nf = info_e.nf
    ns = info_e.ns

    energy_xnames = bundle_gpu.energy.xnames

    nky_r = Vector{Int}(undef, nr)
    sat_params_v = Vector{Any}(undef, nr)
    ky_spectrums = Vector{Any}(undef, nr)
    for r in 1:nr
        it = input_tjlfs[r]
        it.NMODES = 1
        sp = TJLF.get_sat_params(it)
        ks = TJLF.get_ky_spectrum(it, sp.grad_r0)
        sat_params_v[r] = sp
        ky_spectrums[r] = ks
        nky_r[r] = length(ks)
        if it.KY_SPECTRUM === missing || length(it.KY_SPECTRUM) != nky_r[r]
            it.KY_SPECTRUM = collect(ks)
        else
            it.KY_SPECTRUM .= ks
        end
    end

    total_nky = sum(nky_r)
    chunk_starts = Vector{Int}(undef, nr)
    let off = 0
        for r in 1:nr
            chunk_starts[r] = off + 1
            off += nky_r[r]
        end
    end

    # ----- Phase 2: build feature matrix on CPU at the InputTJLF eltype, then -----
    # cast to Float32 for the upload. `_qlnn_fill_xs!` has the parametric
    # constraint `xs::AbstractMatrix{T}` + `input_tjlf::InputTJLF{T}` so we
    # build at T (typically Float64) and convert in the H2D step.
    nfeat = length(energy_xnames)
    xs_all = Matrix{T}(undef, nfeat, total_nky)
    for r in 1:nr
        c0 = chunk_starts[r]
        nk = nky_r[r]
        view_block = view(xs_all, :, c0:c0+nk-1)
        _qlnn_fill_xs!(view_block, input_tjlfs[r],
                       ky_spectrums[r], energy_xnames)
    end
    xs_d = CuArray(Float32.(xs_all))

    # ----- Phase 3: NN forward on GPU. -----
    Y_energy_d   = predict_gpu(bundle_gpu.energy,     xs_d)
    Y_particle_d = predict_gpu(bundle_gpu.particle,   xs_d)
    Y_momentum_d = predict_gpu(bundle_gpu.momentum,   xs_d)
    Y_eig_d      = predict_gpu(bundle_gpu.eigenvalue, xs_d)
    P_unstable_d = bundle_gpu.stability === nothing ? nothing :
                   begin
                       # Stability head outputs `(1, n_samp)` raw activations;
                       # apply σ (logistic) to get unstable probabilities.
                       raw = predict_gpu(bundle_gpu.stability, xs_d)
                       1.0f0 ./ (1.0f0 .+ exp.(-raw))
                   end

    # Pull predictions back to host for QL packing + intensity_sat (CPU code).
    # These tensors are tiny (≤ MB).
    Y_energy   = Float64.(Array(Y_energy_d))
    Y_particle = Float64.(Array(Y_particle_d))
    # Apply the per-bundle toroidal-stress sign convention (see QLNNbundle).
    # The CPU path applies the same factor in `_run_qlnn_predict`.
    Y_momentum = bundle_gpu.momentum_sign .* Float64.(Array(Y_momentum_d))
    Y_eig      = Float64.(Array(Y_eig_d))
    P_unstable = P_unstable_d === nothing ? nothing : vec(Float64.(Array(P_unstable_d)))

    # ----- Optional width head (chained EV inputs). ---------------------
    width_predictions = Vector{Union{Nothing,Vector{Float64}}}(undef, nr)
    width_uses_chained = bundle_gpu.width !== nothing &&
        _qlnn_width_uses_chained_eig(bundle_gpu.width)
    if width_uses_chained
        width_xnames = bundle_gpu.width.xnames
        for r in 1:nr
            c0 = chunk_starts[r]
            nk = nky_r[r]
            cols = c0:c0+nk-1
            ks = ky_spectrums[r]
            γ_phys, ω_phys = _qlnn_recover_eig(view(Y_eig, :, cols), ks,
                                               bundle_gpu.eigenvalue.normalize_by_ky)
            xs_w = Matrix{T}(undef, length(width_xnames), nk)
            _qlnn_fill_xs_with_eig!(xs_w, input_tjlfs[r], ks,
                                    width_xnames, γ_phys, ω_phys)
            xs_w_d = CuArray(Float32.(xs_w))
            y_w_d  = predict_gpu(bundle_gpu.width, xs_w_d)
            y_w    = Float64.(Array(y_w_d))
            width_predictions[r] = Float64[max(y_w[1, j], 1e-3) for j in 1:nk]
        end
    else
        for r in 1:nr
            width_predictions[r] = nothing
        end
    end

    # ----- Phase 4: per-radial CPU pack -> intensity_sat -> GPU phi2 -----
    out = Vector{NamedTuple}(undef, nr)
    for r in 1:nr
        c0 = chunk_starts[r]
        nk = nky_r[r]
        cols = c0:c0+nk-1
        ks = ky_spectrums[r]
        it = input_tjlfs[r]
        sat_params = sat_params_v[r]

        mask = nothing
        if P_unstable !== nothing
            mask = Bool[P_unstable[c0+j-1] >= stability_threshold for j in 1:nk]
        end

        QL = zeros(Float64, nf, ns, 1, nk, 5)
        _qlnn_pack_qlweight!(QL, view(Y_energy,   :, cols),
                             info_e, _QLNN_TARGET_TYPE_IDX[:energy],   ks, mask)
        _qlnn_pack_qlweight!(QL, view(Y_particle, :, cols),
                             info_p, _QLNN_TARGET_TYPE_IDX[:particle], ks, mask)
        _qlnn_pack_qlweight!(QL, view(Y_momentum, :, cols),
                             info_m, _QLNN_TARGET_TYPE_IDX[:momentum], ks, mask)

        γ_phys, ω_phys = _qlnn_recover_eig(view(Y_eig, :, cols), ks,
                                           bundle_gpu.eigenvalue.normalize_by_ky)

        Γ = zeros(Float64, 1, nk)
        for k in 1:nk
            γ = γ_phys[k]
            if mask !== nothing && !mask[k]
                γ = 0.0
            end
            Γ[1, k] = γ
        end

        if width_predictions[r] !== nothing
            wpred = width_predictions[r]
            if it.WIDTH_SPECTRUM === missing || length(it.WIDTH_SPECTRUM) != nk
                it.WIDTH_SPECTRUM = collect(Float64, wpred)
            else
                it.WIDTH_SPECTRUM .= wpred
            end
        else
            if it.WIDTH_SPECTRUM === missing || length(it.WIDTH_SPECTRUM) != nk
                it.WIDTH_SPECTRUM = fill(Float64(it.WIDTH), nk)
            else
                fill!(it.WIDTH_SPECTRUM, Float64(it.WIDTH))
            end
        end

        sat_rule = it.SAT_RULE
        zonal_kwargs = if sat_rule in (2, 3)
            most_unstable = Γ[1, :]
            vzf, kymax, jmax = TJLF.get_zonal_mixing(it, sat_params, most_unstable)
            (; vzf_out_param=vzf, kymax_out_param=kymax, jmax_out_param=jmax)
        else
            NamedTuple()
        end

        params = TJLF.intensity_sat(
            it, sat_params, Γ, QL, 2.0, true;
            zonal_kwargs...
        )

        phinorm  = Matrix{Float64}(params.phinorm)
        kx_width = Vector{Float64}(params.kx_width)
        kx0_e    = Vector{Float64}(params.kx0_e)
        ax       = Float64(params.ax)
        ay       = Float64(params.ay)
        exp_ax   = Int(params.exp_ax)

        ky_vec = collect(Float64, ks)

        # Build phi2 on the GPU.
        kx_d, phi2_d = _build_phi2_gpu(phinorm, kx_width, kx0_e, ax, ay, exp_ax, ky_vec;
                                       n_kx=n_kx, kx_max_sigma=kx_max_sigma, kx=kx)

        # phi_amp = sqrt(sum(phi2; dims=3)) — on device, with NMODES=1.
        nkx_b = size(phi2_d, 1)
        phi_amp_d = sqrt.(max.(reshape(sum(phi2_d; dims=3), nkx_b, length(ky_vec)), 0.0))

        # CPU mirror of fs (so notebook code that reads fs.kx, fs.ky, fs.phi2
        # keeps working — used by `assemble_band_inputs` to fill bands_meta).
        kx_h    = Array(kx_d)
        phi2_h  = Array(phi2_d)
        fs = TurbulentTransport.TJLFFluctuationSpectra{Float64}(
            kx_h, ky_vec, phi2_h, copy(phinorm), kx_width, kx0_e,
            ax, ay, exp_ax, sat_rule, nothing,
        )

        eig_arr = zeros(Float64, 1, nk, 2)
        for k in 1:nk
            eig_arr[1, k, 1] = γ_phys[k]
            eig_arr[1, k, 2] = ω_phys[k]
        end

        # GPU-resident ω, ky for `build_frame_rasters_gpu` consumption.
        ω_ky_d = CuArray(Vector{Float64}(ω_phys))
        ky_d   = CuArray(ky_vec)

        out[r] = (; fs, eigenvalue=eig_arr,
                  phi_amp_d, ω_ky_d, kx_d, ky_d)
    end

    return out
end

# ==========================================================================
#  __init__
# ==========================================================================

function __init__()
    @debug "QLNN_CUDA_Ext loaded — qlnn_to_gpu / qlnn_fluctuation_spectra_gpu live"
    return nothing
end

end # module QLNN_CUDA_Ext
