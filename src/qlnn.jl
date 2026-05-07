# ============================================================================
#  QLNN inference + TJLF saturation-rule integration
#
#  Lightweight inference path for QLNN (Quasi-Linear NN) bundles. Loads BSON
#  files saved by the (external) training pipeline and runs them through the
#  TJLF saturation rule (`intensity_sat` + `sum_ky_spectrum`) so that fluxes
#  still respect `SAT_RULE` / `ALPHA_ZF` / `KYGRID_MODEL` from the InputTJLF.
#
#  The training pipeline lives in a separate package and pulls in heavy
#  training-time deps (Zygote, Plots, ChainRulesCore, ...). We don't want
#  those leaking into FUSE's transport hot-path, so we keep only the minimal
#  struct surface and BSON load logic here. The on-disk dict layout is the
#  canonical one written by the trainer, so the same .bson files load both
#  here and there.
# ============================================================================

import LinearAlgebra: BLAS

#= ============================================================ =#
#  Bundle layout / file naming on disk
#= ============================================================ =#

# Canonical species/field ordering used by every QLNN regressor (frozen at
# training time). These are also the ordering TJLF expects when packing
# QL_weights for `sum_ky_spectrum` (electrons first; matches `get_sat_params`).
const _QLNN_SPECIES = ("e", "D", "C")
const _QLNN_FIELDS = ("phi", "apar", "bpar")  # bpar is optional (drop_bpar=true)

# regressor file names inside a bundle dir (e.g. models/QLNN/<name>.bson)
const _QLNN_REGRESSOR_FILES = (
    energy = "energy_regressor.bson",
    particle = "particle_regressor.bson",
    momentum = "momentum_regressor.bson",
    eigenvalue = "eigenvalue_regressor.bson",
)
const _QLNN_STABILITY_FILE = "stability_classifier.bson"

# Map regressor target Symbol -> TJLF type-axis index.
# (TJLF type axis is: 1=particle, 2=energy, 3=tor stress, 4=par stress, 5=exchange.)
const _QLNN_TARGET_TYPE_IDX = Dict(:particle => 1, :energy => 2, :momentum => 3)

#= ============================================================ =#
#  Structs (minimal inference-side mirror of the trained model)
#= ============================================================ =#

abstract type AbstractQLNNmodel end

"""
    QLNNmodel

Minimal inference-side mirror of the trained QL-weight NN model. Holds only
the fields required to run the trained NN forward and unwind the
training-time normalization / `_over_ky` / `normalize_by_ky` transforms.

The `_pooled_chain` field caches a `PooledChain` (built at load time via
`poolify(fluxmodel)`) so that `predict!` can run the Flux chain with zero
heap allocations after warmup — same trick TGLFNN uses in its hot path.
The pool is per-thread (task-local) so threaded ensemble inference is safe.
"""
struct QLNNmodel <: AbstractQLNNmodel
    fluxmodel::Flux.Chain
    target::Symbol            # :energy, :particle, :momentum, :eigenvalue, :stability
    drop_bpar::Bool
    normalize_by_ky::Bool
    xnames::Vector{String}
    ynames::Vector{String}
    xm::Vector{Float32}
    xσ::Vector{Float32}
    ym::Vector{Float32}
    yσ::Vector{Float32}
    # Per-feature [min, max] of the training set in pre-normalized, post-`_log10`
    # space (i.e. directly comparable to `_qlnn_apply_log10!(value)`). Used by
    # `warn_nn_train_bounds`. Filled with ±Inf when the BSON didn't carry an
    # `:xbounds` field, in which case the bounds check is a silent no-op.
    xbounds::Matrix{Float32}
    _pooled_chain::PooledChain
end

function Base.show(io::IO, ::MIME"text/plain", m::QLNNmodel)
    println(io, "QLNNmodel(target=$(m.target), nx=$(length(m.xnames)), ny=$(length(m.ynames)), normalize_by_ky=$(m.normalize_by_ky))")
end

"""
    QLNNensemble

Wraps a vector of `QLNNmodel`s; `predict` averages the per-member outputs
(mean of physical predictions for regressors; for classifiers we go through
[`predict_unstable_prob`](@ref) which averages probabilities post-σ).
"""
struct QLNNensemble <: AbstractQLNNmodel
    models::Vector{QLNNmodel}
end

function Base.show(io::IO, mime::MIME"text/plain", ens::QLNNensemble)
    println(io, "QLNNensemble (n=$(length(ens.models)))")
    return show(io, mime, ens.models[1])
end

# Forward most metadata accessors to the first member (xnames/target are the
# same across an ensemble by construction).
function Base.getproperty(ens::QLNNensemble, field::Symbol)
    if field === :models
        return getfield(ens, field)
    elseif field === :fluxmodel
        error("QLNNensemble: cannot access fluxmodel directly; iterate `ens.models`")
    else
        return getfield(ens.models[1], field)
    end
end

"""
    QLNNbundle

Bundle of regressors + (optional) stability classifier loaded from a
single directory. `stability === nothing` means the directory had no
`stability_classifier.bson`, in which case `run_qlnn` skips the hard gate.
"""
struct QLNNbundle
    energy::AbstractQLNNmodel
    particle::AbstractQLNNmodel
    momentum::AbstractQLNNmodel
    eigenvalue::AbstractQLNNmodel
    stability::Union{Nothing,AbstractQLNNmodel}
    dir::String
end

function Base.show(io::IO, ::MIME"text/plain", b::QLNNbundle)
    println(io, "QLNNbundle (dir=$(b.dir))")
    println(io, "  energy:     ", _summary_line(b.energy))
    println(io, "  particle:   ", _summary_line(b.particle))
    println(io, "  momentum:   ", _summary_line(b.momentum))
    println(io, "  eigenvalue: ", _summary_line(b.eigenvalue))
    println(io, "  stability:  ", b.stability === nothing ? "<none>" : _summary_line(b.stability))
end

_summary_line(m::QLNNmodel) = "QLNNmodel(target=$(m.target))"
_summary_line(m::QLNNensemble) = "QLNNensemble(n=$(length(m.models)), target=$(m.models[1].target))"
_summary_line(::Nothing) = "<none>"

#= ============================================================ =#
#  BSON loading
#= ============================================================ =#

# Build a QLNNmodel from the on-disk dict layout written by the trainer.
# We pull the inference-relevant fields plus `xbounds` (used by
# `warn_nn_train_bounds`), and ignore `name`, `date`, `source`, `ybounds` so
# the inference-side struct stays small.
function _qlnn_dict2mod(dict::AbstractDict)
    fluxmodel = Flux.fmap(Flux.f64, dict[:fluxmodel])
    target = Symbol(dict[:target])
    drop_bpar = haskey(dict, :drop_bpar) ? Bool(dict[:drop_bpar]) : true
    normalize_by_ky = haskey(dict, :normalize_by_ky) ? Bool(dict[:normalize_by_ky]) : false
    xnames = String[String(x) for x in dict[:xnames]]
    ynames = String[String(y) for y in dict[:ynames]]
    xm = Float32.(vec(dict[:xm]))
    xσ = Float32.(vec(dict[:xσ]))
    ym = Float32.(vec(dict[:ym]))
    yσ = Float32.(vec(dict[:yσ]))
    nfeat = length(xnames)
    xbounds = _qlnn_load_xbounds(dict, nfeat)
    pooled = PooledChain(poolify(fluxmodel))
    return QLNNmodel(fluxmodel, target, drop_bpar, normalize_by_ky, xnames, ynames, xm, xσ, ym, yσ, xbounds, pooled)
end

# Load training-time per-feature [min, max] bounds from the trainer dict.
# Accepts either a `(nfeat, 2)` matrix or a `2×nfeat` matrix and normalizes
# to `(nfeat, 2)`. Returns ±Inf bounds if the dict has no `:xbounds` key, so
# the bounds check downstream becomes a silent no-op for old bundles.
function _qlnn_load_xbounds(dict::AbstractDict, nfeat::Int)
    if !haskey(dict, :xbounds)
        return hcat(fill(Float32(-Inf), nfeat), fill(Float32(Inf), nfeat))
    end
    xb_raw = dict[:xbounds]
    xb = Float32.(xb_raw)
    if size(xb) == (nfeat, 2)
        return xb
    elseif size(xb) == (2, nfeat)
        return permutedims(xb)
    else
        error("QLNN: unexpected `:xbounds` shape $(size(xb_raw)); expected ($nfeat, 2) or (2, $nfeat).")
    end
end

function _qlnn_dict2ens(dict::AbstractDict)
    # Ensemble dicts use Integer keys. Sort so the iteration order is deterministic.
    keys_sorted = sort(collect(keys(dict)))
    members = QLNNmodel[_qlnn_dict2mod(dict[k]) for k in keys_sorted]
    return QLNNensemble(members)
end

"""
    loadqlnnmodel(filename) -> QLNNmodel | QLNNensemble

Resolve `filename` (relative to a registered model search path or built-in
`models/` dir) and load a single regressor/classifier `.bson`. Dispatches
to ensemble vs single-model layout based on the dict key type, matching
[`loadmodel`](@ref) for TGLFNN bundles.
"""
function loadqlnnmodel(filename::AbstractString)
    fullpath = resolve_model_path(filename; extensions=[".bson"])
    savedict = BSON.load(fullpath, @__MODULE__)
    if typeof(first(keys(savedict))) <: Integer
        return _qlnn_dict2ens(savedict)
    else
        return _qlnn_dict2mod(savedict)
    end
end

"Memoize-wrapped variant for repeated lookups in the transport hot path."
Memoize.@memoize function loadqlnnmodelonce(filename::String)
    return loadqlnnmodel(filename)
end

"""
    _resolve_qlnn_dir(name) -> String

Resolve a QLNN bundle directory (`name`) by:
1. Treating it as an absolute or already-existing relative path,
2. Searching the registered model search paths in order,
3. Falling back to `<TurbulentTransport>/models/<name>`.
"""
function _resolve_qlnn_dir(name::AbstractString)
    isdir(name) && return String(name)
    for search_dir in _MODEL_SEARCH_PATHS
        candidate = joinpath(search_dir, name)
        isdir(candidate) && return candidate
    end
    builtin = joinpath(dirname(@__DIR__), "models", name)
    isdir(builtin) && return builtin
    error("QLNN bundle directory '$name' not found. Searched provider paths " *
          "and built-in `models/` directory.")
end

"""
    loadqlnnbundle(name::AbstractString) -> QLNNbundle

Load the four regressor BSONs (energy/particle/momentum/eigenvalue) and the
optional stability classifier from a directory named `name` (resolved via
`_resolve_qlnn_dir`). Each individual file is memoized via
`loadqlnnmodelonce`, so subsequent reloads are cheap.
"""
function loadqlnnbundle(name::AbstractString)
    base = _resolve_qlnn_dir(name)
    energy = loadqlnnmodelonce(joinpath(base, _QLNN_REGRESSOR_FILES.energy))
    particle = loadqlnnmodelonce(joinpath(base, _QLNN_REGRESSOR_FILES.particle))
    momentum = loadqlnnmodelonce(joinpath(base, _QLNN_REGRESSOR_FILES.momentum))
    eigval = loadqlnnmodelonce(joinpath(base, _QLNN_REGRESSOR_FILES.eigenvalue))
    stab_path = joinpath(base, _QLNN_STABILITY_FILE)
    stability = isfile(stab_path) ? loadqlnnmodelonce(stab_path) : nothing
    return QLNNbundle(energy, particle, momentum, eigval, stability, base)
end

"Memoize-wrapped bundle loader to avoid repeated path resolution."
Memoize.@memoize function loadqlnnbundleonce(name::String)
    return loadqlnnbundle(name)
end

#= ============================================================ =#
#  Forward inference (predict + classifier helper)
#= ============================================================ =#

# Apply training-time `_log10` transform to rows whose xname ends in `_log10`.
# The training data was pre-transformed with `log10(...)` for those features
# (they span many orders of magnitude — BETAE, DEBYE, XNUE, ...); at inference
# we get InputTJLF values in linear units so we have to apply the same
# transform here (matches `flux_array!` in tglf_nn.jl).
#
# Operates in place on `xx` and is generic in eltype so that ForwardDiff `Dual`s
# flow through unchanged (`log10(::Dual)` is defined).
function _qlnn_apply_log10!(xx::AbstractMatrix, xnames::Vector{String})
    @inbounds for i in eachindex(xnames)
        if endswith(xnames[i], log_suffix)
            for j in axes(xx, 2)
                xx[i, j] = log10(xx[i, j])
            end
        end
    end
    return xx
end

function _qlnn_apply_log10!(xx::AbstractVector, xnames::Vector{String})
    @inbounds for i in eachindex(xnames)
        if endswith(xnames[i], log_suffix)
            xx[i] = log10(xx[i])
        end
    end
    return xx
end

# Mirror of TGLFNN's `flux_array!` train-bounds check, but scoped to the QLNN
# pipeline and run once per `run_qlnn` call (all bundle members share xnames,
# so scanning each head would produce 4–5× duplicate warnings).
#
# `xs_all` is the linear-space, batched feature matrix for one whole call.
# Per-feature we apply the same `_log10` transform the trainer applied (rows
# whose name ends in `_log10`) and report the most extreme min/max offender
# against `model.xbounds[i, 1:2]`. NaN/Inf are hard errors — they would also
# blow up the forward pass, so we'd rather fail with a precise message here.
#
# Bounds are stored as Float32. Comparisons promote, so `T = ForwardDiff.Dual`
# inputs work too: `Dual <: Real` defines `<` against Float32 by reading the
# value field. The warning's printed value is whatever `xs_all[i, jext]` is —
# `Float64` or a `Dual` — same convention as `flux_array!` in `tglf_nn.jl`.
_qlnn_check_xs_bounds(xs_all::AbstractMatrix, ens::QLNNensemble) =
    _qlnn_check_xs_bounds(xs_all, ens.models[1])

function _qlnn_check_xs_bounds(xs_all::AbstractMatrix, model::QLNNmodel)
    xbounds = model.xbounds
    xnames = model.xnames
    @assert size(xs_all, 1) == length(xnames) "_qlnn_check_xs_bounds: xs_all rows ($(size(xs_all,1))) ≠ length(xnames)=$(length(xnames))"
    nfeat, nsamp = size(xs_all)
    nsamp == 0 && return nothing
    @inbounds for i in 1:nfeat
        lo = xbounds[i, 1]
        hi = xbounds[i, 2]
        # ±Inf bounds (legacy bundles without :xbounds) → skip silently.
        if !isfinite(lo) && !isfinite(hi)
            continue
        end
        log_row = endswith(xnames[i], log_suffix)
        # Track value-space min/max so a single warning per feature reports
        # the worst offender across the whole batch instead of just the first.
        v0 = xs_all[i, 1]
        vmin_xform = log_row ? log10(v0) : v0
        vmax_xform = vmin_xform
        jmin = 1
        jmax = 1
        for j in 1:nsamp
            v = xs_all[i, j]
            if isnan(v) || isinf(v)
                error("QLNN input $(xnames[i]) = $v at sample $j is not allowed (NaN/Inf).")
            end
            v_xform = log_row ? log10(v) : v
            if v_xform < vmin_xform
                vmin_xform = v_xform
                jmin = j
            end
            if v_xform > vmax_xform
                vmax_xform = v_xform
                jmax = j
            end
        end
        if vmin_xform < lo
            @warn "Extrapolation $(xnames[i])=$(xs_all[i, jmin]) is below training bound of $lo (transformed value $vmin_xform)"
        elseif vmax_xform > hi
            @warn "Extrapolation $(xnames[i])=$(xs_all[i, jmax]) is above training bound of $hi (transformed value $vmax_xform)"
        end
    end
    return nothing
end

"""
    predict!(out, model::QLNNmodel, x) -> out

In-place forward pass of one regressor / classifier on `x`. Mirrors TGLFNN's
`flux_array!`:

1. Acquire a thread-local pool buffer for the normalized input.
2. Apply `_log10` per-row + `(x - xm) / xσ` in a single fused loop.
3. Run the pre-built `PooledChain` directly into `out` (no GEMM allocs).
4. Denormalize `out .* yσ .+ ym` in place.

Inputs are row-major `(n_features, n_samples)` matrices (or a length-`nf`
vector for the single-sample variant). Eltype `T = eltype(x)` is preserved
end-to-end so ForwardDiff `Dual` partials flow through (the pool acquires
`T`-typed buffers via `acquire!(pool, T, ...)`).
"""
@with_pool pool function predict!(out::AbstractMatrix{T}, model::QLNNmodel,
                                  x::AbstractMatrix{T}) where {T<:Real}
    @assert size(out, 1) == length(model.ynames) "predict!: output rows ($(size(out,1))) ≠ length(ynames)=$(length(model.ynames))"
    @assert size(out, 2) == size(x, 2) "predict!: output cols ($(size(out,2))) ≠ size(x,2)=$(size(x,2))"
    nfeat, nsamp = size(x)
    xx = acquire!(pool, T, nfeat, nsamp)
    xnames = model.xnames
    xm = model.xm
    xσ = model.xσ
    @inbounds for i in 1:nfeat
        log_row = endswith(xnames[i], log_suffix)
        xm_i = xm[i]
        xσ_i = xσ[i]
        for j in 1:nsamp
            v = x[i, j]
            if log_row
                v = log10(v)
            end
            xx[i, j] = (v - xm_i) / xσ_i
        end
    end
    model._pooled_chain(out, xx)
    ym = model.ym
    yσ = model.yσ
    @inbounds for i in axes(out, 1)
        ym_i = ym[i]
        yσ_i = yσ[i]
        for j in axes(out, 2)
            out[i, j] = out[i, j] * yσ_i + ym_i
        end
    end
    return out
end

@with_pool pool function predict!(out::AbstractVector{T}, model::QLNNmodel,
                                  x::AbstractVector{T}) where {T<:Real}
    @assert length(out) == length(model.ynames) "predict!: output length ($(length(out))) ≠ length(ynames)=$(length(model.ynames))"
    nfeat = length(x)
    xx = acquire!(pool, T, nfeat)
    xnames = model.xnames
    xm = model.xm
    xσ = model.xσ
    @inbounds for i in 1:nfeat
        v = x[i]
        if endswith(xnames[i], log_suffix)
            v = log10(v)
        end
        xx[i] = (v - xm[i]) / xσ[i]
    end
    model._pooled_chain(out, xx)
    @inbounds for i in eachindex(out)
        out[i] = out[i] * model.yσ[i] + model.ym[i]
    end
    return out
end

"""
    predict(model::QLNNmodel, x) -> AbstractArray

Allocating wrapper around `predict!`. Output shape is `length(ynames)` for
a vector input, `(length(ynames), n_samples)` for a matrix input. Eltype
follows `eltype(x)` so ForwardDiff `Dual` inputs propagate through.
"""
function predict(model::QLNNmodel, x::AbstractMatrix{T}) where {T<:Real}
    out = Matrix{T}(undef, length(model.ynames), size(x, 2))
    return predict!(out, model, x)
end

function predict(model::QLNNmodel, x::AbstractVector{T}) where {T<:Real}
    out = Vector{T}(undef, length(model.ynames))
    return predict!(out, model, x)
end

# Threading threshold: only fan out across threads when there's enough work
# to amortize task-spawn overhead. Single-threaded Julia (or single-member
# "ensembles") falls through to the sequential path.
_qlnn_should_thread(nmodels::Int) = nmodels > 1 && Threads.nthreads() > 1

# ---------------------------------------------------------------------------
#  Thread-indexed pool array for zero-alloc threaded ensemble inference
# ---------------------------------------------------------------------------
# One AdaptiveArrayPool per Julia thread, populated at runtime via
# `_qlnn_ensure_thread_pools_capacity!`. The threaded ensemble path uses
# `_qlnn_get_thread_pool()` (which routes via `Threads.threadid()`) instead
# of `get_task_local_pool()`, so pool buffers survive across
# `Threads.@threads` calls (task-local pools are lost when tasks end, causing
# cold-start allocations on every `run_qlnn` call).
#
# Why runtime, not precompile-time: a module-load `[AdaptiveArrayPool() for
# _ in 1:Threads.nthreads()]` evaluates `Threads.nthreads()` *during
# precompilation*, which is single-threaded — so the cached `.ji` ships a
# 1-element vector and indexing it from any other thread raises
# `BoundsError`. We initialize lazily instead, with a one-time resize during
# `TurbulentTransport.__init__()` (see `TurbulentTransport.jl`) plus a
# locked grow-on-demand fallback in case `Threads.threadid()` exceeds the
# count seen at init (interactive threadpools, dynamic spawn, ...).
#
# Safety: each OS thread holds exactly one entry and runs at most one task
# at a time (Julia's `:static` @threads scheduler). Not safe for nested
# @threads or @async on the same thread — acceptable given TurbulentTransport's
# usage patterns.
const _QLNN_THREAD_POOLS = AdaptiveArrayPool[]
const _QLNN_THREAD_POOLS_LOCK = ReentrantLock()

# One-shot capacity bump used by both __init__ and the lazy fallback.
function _qlnn_ensure_thread_pools_capacity!(n::Int)
    Base.@lock _QLNN_THREAD_POOLS_LOCK begin
        while length(_QLNN_THREAD_POOLS) < n
            push!(_QLNN_THREAD_POOLS, AdaptiveArrayPool())
        end
    end
    return _QLNN_THREAD_POOLS
end

# Public init hook: called from `TurbulentTransport.__init__()` so pools are
# pre-allocated for the runtime thread count. Idempotent — safe to call from
# tests / repl after a Threads.nthreads() bump.
#
# We size to `Threads.maxthreadid()`, not `Threads.nthreads()`: on Julia 1.9+
# `Threads.threadid()` can return any id in 1:maxthreadid() because of the
# interactive / foreign threadpools, but `Threads.nthreads()` only counts the
# default pool. Indexing thread-keyed buffers with the latter would BoundsError
# on tasks scheduled into the interactive pool.
function init_qlnn_thread_pools!()
    return _qlnn_ensure_thread_pools_capacity!(Threads.maxthreadid())
end

# Hot-path getter: read without lock when the pool is already big enough,
# fall back to the locked grow path only when first hitting a higher tid.
@inline function _qlnn_get_thread_pool()
    tid = Threads.threadid()
    @inbounds tid <= length(_QLNN_THREAD_POOLS) && return _QLNN_THREAD_POOLS[tid]
    _qlnn_ensure_thread_pools_capacity!(tid)
    return @inbounds _QLNN_THREAD_POOLS[tid]
end

# ---------------------------------------------------------------------------
#  Helpers for the threaded ensemble path
# ---------------------------------------------------------------------------
# In-place log10 + z-score normalization into a pre-allocated buffer `xx`.
# `x` is read-only (shared across threads); `xx` is the thread-local scratch.
function _qlnn_fill_normalized!(xx::AbstractMatrix{T}, model::QLNNmodel,
                                x::AbstractMatrix{T}) where {T<:Real}
    xnames = model.xnames
    xm = model.xm
    xσ = model.xσ
    nfeat, nsamp = size(x)
    @inbounds for i in 1:nfeat
        log_row = endswith(xnames[i], log_suffix)
        xm_i = T(xm[i])
        xσ_i = T(xσ[i])
        for j in 1:nsamp
            v = x[i, j]
            log_row && (v = log10(v))
            xx[i, j] = (v - xm_i) / xσ_i
        end
    end
    return xx
end

function _qlnn_fill_normalized!(xx::AbstractVector{T}, model::QLNNmodel,
                                x::AbstractVector{T}) where {T<:Real}
    xnames = model.xnames
    xm = model.xm
    xσ = model.xσ
    @inbounds for i in eachindex(x)
        v = x[i]
        endswith(xnames[i], log_suffix) && (v = log10(v))
        xx[i] = (v - T(xm[i])) / T(xσ[i])
    end
    return xx
end

# In-place denormalization: `out[i,j] = yn[i,j] * yσ[i] + ym[i]`.
# Safe for aliased `out === yn` (same position read then written).
@inline function _qlnn_denormalize_output!(out::AbstractMatrix{T}, yn::AbstractMatrix,
                                           model::QLNNmodel) where {T<:Real}
    ym = model.ym
    yσ = model.yσ
    @inbounds for i in axes(out, 1)
        ym_i = T(ym[i])
        yσ_i = T(yσ[i])
        for j in axes(out, 2)
            out[i, j] = T(yn[i, j]) * yσ_i + ym_i
        end
    end
    return out
end

@inline function _qlnn_denormalize_output!(out::AbstractVector{T}, yn::AbstractVector,
                                           model::QLNNmodel) where {T<:Real}
    ym = model.ym
    yσ = model.yσ
    @inbounds for i in eachindex(out)
        out[i] = T(yn[i]) * T(yσ[i]) + T(ym[i])
    end
    return out
end

# Forward pass through `model._pooled_chain` using the THREAD-indexed pool.
# `xx` is the already-normalized input (written by `_qlnn_fill_normalized!`).
# Result is written into `out` via `copyto!` + pool rewind, same as PooledChain
# but with a persistent thread-local pool instead of a fresh task-local one.
@inline function _qlnn_forward_threaded!(out::AbstractVecOrMat{T},
                                         model::QLNNmodel,
                                         xx::AbstractVecOrMat{T}) where {T<:Real}
    pool = _qlnn_get_thread_pool()
    checkpoint!(pool, T)
    copyto!(out, _forward_with_pool(model._pooled_chain.model, xx, pool))
    rewind!(pool, T)
    return out
end

# QLNN's pooled-chain forward pass is dominated by tiny matmuls (≤64-wide
# layers). With Julia on 1 thread and OpenBLAS spread across many cores, every
# GEMM call pays ~ms of thread-team setup/teardown for ~μs of FLOPs, so the
# default `BLAS.get_num_threads() == nproc` configuration can be 5–10× slower
# than `BLAS=1` even before considering ensemble-level parallelism. Surface a
# one-shot `@info` so users hitting this configuration know how to fix it.
const _QLNN_PERF_HINT_SHOWN = Ref(false)

function _qlnn_maybe_perf_hint()
    if _QLNN_PERF_HINT_SHOWN[]
        return nothing
    end
    if Threads.nthreads() == 1 && BLAS.get_num_threads() > 1
        @info """run_qlnn: tiny-matmul workload running with `Threads.nthreads()=1` \
and `BLAS.get_num_threads()=$(BLAS.get_num_threads())`. For QLNN's ≤64-wide \
layers this typically loses 5–10× to BLAS thread-team overhead. Recommended: \
restart Julia with `JULIA_NUM_THREADS=N` and call `BLAS.set_num_threads(1)` so \
the ensemble forward + per-radial integration parallelize at the Julia level. \
This message is shown once per session."""
        _QLNN_PERF_HINT_SHOWN[] = true
    end
    return nothing
end

# Ensemble: mean of physical predictions (regressors). For the stability
# classifier, callers should use `predict_unstable_prob` instead, which
# averages probabilities (post-σ) rather than logits.
#
# Threading strategy:
#   Threaded path  — use `_qlnn_forward_threaded!` which checkpoints the
#     thread-indexed `_QLNN_THREAD_POOLS[threadid()]` pool, runs
#     `_forward_with_pool`, then rewinds.  Thread-indexed pools are allocated
#     once at module load and reused across all `@threads` calls, so there are
#     zero fresh-allocation cold-starts and zero GC pressure in steady state.
#     PooledChain is NOT used here because it calls `get_task_local_pool()`,
#     which returns a fresh empty pool for every new task spawned by @threads.
#   Single-thread path — use `predict!` + PooledChain so the pool stays warm
#     across repeated calls (e.g., per-radial single-member inference).
function predict(ens::QLNNensemble, x::AbstractMatrix{T}) where {T<:Real}
    nmodels = length(ens.models)
    nouts = length(ens.models[1].ynames)
    nfeat, nsamp = size(x)
    if _qlnn_should_thread(nmodels)
        all_yy = Array{T,3}(undef, nouts, nsamp, nmodels)
        # One normalized-input scratch buffer per MEMBER (indexed by k, not
        # by threadid).  Each k appears exactly once in @threads, so there is
        # no data race and we allocate nmodels (≤20) matrices instead of
        # nthreads (≤256), avoiding ~12 MiB of unused scratch per head call.
        # The forward pass uses _qlnn_forward_threaded! which checkpoints the
        # thread-indexed AdaptiveArrayPool, so intermediate layer buffers are
        # reused across run_qlnn calls with zero GC pressure in steady state.
        input_bufs = [Matrix{T}(undef, nfeat, nsamp) for _ in 1:nmodels]
        Threads.@threads for k in 1:nmodels
            m = ens.models[k]
            xx = input_bufs[k]
            _qlnn_fill_normalized!(xx, m, x)
            slab = view(all_yy, :, :, k)
            _qlnn_forward_threaded!(slab, m, xx)
            _qlnn_denormalize_output!(slab, slab, m)
        end
        out = zeros(T, nouts, nsamp)
        @inbounds for k in 1:nmodels
            out .+= @view(all_yy[:, :, k])
        end
        out ./= nmodels
        return out
    else
        out = zeros(T, nouts, nsamp)
        each = Matrix{T}(undef, nouts, nsamp)
        for k in 1:nmodels
            predict!(each, ens.models[k], x)
            out .+= each
        end
        nmodels > 1 && (out ./= nmodels)
        return out
    end
end

function predict(ens::QLNNensemble, x::AbstractVector{T}) where {T<:Real}
    nmodels = length(ens.models)
    nouts = length(ens.models[1].ynames)
    nfeat = length(x)
    if _qlnn_should_thread(nmodels)
        all_yy = Matrix{T}(undef, nouts, nmodels)
        input_bufs = [Vector{T}(undef, nfeat) for _ in 1:nmodels]
        Threads.@threads for k in 1:nmodels
            m = ens.models[k]
            xx = input_bufs[k]
            _qlnn_fill_normalized!(xx, m, x)
            slab = view(all_yy, :, k)
            _qlnn_forward_threaded!(slab, m, xx)
            _qlnn_denormalize_output!(slab, slab, m)
        end
        out = zeros(T, nouts)
        @inbounds for k in 1:nmodels
            out .+= @view(all_yy[:, k])
        end
        out ./= nmodels
        return out
    else
        out = zeros(T, nouts)
        each = Vector{T}(undef, nouts)
        for k in 1:nmodels
            predict!(each, ens.models[k], x)
            out .+= each
        end
        nmodels > 1 && (out ./= nmodels)
        return out
    end
end

"""
    predict_unstable_prob(classifier, x) -> AbstractArray

For a `:stability`-target classifier, return `P(unstable | x) = σ(logit)`.
Single models: `σ(predict(model, x))`. Ensembles: `mean_k σ(predict_k(x))`
— calibration-correct mean-of-probabilities, matching the convention used at
training time. Uses the same pooled / threaded inference path as `predict`.
"""
function predict_unstable_prob(model::QLNNmodel, x::AbstractArray)
    model.target === :stability ||
        error("predict_unstable_prob: classifier.target must be :stability, got $(model.target)")
    y = predict(model, x)
    @inbounds for i in eachindex(y)
        y[i] = Flux.sigmoid(y[i])
    end
    return y
end

function predict_unstable_prob(ens::QLNNensemble, x::AbstractMatrix{T}) where {T<:Real}
    ens.models[1].target === :stability ||
        error("predict_unstable_prob: classifier.target must be :stability, got $(ens.models[1].target)")
    nmodels = length(ens.models)
    nouts = length(ens.models[1].ynames)
    nsamp = size(x, 2)
    nfeat = size(x, 1)
    if _qlnn_should_thread(nmodels)
        # Mean-of-probabilities: σ each member's logits FIRST, then average.
        # Differs from σ(mean-of-logits) near the decision boundary; matches
        # the convention used at training time.
        all_yy = Array{T,3}(undef, nouts, nsamp, nmodels)
        input_bufs = [Matrix{T}(undef, nfeat, nsamp) for _ in 1:nmodels]
        Threads.@threads for k in 1:nmodels
            m = ens.models[k]
            xx = input_bufs[k]
            _qlnn_fill_normalized!(xx, m, x)
            slab = view(all_yy, :, :, k)
            _qlnn_forward_threaded!(slab, m, xx)
            _qlnn_denormalize_output!(slab, slab, m)
            @inbounds for j in 1:nsamp, i in 1:nouts
                slab[i, j] = Flux.sigmoid(slab[i, j])
            end
        end
        out = zeros(T, nouts, nsamp)
        @inbounds for k in 1:nmodels
            out .+= @view(all_yy[:, :, k])
        end
        out ./= nmodels
        return out
    else
        out = zeros(T, nouts, nsamp)
        each = Matrix{T}(undef, nouts, nsamp)
        for k in 1:nmodels
            predict!(each, ens.models[k], x)
            @inbounds for j in 1:nsamp, i in 1:nouts
                out[i, j] += Flux.sigmoid(each[i, j])
            end
        end
        nmodels > 1 && (out ./= nmodels)
        return out
    end
end

function predict_unstable_prob(ens::QLNNensemble, x::AbstractVector{T}) where {T<:Real}
    ens.models[1].target === :stability ||
        error("predict_unstable_prob: classifier.target must be :stability, got $(ens.models[1].target)")
    nmodels = length(ens.models)
    nouts = length(ens.models[1].ynames)
    nfeat = length(x)
    if _qlnn_should_thread(nmodels)
        all_yy = Matrix{T}(undef, nouts, nmodels)
        input_bufs = [Vector{T}(undef, nfeat) for _ in 1:nmodels]
        Threads.@threads for k in 1:nmodels
            m = ens.models[k]
            xx = input_bufs[k]
            _qlnn_fill_normalized!(xx, m, x)
            slab = view(all_yy, :, k)
            _qlnn_forward_threaded!(slab, m, xx)
            _qlnn_denormalize_output!(slab, slab, m)
            @inbounds for i in eachindex(slab)
                slab[i] = Flux.sigmoid(slab[i])
            end
        end
        out = zeros(T, nouts)
        @inbounds for k in 1:nmodels
            out .+= @view(all_yy[:, k])
        end
        out ./= nmodels
        return out
    else
        out = zeros(T, nouts)
        each = Vector{T}(undef, nouts)
        for k in 1:nmodels
            predict!(each, ens.models[k], x)
            @inbounds for i in eachindex(each)
                out[i] += Flux.sigmoid(each[i])
            end
        end
        nmodels > 1 && (out ./= nmodels)
        return out
    end
end

#= ============================================================ =#
#  InputTJLF -> NN input matrix
#= ============================================================ =#

# Map the (possibly `_log10`-suffixed) xname to the `InputTJLF` field symbol.
# Per-species fields are flat scalars on `InputTGLF` (e.g. `RLTS_1`, `AS_2`,
# `VPAR_SHEAR_3`, ...) but on `InputTJLF` they live as length-`NS` vectors
# (`RLTS`, `AS`, `VPAR_SHEAR`, ...). When training xnames carry the flat
# `<base>_<i>` form, we have to undo the suffix and index into the vector.
# Other species-suffixed fields exist on `InputTGLF` only (`MASS_i`, `ZS_i`,
# ...) but the trained xnames pull from this fixed set.
const _QLNN_SPECIES_VECTOR_FIELDS = Set{Symbol}([
    :ZS, :AS, :MASS, :RLNS, :RLTS, :TAUS, :VPAR, :VPAR_SHEAR
])

# Trailing `_<int>` capture used to recognize species-tagged xnames.
# `^(.+)_(\d+)$` is intentionally non-greedy on the species index so multi-digit
# indices (NS > 9) still work, while the base captures everything before the
# last underscore (e.g. `VPAR_SHEAR_2` → ("VPAR_SHEAR", "2")).
const _QLNN_SPECIES_SUFFIX_RX = r"^(.+)_(\d+)$"

# Strip an optional `_log10` suffix from an xname (the actual log10 transform
# is applied later inside `_qlnn_apply_log10!`; this helper only returns the
# linear-space field name to look up on `InputTJLF`).
function _qlnn_strip_log_suffix(xname::AbstractString)
    return endswith(xname, log_suffix) ? xname[1:end-n_log_suffix] : String(xname)
end

# Resolve a single `xname` (after stripping any `_log10`) to a value on
# `InputTJLF`. Handles three cases:
#   1. plain field         (`BETAE`, `DEBYE`, `RMIN_LOC`, ...)
#   2. species-vector field (`AS_2` → `input_tjlf.AS[2]`, etc.)
#   3. anything else        → error with a precise diagnostic.
# `ky` is handled separately by the caller (it's a per-ky column, not a field).
function _qlnn_lookup_input_tjlf(input_tjlf::TJLF.InputTJLF, xname::AbstractString)
    base = _qlnn_strip_log_suffix(xname)
    sym = Symbol(base)
    if hasfield(typeof(input_tjlf), sym)
        return getfield(input_tjlf, sym)
    end
    m = match(_QLNN_SPECIES_SUFFIX_RX, base)
    if m !== nothing
        base_sym = Symbol(m.captures[1])
        if base_sym in _QLNN_SPECIES_VECTOR_FIELDS && hasfield(typeof(input_tjlf), base_sym)
            arr = getfield(input_tjlf, base_sym)
            arr === missing && return missing
            idx = parse(Int, m.captures[2])
            1 <= idx <= length(arr) ||
                error("QLNN input feature `$xname` indexes species $(idx) but " *
                      "InputTJLF.$(base_sym) only has length $(length(arr)) (NS=$(input_tjlf.NS)).")
            return arr[idx]
        end
    end
    error("QLNN input feature `$xname` does not map to any InputTJLF field " *
          "(tried `$sym`" *
          (m === nothing ? "" : " and `$(m.captures[1])`[$(m.captures[2])]") *
          ").")
end

# Fill an (n_features, nky) view at eltype `T` with QLNN inputs for one
# radial point. Mirrors what `run_tglfnn` does via the `T`-eltype pool.
# - per-ky scalar: only `ky` (lowercase). Populated from `ky_spectrum`.
# - everything else: pulled from `input_tjlf` via `_qlnn_lookup_input_tjlf`,
#   which handles species-vector fields (`AS_2` → `input_tjlf.AS[2]`).
#   `Bool`/`Int` widen to `T` cleanly (`Dual{T,V,N}(::Bool)` gives a
#   zero-partial Dual with value 0 or 1).
#
# Operates in place on `xs` so the batched `run_qlnn` path can pre-allocate
# ONE big `(n_features, sum_r nky_r)` matrix and fill per-radial chunks via
# views — no per-radial `Matrix{T}` allocation.
function _qlnn_fill_xs!(xs::AbstractMatrix{T}, input_tjlf::TJLF.InputTJLF{T},
                        ky_spectrum::AbstractVector,
                        xnames::Vector{String}) where {T<:Real}
    nky = size(xs, 2)
    @assert length(ky_spectrum) == nky "_qlnn_fill_xs!: nky mismatch (xs has $nky cols, ky_spectrum has $(length(ky_spectrum)))"
    @assert size(xs, 1) == length(xnames) "_qlnn_fill_xs!: nfeat mismatch (xs has $(size(xs,1)) rows, xnames has $(length(xnames)))"
    @inbounds for (i, name) in enumerate(xnames)
        if name == "ky"
            for j in 1:nky
                xs[i, j] = ky_spectrum[j]
            end
            continue
        end
        val = _qlnn_lookup_input_tjlf(input_tjlf, name)
        if val === missing
            error("QLNN input field `$name` is `missing` in InputTJLF — cannot build feature row.")
        end
        scalar = T(val)
        for j in 1:nky
            xs[i, j] = scalar
        end
    end
    return xs
end

# Allocating wrapper: returns a freshly allocated `(n_features, nky)` matrix
# at eltype `T`. Used by tests / single-radial debugging; the hot path goes
# through `_qlnn_fill_xs!` directly.
function _qlnn_build_xs(input_tjlf::TJLF.InputTJLF{T}, ky_spectrum::AbstractVector,
                       xnames::Vector{String}) where {T<:Real}
    xs = Matrix{T}(undef, length(xnames), length(ky_spectrum))
    return _qlnn_fill_xs!(xs, input_tjlf, ky_spectrum, xnames)
end

#= ============================================================ =#
#  QL tensor / gamma matrix packing
#= ============================================================ =#

# Parse a QL-weight ynames vector into per-row (field_idx, species_idx,
# strip_over_ky) tuples plus the (nf, ns) tensor shape it implies.
# field_set/species_set follow the canonical `(phi, apar [, bpar])` /
# `(e, D, C)` order from CANONICAL_FIELDS / CANONICAL_SPECIES so that the
# QL tensor we hand to `sum_ky_spectrum` is shaped exactly as TJLF expects.
function _qlnn_parse_qlweight_ynames(ynames::Vector{String})
    fields = String[]
    species = String[]
    over_ky = Bool[]
    for yn in ynames
        sp_str, fld_str, _rest = split(yn, "_"; limit=3)
        push!(species, String(sp_str))
        push!(fields, String(fld_str))
        push!(over_ky, endswith(yn, "_over_ky"))
    end
    field_set = String[]
    for f in fields
        f in field_set || push!(field_set, f)
    end
    # Validate against canonical sets BEFORE we sort/index. Catching a typo or
    # an unsupported species/field here yields a precise error; if we let the
    # `findfirst(==(s), species_set)` line below silently produce `nothing`,
    # the QL packing in `_qlnn_pack_qlweight!` would later throw a generic
    # `MethodError: no method matching setindex!(..., ::Nothing, ...)`.
    bad_fields = setdiff(field_set, _QLNN_FIELDS)
    if !isempty(bad_fields)
        error("QLNN bundle has unsupported field(s) $(collect(bad_fields)) in ynames; " *
              "expected a subset of $(_QLNN_FIELDS).")
    end
    bad_species = setdiff(unique(species), _QLNN_SPECIES)
    if !isempty(bad_species)
        error("QLNN bundle has unsupported species $(collect(bad_species)) in ynames; " *
              "expected a subset of $(_QLNN_SPECIES).")
    end
    # Pin field order to canonical (phi, apar, bpar) so SAT params match TJLF.
    sort!(field_set; by=f -> findfirst(==(f), _QLNN_FIELDS)::Int)
    species_set = collect(_QLNN_SPECIES)
    nf = length(field_set)
    ns = length(species_set)
    row_field = [findfirst(==(f), field_set)::Int for f in fields]
    row_species = [findfirst(==(s), species_set)::Int for s in species]
    return (; nf, ns, field_set, species_set, row_field, row_species, over_ky)
end

# Pack one regressor's `(n_outputs, nky)` predictions into the
# `(nf, ns, 1, nky, ntype)` QL tensor at the given `type_idx`. `mask`
# (length nky, nothing for "no gating") zeroes out predictions for stable
# rows so they don't drive the saturation rule.
#
# QL eltype is preserved (parameterized on `<:Real`) so ForwardDiff `Dual`
# tensors flow through unchanged. `mask`-zeroed entries use `zero(eltype(QL))`
# (a Dual with zero value AND zero partials) so AD never sees a discontinuous
# write.
function _qlnn_pack_qlweight!(QL::AbstractArray{<:Real,5}, y_pred::AbstractMatrix,
                              ynames_info, type_idx::Int,
                              ky_spectrum::AbstractVector,
                              mask::Union{Nothing,AbstractVector{Bool}})
    nf = ynames_info.nf
    ns = ynames_info.ns
    @assert size(QL, 1) >= nf
    @assert size(QL, 2) >= ns
    nky = size(QL, 4)
    @assert size(y_pred, 2) == nky "Predicted nky=$(size(y_pred, 2)) does not match QL tensor nky=$nky"
    Tz = zero(eltype(QL))
    @inbounds for r in axes(y_pred, 1)
        f = ynames_info.row_field[r]
        s = ynames_info.row_species[r]
        invert_ky = ynames_info.over_ky[r]
        for k in 1:nky
            v = y_pred[r, k]
            if invert_ky
                v = v * ky_spectrum[k]
            end
            if mask !== nothing && !mask[k]
                v = Tz
            end
            QL[f, s, 1, k, type_idx] = v
        end
    end
    return QL
end

#= ============================================================ =#
#  run_qlnn: end-to-end (NN forward -> sat-rule integration)
#= ============================================================ =#

"""
    run_qlnn(input_tjlfs::Vector{InputTJLF{T}}; bundle_name="QLNN", warn_nn_train_bounds=false,
             stability_threshold=0.5)

Run the QLNN bundle on a vector of `InputTJLF` (one per radial point) and
return a `Vector{GACODE.FluxSolution{T}}`.

The implementation mirrors `run_tglfnn`'s batched architecture: instead of a
per-radial-point NN forward pass, we concatenate every radial point's
`(n_features, nky_r)` block into a single `(n_features, sum_r nky_r)`
matrix and run each regressor / classifier exactly once over the whole
batch. Small-matrix GEMM is bandwidth-bound, so one forward pass on a
~600-column matrix is roughly an order of magnitude faster than `nr`
forward passes on ~24-column matrices.

Phases:

1. **Per-radial setup (serial).** Set `NMODES = 1` (QLNN predicts a single
   dominant eigenmode per ky), call `TJLF.get_sat_params(input_tjlf)` and
   `TJLF.get_ky_spectrum(input_tjlf, satParams.grad_r0)`, refresh
   `input_tjlf.KY_SPECTRUM`. `SAT_RULE` / `ALPHA_ZF` / `KYGRID_MODEL` /
   `UNITS` therefore drive the same `sum_ky_spectrum` a regular TJLF run
   would hit. Cheap relative to the NN passes; runs serially because each
   iteration mutates its own `InputTJLF`.

2. **Batched feature matrix.** Allocate one `Matrix{T}` of shape
   `(n_features, sum_r nky_r)` and fill per-radial chunks with
   `_qlnn_fill_xs!`. All non-`ky` rows are scalar within a chunk; only the
   `ky` row varies. Eltype is `T`, so `T = ForwardDiff.Dual` partials
   propagate end-to-end (matches the AD contract of the `:TJLF` /
   `:TGLFNN` / `:GKNN` paths).

3. **One NN forward pass per head.** `predict(bundle.energy, xs_all)` etc.
   each fire ONCE. The (optional) stability classifier also runs once on
   `xs_all`; its `P(unstable)` vector is sliced per radial point in
   phase 4 to build a hard gate (rows below `stability_threshold` are
   zeroed in the QL tensor and gamma matrix using `zero(T)` — Dual with
   zero value AND zero partials, so AD never sees a discontinuous write).

4. **Per-radial integration (threaded).** Slice the batched predictions
   back per radial point, pack the `(nf, ns, 1, nky, 5)` QL tensor (TJLF
   type axis 1=particle, 2=energy, 3=tor stress) and the `(1, nky)` gamma
   matrix at eltype `T`, run `TJLF.sum_ky_spectrum`, and reduce to a
   `GACODE.FluxSolution{T}` via the same `Qe / Qi / Γe / Γi / Πi`
   reductions as `run_tjlf`. `sum_ky_spectrum` is per-equilibrium physics
   so it's not batchable across radial points, but the work is independent
   per-`r` and `Threads.@threads` here doesn't compete with the (now
   finished) NN forward pass.

Two training-time `÷ky` scalings are inverted in phase 4 so what reaches
`sum_ky_spectrum` is the raw, un-normalized spectrum the saturation rule
expects:

- QL-weight side: per-output flag, encoded *in the yname suffix* of each
  row. Rows whose `yname` ends in `_over_ky` (typically the electron
  channels) were trained on `QL/ky`; we multiply those rows by `ky`
  row-by-row.
- Eigenvalue side: per-bundle flag, encoded as a single `Bool`
  `normalize_by_ky` saved into `eigenvalue_regressor.bson` alongside
  `xnames`/`ynames`/`xm`/`xσ`/.... `_qlnn_dict2mod` reads it back into
  `QLNNmodel.normalize_by_ky`. When `true`, the eigenvalue head was
  trained on `γ/ky` and we multiply each `γ_k` by `ky_spectrum[k]` before
  storing into `Γ`. When `false`, the output is already `γ` and we use it
  as-is. The choice is therefore fixed at training time and travels with
  the model — no runtime parameter to set.

When `warn_nn_train_bounds=true` we scan the batched feature matrix once
against `bundle.energy.xbounds` and emit one `@warn` per feature whose
worst-offender (across all radii × ky) falls outside the training range.
NaN/Inf values are hard errors. Bundles that didn't carry an `:xbounds`
field fall back to ±Inf bounds, so the check becomes a silent no-op.

The `stability_threshold` is exposed but not yet plumbed through `ActorTGLF`
(future actor parameter).
"""
function run_qlnn(input_tjlfs::Vector{TJLF.InputTJLF{T}};
                  bundle_name::AbstractString = "QLNN",
                  warn_nn_train_bounds::Bool = false,
                  stability_threshold::Real = 0.5) where {T<:Real}
    bundle = loadqlnnbundleonce(String(bundle_name))
    return run_qlnn(input_tjlfs, bundle;
                    warn_nn_train_bounds=warn_nn_train_bounds,
                    stability_threshold=stability_threshold)
end

function run_qlnn(input_tjlfs::Vector{TJLF.InputTJLF{T}}, bundle::QLNNbundle;
                  warn_nn_train_bounds::Bool = false,
                  stability_threshold::Real = 0.5) where {T<:Real}
    _qlnn_maybe_perf_hint()
    # Pre-parse ynames once per bundle (regressors share xnames, but ynames
    # may differ slightly between {energy, particle, momentum}).
    info_e = _qlnn_parse_qlweight_ynames(bundle.energy.ynames)
    info_p = _qlnn_parse_qlweight_ynames(bundle.particle.ynames)
    info_m = _qlnn_parse_qlweight_ynames(bundle.momentum.ynames)

    # All three regressors should agree on (nf, ns).
    if !(info_e.nf == info_p.nf == info_m.nf) || !(info_e.ns == info_p.ns == info_m.ns)
        error("QLNN bundle has inconsistent (nf, ns) across energy/particle/momentum regressors.")
    end
    # All bundle members must share the same xnames vector — the batched
    # feature matrix below is built once from `bundle.energy.xnames` and
    # reused for every regressor + classifier call.
    energy_xnames = bundle.energy.xnames
    for (label, m) in (("particle", bundle.particle), ("momentum", bundle.momentum),
                       ("eigenvalue", bundle.eigenvalue))
        m.xnames == energy_xnames ||
            error("QLNN bundle: $(label) xnames mismatch; expected energy regressor's xnames " *
                  "(length $(length(energy_xnames))), got length $(length(m.xnames)).")
    end
    if bundle.stability !== nothing
        bundle.stability.xnames == energy_xnames ||
            error("QLNN bundle: stability classifier xnames mismatch with energy regressor.")
    end
    nf = info_e.nf
    ns = info_e.ns

    nr = length(input_tjlfs)
    flux_solutions = Vector{GACODE.FluxSolution{T}}(undef, nr)
    if nr == 0
        return flux_solutions
    end

    # === Phase 1: per-radial setup =====================================
    # Set NMODES=1 (QLNN predicts one dominant mode per ky), pull
    # SaturationParameters + ky spectrum from TJLF, refresh KY_SPECTRUM.
    # Mutates each `input_tjlf` in place — safe to do serially because
    # each iteration touches only its own struct.
    #
    # Run the first iteration to probe TJLF's return types, then declare
    # strongly typed storage. Avoids the `Vector{Any}` boxing the previous
    # version paid in the per-radial scratch arrays.
    nky_r = Vector{Int}(undef, nr)
    it1 = input_tjlfs[1]
    it1.NMODES = 1
    sp1 = TJLF.get_sat_params(it1)
    ks1 = TJLF.get_ky_spectrum(it1, sp1.grad_r0)
    sat_params_v = Vector{typeof(sp1)}(undef, nr)
    ky_spectrums = Vector{typeof(ks1)}(undef, nr)
    sat_params_v[1] = sp1
    ky_spectrums[1] = ks1
    nky_r[1] = length(ks1)
    if it1.KY_SPECTRUM === missing || length(it1.KY_SPECTRUM) != nky_r[1]
        it1.KY_SPECTRUM = collect(ks1)
    else
        it1.KY_SPECTRUM .= ks1
    end
    for r in 2:nr
        it = input_tjlfs[r]
        # See note in `_qlnn_integrate_radial`: `intensity_sat` would otherwise
        # iterate `for k in 1:NMODES` (default `NS+2`) past the end of the
        # (nf, ns, 1, nky, 5) QL tensor we pack with `nmodes = 1`.
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

    # === Phase 2: build ONE batched feature matrix =====================
    # Each per-radial chunk shares all features except the `ky` row.
    # `_qlnn_fill_xs!` writes into a view, so this is one big alloc + nr
    # in-place fills (vs. nr separate `Matrix{T}` allocations before).
    nfeat = length(energy_xnames)
    xs_all = Matrix{T}(undef, nfeat, total_nky)
    for r in 1:nr
        c0 = chunk_starts[r]
        nk = nky_r[r]
        @views _qlnn_fill_xs!(xs_all[:, c0:c0+nk-1], input_tjlfs[r],
                              ky_spectrums[r], energy_xnames)
    end

    # Optional training-distribution check. All bundle members share xnames,
    # so we scan once against the energy regressor's bounds (worst-case across
    # all ky/radii so the user sees the most extreme offender per feature).
    if warn_nn_train_bounds
        _qlnn_check_xs_bounds(xs_all, bundle.energy)
    end

    # === Phase 3: ONE NN forward pass per head =========================
    # The big speedup vs. the per-radial version: 4 (or 5 with stability)
    # batched GEMMs total instead of `4 × nr`. `predict` preserves
    # `eltype(xs_all)`, so Duals flow through unchanged.
    Y_energy   = predict(bundle.energy,     xs_all)
    Y_particle = predict(bundle.particle,   xs_all)
    Y_momentum = predict(bundle.momentum,   xs_all)
    Y_eig      = predict(bundle.eigenvalue, xs_all)
    # Stability classifier output is `(1, total_nky)`; `vec(...)` flattens
    # so phase 4 can index `P_unstable[c0+j-1]` per radial chunk. Comparing
    # `Dual >= Real` reads the value part only, so the mask is `Bool` even
    # under AD (intentional — the hard gate is non-differentiable).
    P_unstable = bundle.stability === nothing ? nothing :
                 vec(predict_unstable_prob(bundle.stability, xs_all))

    # === Phase 4: per-radial QL packing + sum_ky_spectrum ==============
    # Independent across radial points and free of NN forward passes, so
    # threading here doesn't contend with BLAS threads inside `predict`.
    eig_norm_by_ky = bundle.eigenvalue.normalize_by_ky
    Threads.@threads for r in 1:nr
        c0 = chunk_starts[r]
        nk = nky_r[r]
        cols = c0:c0+nk-1
        ks = ky_spectrums[r]
        mask = nothing
        if P_unstable !== nothing
            mask = Bool[P_unstable[c0+j-1] >= stability_threshold for j in 1:nk]
        end
        flux_solutions[r] = _qlnn_integrate_radial(
            input_tjlfs[r], sat_params_v[r], ks,
            view(Y_energy,   :, cols),
            view(Y_particle, :, cols),
            view(Y_momentum, :, cols),
            view(Y_eig,      :, cols),
            mask, info_e, info_p, info_m,
            eig_norm_by_ky, nf, ns,
        )
    end
    return flux_solutions
end

# Per-radial-point integration: pack pre-computed predictions into the TJLF
# QL/Γ tensors and call `sum_ky_spectrum`. NN forward passes happen ONCE in
# `run_qlnn` over the batched feature matrix; this function only consumes
# those predictions, so it's pure tensor packing + saturation-rule physics.
#
# Every intermediate (QL tensor, Γ matrix, FluxSolution) is allocated at the
# InputTJLF's eltype `T`, so ForwardDiff `Dual`s propagate through the entire
# pipeline without ever round-tripping through Float64 (which would discard
# partials and silently zero the Jacobian).
function _qlnn_integrate_radial(input_tjlf::TJLF.InputTJLF{T},
                                sat_params,
                                ky_spectrum::AbstractVector,
                                y_energy::AbstractMatrix,
                                y_particle::AbstractMatrix,
                                y_momentum::AbstractMatrix,
                                y_eig::AbstractMatrix,
                                mask::Union{Nothing,AbstractVector{Bool}},
                                info_e, info_p, info_m,
                                eig_norm_by_ky::Bool,
                                nf::Int, ns::Int) where {T<:Real}
    nky = length(ky_spectrum)
    QL = zeros(T, nf, ns, 1, nky, 5)
    _qlnn_pack_qlweight!(QL, y_energy,   info_e, _QLNN_TARGET_TYPE_IDX[:energy],   ky_spectrum, mask)
    _qlnn_pack_qlweight!(QL, y_particle, info_p, _QLNN_TARGET_TYPE_IDX[:particle], ky_spectrum, mask)
    _qlnn_pack_qlweight!(QL, y_momentum, info_m, _QLNN_TARGET_TYPE_IDX[:momentum], ky_spectrum, mask)

    # Gamma matrix (1, nky). Multiply by ky if eigenvalue head was trained on γ/ky.
    Γ = zeros(T, 1, nky)
    Tz = zero(T)
    @inbounds for k in 1:nky
        γ = y_eig[1, k]
        if eig_norm_by_ky
            γ = γ * ky_spectrum[k]
        end
        if mask !== nothing && !mask[k]
            γ = Tz
        end
        Γ[1, k] = γ
    end

    # SAT_RULE / ALPHA_ZF / UNITS / RLNP_CUTOFF / ETG_FACTOR / ... all live
    # on `input_tjlf` and flow into sum_ky_spectrum. `sum_ky_spectrum` is
    # parameterized on `T<:Real`, so it accepts Duals.
    QL_flux_out, _flux_spectrum = TJLF.sum_ky_spectrum(input_tjlf, sat_params, Γ, QL)

    return GACODE.FluxSolution{T}(
        TJLF.Qe(QL_flux_out),
        TJLF.Qi(QL_flux_out),
        TJLF.Γe(QL_flux_out),
        TJLF.Γi(QL_flux_out),
        TJLF.Πi(QL_flux_out),
    )
end

# Single-input convenience wrapper (parallel to run_tglfnn).
function run_qlnn(input_tjlf::TJLF.InputTJLF{T}; kw...) where {T<:Real}
    return run_qlnn([input_tjlf]; kw...)[1]
end

# Convenience: convert InputTGLF -> InputTJLF first. We don't memoize this
# because the caller is typically `ActorTGLF._step`, which already maintains
# a Vector{InputTJLF} for `:QLNN` so width memory survives across iterations.
function run_qlnn(input_tglfs::Vector{TJLF.InputTGLF{T}}; kw...) where {T<:Real}
    input_tjlfs = TJLF.InputTJLF{T}[TJLF.InputTJLF{T}(it) for it in input_tglfs]
    return run_qlnn(input_tjlfs; kw...)
end

function run_qlnn(input_tglf::TJLF.InputTGLF{T}; kw...) where {T<:Real}
    return run_qlnn([input_tglf]; kw...)[1]
end
