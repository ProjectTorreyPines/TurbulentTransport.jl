#= ====================================== =#
#  Pooled Flux Layers using AdaptiveArrayPools
#= ====================================== =#
#
# Dynamic memory pooling for Flux layers - no fixed max_batch required.
# Uses AdaptiveArrayPools.jl for thread-safe, zero-allocation (after warmup) inference.
#
# Usage (simple - auto pool management):
#   model = PooledChain(poolify(flux_chain))
#   result = model(x)  # just call it!
#
# Usage (explicit - fine-grained control):
#   pooled = poolify(flux_chain)
#   y = @with_pool pool pooled(x)  # checkpoint → run → rewind
#
# Note: Pool memory is reused ACROSS @with_pool calls (after warmup),
# not within a single block. Each @with_pool rewinds on exit.

import Flux
import LinearAlgebra: mul!

#= ====================================== =#
#  Activation Detection
#= ====================================== =#

# Known activation functions that can be wrapped
const POOLABLE_ACTIVATIONS = (
    Flux.elu, Flux.relu, Flux.leakyrelu, Flux.selu, Flux.celu, Flux.gelu,
    Flux.sigmoid, Flux.hardsigmoid, Flux.hardtanh, Flux.tanh, Flux.softsign,
    Flux.softplus, Flux.swish, Flux.mish, Flux.lisht, Flux.tanhshrink
)

is_poolable_activation(f) = f in POOLABLE_ACTIVATIONS

#= ====================================== =#
#  PooledActivation
#= ====================================== =#

"""
    PooledActivation{F}

Wraps an activation function to use memory from AdaptiveArrayPools.
No pre-allocated buffer - acquires from global pool at call time.

Requires `@with_pool` block at the outermost call site.
"""
struct PooledActivation{F}
    σ::F
end

# Fallback callable (no pool) — enables use outside @with_pool / PooledChain context
@inline (pa::PooledActivation)(x) = pa.σ.(x)

#= ====================================== =#
#  PooledDense
#= ====================================== =#

"""
    PooledDense{D}

Wraps a `Flux.Dense` layer to use memory from AdaptiveArrayPools.
No pre-allocated buffer - acquires from global pool at call time.

Requires `@with_pool` block at the outermost call site.
"""
struct PooledDense{D<:Flux.Dense}
    dense::D
end

# Fallback callable (no pool) — enables use outside @with_pool / PooledChain context
@inline (pd::PooledDense)(x) = pd.dense(x)

#= ====================================== =#
#  PooledParallelAdd (ResNet skip connection)
#= ====================================== =#

"""
    PooledParallelAdd{T<:Tuple}

Pooled version of `Parallel(+, ...)` that uses in-place addition.

For ResNet-style skip connections: `Parallel(+, chain, identity)`
- Runs first branch, gets pooled output
- Adds remaining branches in-place
- Zero allocation from the `+` operation

Assumes all branches produce same-sized output.
"""
struct PooledParallelAdd{T<:Tuple}
    layers::T
end

@inline function (m::PooledParallelAdd)(x::AbstractVecOrMat)
    # Run first branch - this is our output buffer (already pooled)
    out = m.layers[1](x)

    # Add remaining branches in-place
    @inbounds for i in 2:length(m.layers)
        layer = m.layers[i]
        if layer === identity
            out .+= x
        else
            out .+= layer(x)
        end
    end

    return out
end

#= ====================================== =#
#  Pool-Propagating Forward Pass
#= ====================================== =#

"""
    _forward_with_pool(layer, x, pool)

Forward pass that propagates pool through layers explicitly.
Avoids repeated `get_task_local_pool()` calls by passing pool as argument.

This is an internal function used by `PooledChain` to optimize inference.
"""
# PooledDense - Matrix path
@inline function _forward_with_pool(pd::PooledDense, x::AbstractMatrix, pool)
    d = pd.dense
    Flux._size_check(d, x, 1 => size(d.weight, 2))
    xT = Flux._match_eltype(d, x)
    out = acquire_view!(pool, eltype(xT), size(d.weight, 1), size(xT, 2))
    mul!(out, d.weight, xT)
    return Flux.NNlib.bias_act!(d.σ, out, d.bias)
end

# PooledDense - Vector path (BLAS gemv)
@inline function _forward_with_pool(pd::PooledDense, x::AbstractVector, pool)
    d = pd.dense
    Flux._size_check(d, x, 1 => size(d.weight, 2))
    xT = Flux._match_eltype(d, x)
    out = acquire_view!(pool, eltype(xT), size(d.weight, 1))
    mul!(out, d.weight, xT)
    return Flux.NNlib.bias_act!(d.σ, out, d.bias)
end

# PooledActivation - Matrix path
@inline function _forward_with_pool(pa::PooledActivation, x::AbstractMatrix, pool)
    out = acquire_view!(pool, eltype(x), size(x))
    out .= pa.σ.(x)
    return out
end

# PooledActivation - Vector path
@inline function _forward_with_pool(pa::PooledActivation, x::AbstractVector, pool)
    out = acquire_view!(pool, eltype(x), length(x))
    out .= pa.σ.(x)
    return out
end

# PooledParallelAdd
@inline function _forward_with_pool(layer::PooledParallelAdd, x, pool)
    # Run first branch with pool
    out = _forward_with_pool(layer.layers[1], x, pool)

    # Add remaining branches in-place
    @inbounds for i in 2:length(layer.layers)
        branch = layer.layers[i]
        if branch === identity
            out .+= x
        else
            out .+= _forward_with_pool(branch, x, pool)
        end
    end
    return out
end

# Type-stable chain traversal using tuple recursion (like Flux._applychain)
@inline _chain_forward_with_pool(::Tuple{}, x, pool) = x
@inline _chain_forward_with_pool(layers::Tuple{Any}, x, pool) = _forward_with_pool(layers[1], x, pool)
@inline function _chain_forward_with_pool(layers::Tuple, x, pool)
    return _chain_forward_with_pool(Base.tail(layers), _forward_with_pool(layers[1], x, pool), pool)
end

@inline function _forward_with_pool(chain::Flux.Chain, x, pool)
    return _chain_forward_with_pool(chain.layers, x, pool)
end

# Fallback for identity and unknown layers (don't need pool)
@inline function _forward_with_pool(layer, x, pool)
    return layer(x)
end

#= ====================================== =#
#  Model Poolification
#= ====================================== =#

"""
    poolify(layer)

Recursively convert a Flux model to use pooled layers.

No `max_batch` is required - the pool dynamically handles any batch size.

# Example
```julia
using AdaptiveArrayPools

pooled_model = poolify(model.fluxmodel)

# Inference with any batch size
@with_pool pool begin
    y1 = pooled_model(x_batch_10)
    y2 = pooled_model(x_batch_100)   # works fine
    y3 = pooled_model(x_batch_1000)  # no error
end
```

# Notes
- Requires AdaptiveArrayPools.jl to be loaded
- Must wrap inference calls with `@with_pool` block
- Thread-safe via task-local storage
"""
function poolify(layer)
    if layer isa Flux.Dense
        return PooledDense(layer)
    elseif is_poolable_activation(layer)
        return PooledActivation(layer)
    elseif layer isa Flux.Chain
        return Flux.Chain(map(poolify, layer.layers)...)
    elseif layer isa Flux.Parallel
        pooled_branches = Tuple(map(poolify, layer.layers))
        # Use PooledParallelAdd for + connection (ResNet skip connections)
        if layer.connection === +
            return PooledParallelAdd(pooled_branches)
        else
            return Flux.Parallel(layer.connection, pooled_branches...)
        end
    else
        return layer
    end
end

# Note: poolify(model::TGLFNNmodel) is defined in tglf_nn.jl to avoid circular dependency

#= ====================================== =#
#  PooledChain (Auto-managed pool wrapper)
#= ====================================== =#

"""
    PooledChain{M<:Flux.Chain}

Wrapper that automatically manages `@with_pool` on each call.
Use this when you want a simple callable interface without explicit pool management.

# Usage
```julia
model = PooledChain(poolify(flux_chain))

# Allocating version - returns owned Array
y = model(x)

# In-place version - zero allocation (output first, Julia convention)
model(y, x)  # writes result to y
```

# Notes
- Allocating `model(x)`: copies result via `collect` (only allocation after warmup)
- In-place `model(out, x)`: uses `copyto!`, truly zero-allocation after warmup
- Thread-safe via task-local storage
- Pool intermediates are reused across calls
- Type parameter `M` preserves concrete Chain type for zero-allocation dispatch
"""
struct PooledChain{M<:Flux.Chain}
    model::M
end

# Allocating versions (return owned Array via collect)
# Uses _forward_with_pool to avoid repeated get_task_local_pool() calls
function (pm::PooledChain)(x::AbstractMatrix)
    T = eltype(x)
    pool = get_task_local_pool()
    checkpoint!(pool, T)
    result = collect(_forward_with_pool(pm.model, x, pool))
    rewind!(pool, T)
    return result
end

function (pm::PooledChain)(x::AbstractVector)
    T = eltype(x)
    pool = get_task_local_pool()
    checkpoint!(pool, T)
    result = collect(_forward_with_pool(pm.model, x, pool))
    rewind!(pool, T)
    return result
end

# In-place versions: pm(output, input) following Julia convention (mutated arg first)
# Uses _forward_with_pool to avoid repeated get_task_local_pool() calls
function (pm::PooledChain)(out::AbstractMatrix, x::AbstractMatrix)
    T = eltype(x)
    pool = get_task_local_pool()
    checkpoint!(pool, T)
    copyto!(out, _forward_with_pool(pm.model, x, pool))
    rewind!(pool, T)
    return out
end

function (pm::PooledChain)(out::AbstractVector, x::AbstractVector)
    T = eltype(x)
    pool = get_task_local_pool()
    checkpoint!(pool, T)
    copyto!(out, _forward_with_pool(pm.model, x, pool))
    rewind!(pool, T)
    return out
end