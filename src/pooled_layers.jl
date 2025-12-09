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

@inline function _pooled_activation_forward!(pa::PooledActivation, x::AbstractVecOrMat)
    pool = get_task_local_pool()
    out = acquire!(pool, Float64, size(x))
    out .= pa.σ.(x)
    return out
end

# Matrix input → Matrix output
(pa::PooledActivation)(x::AbstractMatrix) = _pooled_activation_forward!(pa, x)

# Vector input → Vector output
(pa::PooledActivation)(x::AbstractVector) = vec(_pooled_activation_forward!(pa, x))

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

@inline function _pooled_dense_forward!(pd::PooledDense, x::AbstractVecOrMat)
    pool = get_task_local_pool()
    d = pd.dense
    Flux._size_check(d, x, 1 => size(d.weight, 2))
    xT = Flux._match_eltype(d, x)
    out = acquire!(pool, Float64, size(d.weight, 1), size(xT, 2))  # Vector output
    mul!(out, d.weight, xT)
    return Flux.NNlib.bias_act!(d.σ, out, d.bias)
end

# Matrix input → Matrix output
(pd::PooledDense)(x::AbstractMatrix) = _pooled_dense_forward!(pd, x)

# Vector input → Vector output
(pd::PooledDense)(x::AbstractVector) = vec(_pooled_dense_forward!(pd, x))

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
@with_pool function (pm::PooledChain)(x::AbstractMatrix)
    return collect(pm.model(x))::Matrix{Float64}
end

@with_pool function (pm::PooledChain)(x::AbstractVector)
    return collect(pm.model(x))::Vector{Float64}
end

# In-place versions: pm(output, input) following Julia convention (mutated arg first)
@with_pool function (pm::PooledChain)(out::AbstractMatrix, x::AbstractMatrix)
    copyto!(out, pm.model(x))
    return out
end

@with_pool function (pm::PooledChain)(out::AbstractVector, x::AbstractVector)
    copyto!(out, pm.model(x))
    return out
end