#= ====================================== =#
#  Pooled Flux Layers using AdaptiveArrayPools
#= ====================================== =#
#
# Dynamic memory pooling for Flux layers - no fixed max_batch required.
# Uses AdaptiveArrayPools.jl for thread-safe, zero-allocation (after warmup) inference.
#
# Usage:
#   pooled_model = poolify(flux_chain)
#   @with_pool pool begin
#       result = pooled_model(x)  # any batch size works
#   end

import Flux
import LinearAlgebra: mul!
using AdaptiveArrayPools: get_global_pool, acquire!

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
    pool = get_global_pool()
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
    pool = get_global_pool()
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

"""
    poolify(model::TGLFNNmodel)

Convert a TGLFNNmodel's fluxmodel to use pooled layers.
"""
function poolify(model::TGLFNNmodel)
    poolify(model.fluxmodel)
end
