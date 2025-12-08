#= ====================================== =#
#  Buffered Flux Layers for Zero-Allocation Inference
#= ====================================== =#
#
# These wrappers preallocate output buffers to eliminate allocations
# during forward passes. Useful for hot loops and real-time inference.
#
# Usage:
#   buffered_model = bufferize(flux_chain, max_batch_size)
#   y = buffered_model(x)  # allocation-free for batch_size <= max_batch_size

import Flux
import LinearAlgebra: mul!

#= ====================================== =#
#  Buffered Activation
#= ====================================== =#

"""
    BufferedActivation{F, B}

Wraps an activation function with a preallocated output buffer.
Eliminates allocation from `σ.(x)` broadcasting.
"""
struct BufferedActivation{F, B<:AbstractMatrix}
    σ::F
    buffer::B
end

function (ba::BufferedActivation)(x::AbstractVecOrMat)
    out_rows = size(x, 1)
    out_cols = ndims(x) == 1 ? 1 : size(x, 2)
    out = view(ba.buffer, 1:out_rows, 1:out_cols)
    out .= ba.σ.(x)
    return out
end

# Known activation functions that can be wrapped
const BUFFERABLE_ACTIVATIONS = (
    Flux.elu, Flux.relu, Flux.leakyrelu, Flux.selu, Flux.celu, Flux.gelu,
    Flux.sigmoid, Flux.hardsigmoid, Flux.hardtanh, Flux.tanh, Flux.softsign,
    Flux.softplus, Flux.swish, Flux.mish, Flux.lisht, Flux.tanhshrink
)

is_bufferable_activation(f) = f in BUFFERABLE_ACTIVATIONS

#= ====================================== =#
#  Buffered Dense
#= ====================================== =#

"""
    BufferedDense{D, B}

Wraps a `Flux.Dense` layer with a preallocated output buffer.
Uses in-place `mul!` for matrix multiplication.
"""
struct BufferedDense{D<:Flux.Dense, B<:AbstractMatrix}
    dense::D
    buffer::B
end

function (bd::BufferedDense)(x::AbstractVecOrMat)
    d = bd.dense
    Flux._size_check(d, x, 1 => size(d.weight, 2))
    xT = Flux._match_eltype(d, x)

    out_rows = size(d.weight, 1)
    out_cols = ndims(xT) == 1 ? 1 : size(xT, 2)
    out = view(bd.buffer, 1:out_rows, 1:out_cols)

    mul!(out, d.weight, xT)
    return Flux.NNlib.bias_act!(d.σ, out, d.bias)
end

#= ====================================== =#
#  Model Bufferization
#= ====================================== =#

"""
    max_layer_width(layer) -> Int

Scan a Flux model to find the maximum output dimension (width) of any layer.
Used to size activation buffers appropriately to handle the widest layer.
"""
function max_layer_width(layer)
    if layer isa Flux.Dense
        return size(layer.weight, 1)
    elseif layer isa Flux.Chain
        return maximum(max_layer_width(l) for l in layer.layers; init=0)
    elseif layer isa Flux.Parallel
        # Parallel with vcat: sum of all branch outputs
        return sum(max_layer_width(l) for l in layer.layers)
    else
        return 0
    end
end

"""
    _bufferize(layer, max_batch, max_features) -> layer

Recursively replace Flux layers with buffered versions.
Internal function - use `bufferize` instead.
"""
function _bufferize(layer, max_batch::Int, max_features::Int)
    if layer isa Flux.Dense
        buf = Matrix{Float64}(undef, size(layer.weight, 1), max_batch)
        return BufferedDense(layer, buf)
    elseif is_bufferable_activation(layer)
        buf = Matrix{Float64}(undef, max_features, max_batch)
        return BufferedActivation(layer, buf)
    elseif layer isa Flux.Chain
        return Flux.Chain(map(l -> _bufferize(l, max_batch, max_features), layer.layers)...)
    elseif layer isa Flux.Parallel
        return Flux.Parallel(layer.connection, map(l -> _bufferize(l, max_batch, max_features), layer.layers)...)
    else
        return layer
    end
end

"""
    bufferize(model, max_batch::Int) -> model

Convert a Flux model to use preallocated buffers for inference.

# Arguments
- `model`: A Flux model (Chain, Dense, etc.)
- `max_batch`: Maximum expected batch size

# Returns
A new model with `Dense` and activation layers wrapped in buffered versions.
Forward passes will be allocation-free for batch sizes ≤ `max_batch`.

# Example
```julia
model = Chain(Dense(10 => 32, relu), Dense(32 => 5))
buffered = bufferize(model, 100)
y = buffered(x)  # allocation-free for size(x, 2) <= 100
```

# Notes
- Original model is not modified
- Buffered model is NOT thread-safe (shared buffers)
- For multi-threaded use, create separate buffered models per thread
"""
function bufferize(model, max_batch::Int)
    max_features = max_layer_width(model)
    return _bufferize(model, max_batch, max_features)
end
