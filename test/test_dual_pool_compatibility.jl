#=
Minimal test demonstrating what AdaptiveArrayPools needs to support
for ForwardDiff.Dual number compatibility.

Currently, `acquire!` and `acquire_view!` only work with concrete float types
(Float64, Float32). To propagate ForwardDiff.Dual numbers through pooled
neural network inference, the pool needs to handle Dual element types.

Run standalone:
    julia --project -e 'using TurbulentTransport, Test; include("test/test_dual_pool_compatibility.jl")'
=#

using Test
using ForwardDiff
using ForwardDiff: Dual, Partials
using AdaptiveArrayPools: AdaptiveArrayPool, acquire!, acquire_view!, @with_pool, checkpoint!, rewind!

@testset "AdaptiveArrayPools Dual compatibility" begin

    # ForwardDiff.Dual is the element type used for automatic differentiation.
    # A Dual{T,V,N} carries a value of type V plus N partial derivatives.
    # For our use case: T=ForwardDiff.Tag, V=Float64, N varies (typically 1-40).
    D = Dual{Nothing,Float64,3}  # 3 partials — representative

    @testset "acquire! with Dual element type" begin
        # This is the core operation: acquire a buffer of Dual numbers from the pool.
        # Currently fails because the pool only handles concrete float types.
        @with_pool pool begin
            # Vector
            v = acquire!(pool, D, 10)
            @test v isa Vector{D}
            @test length(v) == 10

            # Matrix
            m = acquire!(pool, D, 5, 3)
            @test m isa Matrix{D}
            @test size(m) == (5, 3)

            # 3D array (used in ensemble inference)
            a = acquire!(pool, D, 4, 3, 2)
            @test a isa Array{D,3}
            @test size(a) == (4, 3, 2)
        end
    end

    @testset "acquire_view! with Dual element type" begin
        # acquire_view! returns a view into pooled memory (zero-copy).
        # Used in _forward_with_pool for intermediate layer outputs.
        @with_pool pool begin
            v = acquire_view!(pool, D, 10)
            @test v isa AbstractVector{D}
            @test length(v) == 10

            m = acquire_view!(pool, D, 5, 3)
            @test m isa AbstractMatrix{D}
            @test size(m) == (5, 3)
        end
    end

    @testset "Dual buffers are writable and preserve partials" begin
        # The acquired buffers must support in-place operations with Dual arithmetic.
        # This is what happens inside PooledDense: mul!(out, weight, x) then bias_act!
        @with_pool pool begin
            out = acquire!(pool, D, 3)

            # Write Dual values with known partials
            out[1] = Dual{Nothing}(1.0, Partials((2.0, 3.0, 4.0)))
            out[2] = Dual{Nothing}(5.0, Partials((6.0, 7.0, 8.0)))
            out[3] = Dual{Nothing}(9.0, Partials((10.0, 11.0, 12.0)))

            # Verify values and partials survive
            @test ForwardDiff.value(out[1]) == 1.0
            @test ForwardDiff.partials(out[1]) == Partials((2.0, 3.0, 4.0))

            # In-place broadcast (activation function application)
            out .= tanh.(out)
            @test ForwardDiff.value(out[1]) ≈ tanh(1.0)
        end
    end

    @testset "mul! with Dual arrays (PooledDense forward pass)" begin
        # This is the actual operation in _forward_with_pool(::PooledDense, ...):
        #   out = acquire_view!(pool, T, nout)
        #   mul!(out, weight, x)
        # where weight is Float64 but x and out are Dual.
        using LinearAlgebra: mul!

        W = randn(3, 4)  # Float64 weight matrix (frozen model weights)
        b = randn(3)     # Float64 bias

        # Input is Dual (carrying derivative information)
        x_dual = [Dual{Nothing}(randn(), Partials(Tuple(randn(3)))) for _ in 1:4]

        @with_pool pool begin
            out = acquire_view!(pool, D, 3)
            mul!(out, W, x_dual)
            out .+= b  # bias addition (broadcasts Float64 + Dual → Dual)
            out .= tanh.(out)  # activation

            # Verify result matches non-pooled computation
            y_ref = tanh.(W * x_dual .+ b)
            @test all(out .≈ y_ref)
        end
    end

    @testset "checkpoint/rewind cycle with Dual buffers" begin
        # PooledChain does checkpoint! → forward pass → rewind! on every call.
        # Dual buffers must survive within a checkpoint and be reclaimable after rewind.
        @with_pool pool begin
            checkpoint!(pool)

            v1 = acquire!(pool, D, 10)
            v1 .= Dual{Nothing}(1.0, Partials((1.0, 0.0, 0.0)))
            @test ForwardDiff.value(v1[1]) == 1.0

            v2 = acquire!(pool, D, 5, 3)
            @test size(v2) == (5, 3)

            rewind!(pool)

            # After rewind, pool memory is reclaimed — new acquire reuses it
            v3 = acquire!(pool, D, 10)
            @test length(v3) == 10
        end
    end

    @testset "End-to-end: PooledChain with Dual input" begin
        # Integration test: run a full PooledChain forward pass with Dual inputs.
        # This is the actual use case — ForwardDiff pushes Dual numbers through
        # the neural network to compute gradients of the transport fluxes.
        using TurbulentTransport: PooledChain, poolify, loadmodel, generate_valid_input

        ensemble = loadmodel("sat2_em_d3d_azf-1")
        model = ensemble.models[1]
        pm = PooledChain(poolify(model.fluxmodel))

        # Create Dual input (simulating ForwardDiff jacobian computation)
        x_float = generate_valid_input(model)
        N = length(x_float)
        x_dual = [Dual{Nothing}(x_float[i], Partials(Tuple(Float64(j == i) for j in 1:N))) for i in 1:N]

        # This should work through the pool natively (AdaptiveArrayPools supports Dual)
        y_dual = pm(x_dual)

        # Extract primal values — must match Float64 result
        y_float = pm(x_float)
        y_values = ForwardDiff.value.(y_dual)
        @test y_values ≈ y_float

        # Partials must be non-trivial (the NN is not constant)
        J = stack(ForwardDiff.partials.(y_dual))  # (N_partials, N_outputs) Jacobian
        @test !all(iszero, J)
    end
end
