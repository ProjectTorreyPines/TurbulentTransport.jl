# Test PooledChain with AdaptiveArrayPools
using TurbulentTransport: poolify, PooledChain, PooledDense, PooledActivation, PooledParallelAdd

@testset "PooledChain" begin
    # Sample models to test
    TEST_MODELS = [
        "sat0_em_d3d",
        "sat2_em_d3d_azf-1",
        "sat3_em_d3d_azf-1",
    ]

    @testset "PooledChain correctness" begin
        @testset "Model: $model_name" for model_name in TEST_MODELS
            ensemble = loadmodel(model_name)
            model = ensemble.models[1]

            # Create PooledChain (auto-managed pool)
            pm = PooledChain(poolify(model.fluxmodel))

            # Generate valid test input
            x_vec = generate_valid_input(model)
            x_mat = generate_valid_input_matrix(model, 10)

            @testset "Vector input" begin
                y_orig = model.fluxmodel(x_vec)
                y_pooled = pm(x_vec)  # No @with_pool needed!

                @test y_orig isa AbstractVector
                @test y_pooled isa AbstractVector
                @test size(y_orig) == size(y_pooled)
                @test isapprox(y_orig, y_pooled; rtol=REGRESSION_RTOL)
            end

            @testset "Matrix input (batch)" begin
                y_orig = model.fluxmodel(x_mat)
                y_pooled = pm(x_mat)

                @test size(y_orig) == size(y_pooled)
                @test isapprox(y_orig, y_pooled; rtol=REGRESSION_RTOL)
            end
        end
    end

    @testset "PooledChain with varying batch sizes" begin
        ensemble = loadmodel("sat2_em_d3d_azf-1")
        model = ensemble.models[1]
        pm = PooledChain(poolify(model.fluxmodel))

        # No max_batch limit - works with any size
        @testset "Batch size: $batch_size" for batch_size in [1, 5, 10, 50, 100, 500]
            x = generate_valid_input_matrix(model, batch_size)

            y_orig = model.fluxmodel(x)
            y_pooled = pm(x)

            @test size(y_orig) == size(y_pooled)
            @test isapprox(y_orig, y_pooled; rtol=REGRESSION_RTOL)
        end
    end

    @testset "PooledChain multiple sequential calls" begin
        ensemble = loadmodel("sat2_em_d3d_azf-1")
        model = ensemble.models[1]
        pm = PooledChain(poolify(model.fluxmodel))

        x1 = generate_valid_input_matrix(model, 10)
        x2 = generate_valid_input_matrix(model, 20)
        x3 = generate_valid_input_matrix(model, 5)

        # Each call auto-manages pool (checkpoint → run → rewind)
        # Pool memory reused across calls
        r1 = pm(x1)
        r2 = pm(x2)
        r3 = pm(x3)

        @test size(r1) == (length(model.ynames), 10)
        @test size(r2) == (length(model.ynames), 20)
        @test size(r3) == (length(model.ynames), 5)

        # Verify correctness
        @test isapprox(r1, model.fluxmodel(x1); rtol=REGRESSION_RTOL)
        @test isapprox(r2, model.fluxmodel(x2); rtol=REGRESSION_RTOL)
        @test isapprox(r3, model.fluxmodel(x3); rtol=REGRESSION_RTOL)
    end

    @testset "PooledChain convenience constructor" begin
        ensemble = loadmodel("sat2_em_d3d_azf-1")
        model = ensemble.models[1]

        # Direct from TGLFNNmodel
        pm = PooledChain(model)

        x = generate_valid_input_matrix(model, 10)
        y_orig = model.fluxmodel(x)
        y_pooled = pm(x)

        @test isapprox(y_orig, y_pooled; rtol=REGRESSION_RTOL)
    end

    @testset "Bit-wise identical accuracy" begin
        ensemble = loadmodel("sat2_em_d3d_azf-1")
        model = ensemble.models[1]

        # Poolify the model's fluxmodel
        pm = PooledChain(model.fluxmodel)

        x = generate_valid_input_matrix(model, 10)

        # Test Matrix path (where input x is Matrix)
        @test model.fluxmodel(x) == pm(x) # bit-wise identical

        # Test vector path (where input x is Vector)
        x_vec = x[:, 1]
        @test model.fluxmodel(x_vec) == pm(x_vec) # bit-wise identical
    end

    @testset "PooledChain minimal allocation (output only)" begin
        # After warmup, only the output array is allocated (via collect).
        # Pool intermediates are reused — no GC pressure from intermediate layers.
        # Julia's array allocation costs ~2 allocations + output_size bytes + overhead.
        ensemble = loadmodel("sat2_em_d3d_azf-1")
        model = ensemble.models[1]
        pm = PooledChain(model)
        nouts = length(model.ynames)

        # Test with various batch sizes
        @testset "Batch size: $batch_size" for batch_size in [1, 10, 50]
            x = generate_valid_input_matrix(model, batch_size)

            # Warmup call - this allocates pool memory
            pm(x)

            # Subsequent calls allocate only output array (nouts × batch_size Float64s)
            # Plus fixed overhead (~128 bytes for array header)
            allocs = @allocated pm(x)
            expected_output_bytes = sizeof(Float64) * nouts * batch_size + 128
            @test allocs <= expected_output_bytes
        end

        # Test vector input
        @testset "Vector input" begin
            x_vec = generate_valid_input(model)

            # Warmup
            pm(x_vec)

            # Should allocate only output vector (nouts Float64s) + overhead
            allocs = @allocated pm(x_vec)
            expected_output_bytes = sizeof(Float64) * nouts + 128
            @test allocs <= expected_output_bytes
        end
    end

    @testset "PooledChain in-place: pm(out, x) zero allocation" begin
        ensemble = loadmodel("sat2_em_d3d_azf-1")
        model = ensemble.models[1]
        pm = PooledChain(model)
        nouts = length(model.ynames)
        nins = length(model.xnames)

        # Matrix version
        @testset "Matrix in-place" begin
            batch_size = 10
            x = generate_valid_input_matrix(model, batch_size)
            out = Matrix{Float64}(undef, nouts, batch_size)

            # Warmup
            pm(out, x)

            # Verify correctness
            y_alloc = pm(x)
            @test isapprox(out, y_alloc; rtol=REGRESSION_RTOL)

            # Zero allocation after warmup
            allocs = @allocated pm(out, x)
            @test allocs == 0
        end

        # Vector version
        @testset "Vector in-place" begin
            x_vec = generate_valid_input(model)
            out_vec = Vector{Float64}(undef, nouts)

            # Warmup
            pm(out_vec, x_vec)

            # Verify correctness
            y_alloc = pm(x_vec)
            @test isapprox(out_vec, y_alloc; rtol=REGRESSION_RTOL)

            # Zero allocation after warmup
            allocs = @allocated pm(out_vec, x_vec)
            @test allocs == 0
        end

        # View as output (Matrix)
        @testset "Matrix view as output" begin
            batch_size = 10
            x = generate_valid_input_matrix(model, batch_size)

            # Output view into larger buffer
            buffer = Matrix{Float64}(undef, nouts + 2, batch_size + 2)
            out_view = @view buffer[2:nouts+1, 2:batch_size+1]

            # Warmup
            pm(out_view, x)

            # Verify correctness
            y_alloc = pm(x)
            @test isapprox(out_view, y_alloc; rtol=REGRESSION_RTOL)

            # Zero allocation after warmup
            allocs = @allocated pm(out_view, x)
            @test allocs == 0
        end

        # View as output (Vector)
        @testset "Vector view as output" begin
            x_vec = generate_valid_input(model)

            # Output view into larger buffer
            buffer = Vector{Float64}(undef, nouts + 2)
            out_view = @view buffer[2:nouts+1]

            # Warmup
            pm(out_view, x_vec)

            # Verify correctness
            y_alloc = pm(x_vec)
            @test isapprox(out_view, y_alloc; rtol=REGRESSION_RTOL)

            # Zero allocation after warmup
            allocs = @allocated pm(out_view, x_vec)
            @test allocs == 0
        end

        # View as input (Matrix)
        @testset "Matrix view as input" begin
            batch_size = 10

            # Input view from larger buffer
            x_buffer = Matrix{Float64}(undef, nins + 2, batch_size + 2)
            x_core = generate_valid_input_matrix(model, batch_size)
            x_buffer[2:nins+1, 2:batch_size+1] .= x_core
            x_view = @view x_buffer[2:nins+1, 2:batch_size+1]

            out = Matrix{Float64}(undef, nouts, batch_size)

            # Warmup
            pm(out, x_view)

            # Verify correctness (compare with non-view input)
            y_expected = pm(x_core)
            @test isapprox(out, y_expected; rtol=REGRESSION_RTOL)

            # Zero allocation after warmup
            allocs = @allocated pm(out, x_view)
            @test allocs == 0
        end

        # View to view (Matrix)
        @testset "View to view (Matrix)" begin
            batch_size = 10

            # Input view
            x_buffer = Matrix{Float64}(undef, nins + 2, batch_size + 2)
            x_core = generate_valid_input_matrix(model, batch_size)
            x_buffer[2:nins+1, 2:batch_size+1] .= x_core
            x_view = @view x_buffer[2:nins+1, 2:batch_size+1]

            # Output view
            out_buffer = Matrix{Float64}(undef, nouts + 2, batch_size + 2)
            out_view = @view out_buffer[2:nouts+1, 2:batch_size+1]

            # Warmup
            pm(out_view, x_view)

            # Verify correctness
            y_expected = pm(x_core)
            @test isapprox(out_view, y_expected; rtol=REGRESSION_RTOL)

            # Zero allocation after warmup
            allocs = @allocated pm(out_view, x_view)
            @test allocs == 0
        end

        # View to view (Vector)
        @testset "View to view (Vector)" begin
            # Input view
            x_buffer = Vector{Float64}(undef, nins + 2)
            x_core = generate_valid_input(model)
            x_buffer[2:nins+1] .= x_core
            x_view = @view x_buffer[2:nins+1]

            # Output view
            out_buffer = Vector{Float64}(undef, nouts + 2)
            out_view = @view out_buffer[2:nouts+1]

            # Warmup
            pm(out_view, x_view)

            # Verify correctness
            y_expected = pm(x_core)
            @test isapprox(out_view, y_expected; rtol=REGRESSION_RTOL)

            # Zero allocation after warmup
            allocs = @allocated pm(out_view, x_view)
            @test allocs == 0
        end
    end

    @testset "PooledDense and PooledActivation types" begin
        ensemble = loadmodel("sat2_em_d3d_azf-1")
        model = ensemble.models[1]
        pooled = poolify(model.fluxmodel)

        # Check that PooledDense layers exist
        first_layer = pooled.layers[1]
        @test first_layer isa PooledDense

        # Check that the wrapped dense is preserved
        @test first_layer.dense isa Flux.Dense
    end

    @testset "PooledParallelAdd for ResNet blocks" begin
        ensemble = loadmodel("sat2_em_d3d_azf-1")
        model = ensemble.models[1]
        pm = PooledChain(poolify(model.fluxmodel))

        # Count PooledParallelAdd instances (should be 5 for this model)
        function count_parallel_add(layer, count=Ref(0))
            if layer isa PooledParallelAdd
                count[] += 1
            elseif layer isa Flux.Chain
                for l in layer.layers
                    count_parallel_add(l, count)
                end
            end
            return count[]
        end

        n_parallel = count_parallel_add(pm.model)
        @test n_parallel == 5  # 5 ResNet blocks

        # Verify PooledParallelAdd works correctly
        x = generate_valid_input_matrix(model, 1)
        y_orig = model.fluxmodel(x)
        y_pooled = pm(x)
        @test isapprox(y_orig, y_pooled; rtol=REGRESSION_RTOL)
    end

    @testset "PooledChain all available models" begin
        all_models = available_models()

        # Filter to TGLF-NN .bson files only (skip FINN/ModeID models, directories/symlinks)
        bson_models = filter(m -> !contains(m, "/") && !startswith(m, "finn_") && !startswith(m, "modeid_") && (endswith(m, ".bson") || !contains(m, ".")), all_models)

        # Test a subset for CI speed (full test can be run manually)
        test_subset = first(bson_models, 10)

        @testset "Model: $model_name" for model_name in test_subset
            ensemble = try
                loadmodel(model_name)
            catch e
                @warn "Failed to load model: $model_name" exception=e
                continue
            end

            model = ensemble.models[1]
            pm = PooledChain(model)

            x = generate_valid_input(model)
            y_orig = model.fluxmodel(x)
            y_pooled = pm(x)

            @test isapprox(y_orig, y_pooled; rtol=REGRESSION_RTOL)
        end
    end

    @testset "Thread safety (task-local pools)" begin
        # PooledChain uses task-local storage, so concurrent access should be safe
        # Each thread gets its own pool → no data races
        if Threads.nthreads() > 1
            ensemble = loadmodel("sat2_em_d3d_azf-1")
            model = ensemble.models[1]
            pm = PooledChain(model)
            x = generate_valid_input(model)

            # Warmup
            pm(x)

            # Concurrent access from multiple threads
            nthreads = Threads.nthreads()
            results = Vector{Vector{Float64}}(undef, nthreads)

            Threads.@threads for i in 1:nthreads
                results[i] = pm(x)
            end

            # All threads should produce identical results
            @test all(isapprox(r, results[1]; rtol=REGRESSION_RTOL) for r in results)
        else
            @info "Skipping thread safety test (single-threaded Julia)"
        end
    end
end
