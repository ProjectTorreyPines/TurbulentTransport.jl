# Test poolify with AdaptiveArrayPools
using TurbulentTransport.AdaptiveArrayPools: @with_pool
using TurbulentTransport: poolify, PooledDense, PooledActivation, PooledParallelAdd

@testset "Poolify" begin
        # Sample models to test
        TEST_MODELS = [
            "sat0_em_d3d",
            "sat2_em_d3d_azf-1",
            "sat3_em_d3d_azf-1",
        ]

        @testset "poolify correctness" begin
            @testset "Model: $model_name" for model_name in TEST_MODELS
                ensemble = loadmodel(model_name)
                model = ensemble.models[1]

                # Poolify the fluxmodel
                pooled = poolify(model.fluxmodel)

                # Generate valid test input
                x_vec = generate_valid_input(model)
                x_mat = generate_valid_input_matrix(model, 10)

                @testset "Vector input" begin
                    y_orig = model.fluxmodel(x_vec)

                    # Must use @with_pool for pooled model
                    y_pooled = @with_pool pool pooled(x_vec)

                    @test y_orig isa AbstractVector
                    @test y_pooled isa AbstractVector
                    @test size(y_orig) == size(y_pooled)
                    @test y_orig == y_pooled  # bit-exact same result
                end

                @testset "Matrix input (batch)" begin
                    y_orig = model.fluxmodel(x_mat)

                    y_pooled = @with_pool pool pooled(x_mat)

                    @test size(y_orig) == size(y_pooled)
                    @test y_orig == y_pooled  # bit-exact same result
                end
            end
        end

        @testset "poolify with varying batch sizes (no max_batch limit)" begin
            ensemble = loadmodel("sat2_em_d3d_azf-1")
            model = ensemble.models[1]
            pooled = poolify(model.fluxmodel)

            # Unlike bufferize, poolify has no max_batch limit
            @testset "Batch size: $batch_size" for batch_size in [1, 5, 10, 50, 100, 500]
                x = generate_valid_input_matrix(model, batch_size)

                y_orig = model.fluxmodel(x)
                y_pooled = @with_pool pool pooled(x)

                @test size(y_orig) == size(y_pooled)
                @test y_orig == y_pooled
            end
        end

        @testset "poolify multiple calls in single @with_pool block" begin
            ensemble = loadmodel("sat2_em_d3d_azf-1")
            model = ensemble.models[1]
            pooled = poolify(model.fluxmodel)

            x1 = generate_valid_input_matrix(model, 10)
            x2 = generate_valid_input_matrix(model, 20)
            x3 = generate_valid_input_matrix(model, 5)

            # Multiple calls should reuse pool memory
            results = @with_pool pool begin
                r1 = copy(pooled(x1))  # copy to keep result after pool rewind
                r2 = copy(pooled(x2))
                r3 = copy(pooled(x3))
                (r1, r2, r3)
            end

            @test size(results[1]) == (length(model.ynames), 10)
            @test size(results[2]) == (length(model.ynames), 20)
            @test size(results[3]) == (length(model.ynames), 5)
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
            pooled = poolify(model.fluxmodel)

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

            n_parallel = count_parallel_add(pooled)
            @test n_parallel == 5  # 5 ResNet blocks

            # Verify PooledParallelAdd works correctly
            x = generate_valid_input_matrix(model, 1)
            y_orig = model.fluxmodel(x)
            y_pooled = @with_pool pool pooled(x)
            @test y_orig == y_pooled
        end
end
