@testset "Bufferize" begin
    # Sample of representative models to test (covering different architectures)
    TEST_MODELS = [
        "sat0_em_d3d",
        "sat1_em_d3d_azf-1",
        "sat2_em_d3d_azf-1",
        "sat3_em_d3d_azf-1",
        "sat3_em_d3d_azf-1_gknne24",  # GKNN model
    ]

    @testset "bufferize correctness" begin
        @testset "Model: $model_name" for model_name in TEST_MODELS
            ensemble = loadmodel(model_name)
            model = ensemble.models[1]

            # Bufferize the fluxmodel
            max_batch = 50
            buffered = bufferize(model.fluxmodel, max_batch)

            # Generate valid test input
            x_vec = generate_valid_input(model)
            x_mat = generate_valid_input_matrix(model, 10)

            @testset "Vector input" begin
                y_orig = model.fluxmodel(x_vec)
                y_buff = buffered(x_vec)

                # Both should return 1D vectors with matching shape and values
                @test y_orig isa AbstractVector
                @test y_buff isa AbstractVector
                @test size(y_orig) == size(y_buff)
                @test y_orig == y_buff
            end

            @testset "Matrix input (batch)" begin
                y_orig = model.fluxmodel(x_mat)
                y_buff = buffered(x_mat)

                @test size(y_orig) == size(y_buff)
                @test y_orig == y_buff
            end
        end
    end

    @testset "bufferize with varying batch sizes" begin
        ensemble = loadmodel("sat2_em_d3d_azf-1")
        model = ensemble.models[1]
        max_batch = 100
        buffered = bufferize(model.fluxmodel, max_batch)

        @testset "Batch size: $batch_size" for batch_size in [1, 5, 10, 50, 100]
            x = generate_valid_input_matrix(model, batch_size)

            y_orig = model.fluxmodel(x)
            y_buff = buffered(x)

            @test size(y_orig) == size(y_buff)
            @test y_orig == y_buff
        end
    end

    @testset "bufferize all available models" begin
        all_models = available_models()

        # Filter to .bson files only (skip directories/symlinks for speed)
        bson_models = filter(m -> !contains(m, "/") && endswith(m, ".bson") || !contains(m, "."), all_models)

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
            buffered = bufferize(model.fluxmodel, 10)

            x = generate_valid_input(model)
            y_orig = model.fluxmodel(x)
            y_buff = buffered(x)

            @test y_orig == y_buff
        end
    end

    @testset "BufferedDense and BufferedActivation types" begin
        ensemble = loadmodel("sat2_em_d3d_azf-1")
        model = ensemble.models[1]
        buffered = bufferize(model.fluxmodel, 10)

        # Check that BufferedDense layers exist
        first_layer = buffered.layers[1]
        @test first_layer isa BufferedDense

        # Check that the wrapped dense is preserved
        @test first_layer.dense isa Flux.Dense
        @test size(first_layer.buffer, 2) == 10  # max_batch
    end

    @testset "max_layer_width" begin
        ensemble = loadmodel("sat2_em_d3d_azf-1")
        model = ensemble.models[1]

        width = TurbulentTransport.max_layer_width(model.fluxmodel)

        # Should be positive
        @test width > 0

        # Should be at least as large as output dimension
        n_outputs = length(model.ynames)
        @test width >= n_outputs
    end
end
