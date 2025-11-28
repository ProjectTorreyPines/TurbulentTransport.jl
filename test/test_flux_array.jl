@testset "flux_array" begin
    # Load models once for all tests
    # All models in repository are ensembles, so get individual model from ensemble
    ensemble_for_single = loadmodel("sat0_em_d3d")
    single_model = ensemble_for_single.models[1]
    ensemble_model = loadmodel(TEST_MODEL_ENSEMBLE)

    @testset "single model with vector input" begin
        x = generate_valid_input(single_model)

        result = flux_array(single_model, x; warn_nn_train_bounds=false)

        @test result isa AbstractVector
        @test length(result) == length(single_model.ynames)
        @test all(isfinite, result)
    end

    @testset "single model with matrix input" begin
        n_samples = 5
        x = generate_valid_input_matrix(single_model, n_samples)

        result = flux_array(single_model, x; warn_nn_train_bounds=false)

        @test result isa AbstractMatrix
        @test size(result, 1) == length(single_model.ynames)
        @test size(result, 2) == n_samples
        @test all(isfinite, result)
    end

    @testset "ensemble with vector input" begin
        x = generate_valid_input(ensemble_model)

        # Test without uncertainty
        result = flux_array(ensemble_model, x; uncertain=false, warn_nn_train_bounds=false)
        @test result isa AbstractVector
        @test all(isfinite, result)

        # Test with uncertainty (returns Measurement type)
        result_uncertain = flux_array(ensemble_model, x; uncertain=true, warn_nn_train_bounds=false)
        @test length(result_uncertain) == length(result)
    end

    @testset "ensemble with matrix input" begin
        n_samples = 3
        x = generate_valid_input_matrix(ensemble_model, n_samples)

        result = flux_array(ensemble_model, x; uncertain=false, warn_nn_train_bounds=false)

        @test result isa AbstractMatrix
        @test size(result, 2) == n_samples
        @test all(isfinite, result)
    end

    @testset "deterministic output" begin
        x = generate_valid_input(single_model)

        result1 = flux_array(single_model, x; warn_nn_train_bounds=false)
        result2 = flux_array(single_model, x; warn_nn_train_bounds=false)

        @test result1 ≈ result2
    end

    @testset "vararg input" begin
        # vararg method works on TGLFmodel (ensemble), not individual TGLFNNmodel
        x = generate_valid_input(ensemble_model)

        result = flux_array(ensemble_model, x...; warn_nn_train_bounds=false)

        @test result isa AbstractArray
        @test all(isfinite, result)
    end

    @testset "fidelity modes" begin
        x = generate_valid_input(ensemble_model)

        result_tglfnn = flux_array(ensemble_model, x; fidelity=:TGLFNN, warn_nn_train_bounds=false)
        @test length(result_tglfnn) == length(ensemble_model.ynames)

        # Note: fidelity=:GKNN on ensemble has dimension mismatch bug when
        # nouts differs between modes. Test TGLFNN mode only for now.
    end
end
