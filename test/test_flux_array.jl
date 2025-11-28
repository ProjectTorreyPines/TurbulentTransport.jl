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
        @test result == EXPECTED_FLUX_ARRAY_SINGLE
    end

    @testset "single model with matrix input" begin
        n_samples = 5
        x = generate_valid_input_matrix(single_model, n_samples)

        result = flux_array(single_model, x; warn_nn_train_bounds=false)

        @test result isa AbstractMatrix
        @test size(result, 1) == length(single_model.ynames)
        @test size(result, 2) == n_samples
        @test all(isfinite, result)
        # All columns should be identical since input is identical
        @test result[:, 1] == EXPECTED_FLUX_ARRAY_SINGLE
        for j in 2:n_samples
            @test result[:, j] == result[:, 1]
        end
    end

    @testset "ensemble with vector input" begin
        x = generate_valid_input(ensemble_model)

        # Test without uncertainty
        result = flux_array(ensemble_model, x; uncertain=false, warn_nn_train_bounds=false)
        @test result isa AbstractVector
        @test all(isfinite, result)
        @test result == EXPECTED_FLUX_ARRAY_ENSEMBLE

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
        # All columns should be identical since input is identical
        @test result[:, 1] == EXPECTED_FLUX_ARRAY_ENSEMBLE
        for j in 2:n_samples
            @test result[:, j] == result[:, 1]
        end
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
        @test vec(result) == EXPECTED_FLUX_ARRAY_ENSEMBLE
    end

    @testset "fidelity modes" begin
        x = generate_valid_input(ensemble_model)

        result_tglfnn = flux_array(ensemble_model, x; fidelity=:TGLFNN, warn_nn_train_bounds=false)
        @test length(result_tglfnn) == length(ensemble_model.ynames)
        @test result_tglfnn == EXPECTED_FLUX_ARRAY_ENSEMBLE

        # Note: fidelity=:GKNN on ensemble has dimension mismatch bug when
        # nouts differs between modes. Test TGLFNN mode only for now.
    end

    @testset "error on NaN input" begin
        x = generate_valid_input(single_model)
        x[1] = NaN

        @test_throws ErrorException flux_array(single_model, x; warn_nn_train_bounds=true)
    end

    @testset "error on Inf input" begin
        x = generate_valid_input(single_model)
        x[1] = Inf

        @test_throws ErrorException flux_array(single_model, x; warn_nn_train_bounds=true)
    end

    @testset "error on negative input for log10 field" begin
        x = generate_valid_input(single_model)
        # Find a log10 field and set it to negative (which will produce NaN after log10)
        for (i, name) in enumerate(single_model.xnames)
            if contains(name, "_log10")
                x[i] = -1.0  # log10(-1) throws DomainError
                break
            end
        end

        # Julia's log10 throws DomainError for negative inputs
        @test_throws DomainError flux_array(single_model, x; warn_nn_train_bounds=true)
    end

    @testset "GKNN fidelity on single model" begin
        x = generate_valid_input(single_model)

        # GKNN mode returns raw NN output (no denormalization)
        result_gknn = flux_array(single_model, x; fidelity=:GKNN, warn_nn_train_bounds=false)

        @test result_gknn isa AbstractVector
        @test length(result_gknn) == length(single_model.ynames)
        @test all(isfinite, result_gknn)
        @test result_gknn == EXPECTED_FLUX_ARRAY_SINGLE_GKNN

        # GKNN and TGLFNN should produce different results
        result_tglfnn = flux_array(single_model, x; fidelity=:TGLFNN, warn_nn_train_bounds=false)
        @test result_gknn != result_tglfnn
    end

    @testset "GKNN fidelity on single model with matrix input" begin
        n_samples = 3
        x = generate_valid_input_matrix(single_model, n_samples)

        result_gknn = flux_array(single_model, x; fidelity=:GKNN, warn_nn_train_bounds=false)

        @test result_gknn isa AbstractMatrix
        @test size(result_gknn, 1) == length(single_model.ynames)
        @test size(result_gknn, 2) == n_samples
        @test all(isfinite, result_gknn)
        @test result_gknn[:, 1] == EXPECTED_FLUX_ARRAY_SINGLE_GKNN
    end

    @testset "log10 transformation applied correctly" begin
        # Verify that log10 fields are transformed before NN input
        x = generate_valid_input(single_model)

        # Find log10 fields and verify their indices
        log10_indices = Int[]
        for (i, name) in enumerate(single_model.xnames)
            if contains(name, "_log10")
                push!(log10_indices, i)
            end
        end

        @test !isempty(log10_indices)  # Model should have log10 fields

        # Verify values at bounds work correctly
        for idx in log10_indices
            # Lower bound in transformed space
            lower_transformed = single_model.xbounds[idx, 1]
            # Original value that maps to lower bound
            x_lower = copy(x)
            x_lower[idx] = 10.0^lower_transformed

            result = flux_array(single_model, x_lower; warn_nn_train_bounds=false)
            @test all(isfinite, result)
        end
    end

    @testset "normalization verification" begin
        # The difference between GKNN and TGLFNN modes shows normalization
        # TGLFNN: yy = yy * yσ + ym (denormalized)
        # GKNN: yy = raw output (normalized)
        x = generate_valid_input(single_model)

        result_raw = flux_array(single_model, x; fidelity=:GKNN, warn_nn_train_bounds=false)
        result_denorm = flux_array(single_model, x; fidelity=:TGLFNN, warn_nn_train_bounds=false)

        # Manually apply denormalization: yy * yσ + ym
        manual_denorm = result_raw .* single_model.yσ .+ single_model.ym
        @test manual_denorm ≈ result_denorm
    end

    @testset "matrix input with varying values" begin
        n_samples = 4
        base_input = generate_valid_input(single_model)

        # Create matrix with slightly different values per column
        x = zeros(Float64, length(base_input), n_samples)
        for j in 1:n_samples
            # Scale inputs by different factors
            scale = 0.9 + 0.05 * j  # 0.95, 1.0, 1.05, 1.1
            x[:, j] = base_input .* scale
        end

        result = flux_array(single_model, x; warn_nn_train_bounds=false)

        @test size(result, 2) == n_samples
        # Different inputs should give different outputs
        @test result[:, 1] != result[:, n_samples]
        @test all(isfinite, result)
    end

    @testset "Float32 input type" begin
        x = Float32.(generate_valid_input(single_model))

        result = flux_array(single_model, x; warn_nn_train_bounds=false)

        @test all(isfinite, result)
        # Output type depends on model (Float64), not input
        @test result isa AbstractVector
    end

    @testset "ensemble uncertainty returns Measurement type" begin
        x = generate_valid_input(ensemble_model)

        result_uncertain = flux_array(ensemble_model, x; uncertain=true, warn_nn_train_bounds=false)

        # Check that results have uncertainty
        @test all(r -> r isa Measurements.Measurement, result_uncertain)

        # Extract values and check they're finite
        values = Measurements.value.(result_uncertain)
        @test all(isfinite, values)

        # Uncertainties should be non-negative
        uncertainties = Measurements.uncertainty.(result_uncertain)
        @test all(u -> u >= 0, uncertainties)
    end

    @testset "warn_nn_train_bounds parameter" begin
        x = generate_valid_input(single_model)

        # With bounds checking disabled, no warnings expected
        result_no_warn = flux_array(single_model, x; warn_nn_train_bounds=false)
        @test all(isfinite, result_no_warn)

        # With bounds checking enabled, should still work for valid input
        result_with_warn = flux_array(single_model, x; warn_nn_train_bounds=true)
        @test result_no_warn == result_with_warn
    end
end
