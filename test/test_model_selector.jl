@testset "model_selector" begin
    @testset "basic functionality with Vector{InputTGLF}" begin
        # Load sample inputs
        inputs = create_regression_inputs()

        # Run model_selector with a known model
        # Note: rho values automatically extracted from RMIN_LOC in each InputTGLF
        results = model_selector(
            inputs;
            filter_sat_rule=:sat3,
            electromagnetic=true,
            max_models=3,
            verbose=false
        )

        # Test return structure
        @test results isa NamedTuple
        @test haskey(results, :rho_grid)
        @test haskey(results, :rankings)
        @test haskey(results, :input_tglfs)
        @test haskey(results, :all_results)

        # Test rho_grid (should contain RMIN_LOC from each InputTGLF)
        @test length(results.rho_grid) == 3
        @test length(results.rankings) == 3

        # Test rankings structure
        for ranking in results.rankings
            @test haskey(ranking, :rho)
            @test haskey(ranking, :top_models)
            @test haskey(ranking, :confidences)
            @test haskey(ranking, :flux_outputs)

            # Should have at most 3 models (max_models=3)
            @test length(ranking.top_models) <= 3
            @test length(ranking.top_models) == length(ranking.confidences)
            @test length(ranking.top_models) == length(ranking.flux_outputs)

            # Confidences should be sorted (lower is better)
            @test issorted(ranking.confidences)

            # Each model name should start with "sat3_em_" (filtered)
            for model in ranking.top_models
                @test startswith(model, "sat3_em_")
                @test !occursin("gknn", model)
                @test !occursin("qlnn", model)
                @test !occursin("edge", model)
            end

            # Flux outputs should be FluxSolution with Measurements
            for flux_sol in ranking.flux_outputs
                @test flux_sol isa GACODE.FluxSolution
                @test flux_sol.PARTICLE_FLUX_e isa Measurements.Measurement
                @test flux_sol.ENERGY_FLUX_e isa Measurements.Measurement
                @test flux_sol.ENERGY_FLUX_i isa Measurements.Measurement
            end
        end

        # Test all_results dict
        @test results.all_results isa Dict{String, Tuple{Bool, Any}}

        # Count successful and failed models
        n_success = count(x -> x[1], values(results.all_results))
        n_failed = count(x -> !x[1], values(results.all_results))

        @test n_success > 0  # Should have at least some successful models
        @test n_success + n_failed == length(results.all_results)
    end

    @testset "single InputTGLF" begin
        # Test single InputTGLF method
        input = load_sample_input()
        TurbulentTransport.apply_presets!(input)

        # Note: rho value automatically extracted from RMIN_LOC
        result = model_selector(
            input;
            filter_sat_rule=:sat3,
            electromagnetic=true,
            max_models=3,
            verbose=false
        )

        # Test simplified return structure for single input
        @test result isa NamedTuple
        @test haskey(result, :rho)
        @test haskey(result, :top_models)
        @test haskey(result, :confidences)
        @test haskey(result, :flux_outputs)
        @test haskey(result, :input_tglf)
        @test haskey(result, :all_results)

        # Test values
        @test result.rho == input.RMIN_LOC  # Should match RMIN_LOC from input
        @test length(result.top_models) <= 3
        @test length(result.top_models) == length(result.confidences)
        @test length(result.top_models) == length(result.flux_outputs)
        @test issorted(result.confidences)

        # Test input_tglf is preserved
        @test result.input_tglf === input
    end

    @testset "electromagnetic filtering" begin
        input = load_sample_input()
        TurbulentTransport.apply_presets!(input)

        # Test EM filtering (default)
        results_em = model_selector(
            [input];
            filter_sat_rule=:sat3,
            electromagnetic=true,
            max_models=3,
            verbose=false
        )

        # Test ES filtering
        results_es = model_selector(
            [input];
            filter_sat_rule=:sat3,
            electromagnetic=false,
            max_models=3,
            verbose=false
        )

        # Both should succeed but with different models
        @test length(results_em.rankings) == 1
        @test length(results_es.rankings) == 1

        # EM models should have "_em_" in their names
        for model in results_em.rankings[1].top_models
            @test occursin("_em_", model)
            @test !occursin("_es_", model)
        end

        # ES models should have "_es_" in their names
        for model in results_es.rankings[1].top_models
            @test occursin("_es_", model)
            @test !occursin("_em_", model)
        end
    end

    @testset "saturation rule filtering" begin
        input = load_sample_input()
        TurbulentTransport.apply_presets!(input)

        # Test sat3 filtering
        results_sat3 = model_selector(
            [input];
            filter_sat_rule=:sat3,
            electromagnetic=true,
            max_models=3,
            verbose=false
        )

        # All models should start with "sat3"
        for model in results_sat3.rankings[1].top_models
            @test startswith(model, "sat3")
        end

        # Test sat1 filtering (if models exist)
        sat1_models = filter(m -> startswith(m, "sat1"), available_models())
        if !isempty(sat1_models)
            results_sat1 = model_selector(
                [input];
                filter_sat_rule=:sat1,
                electromagnetic=true,
                max_models=3,
                verbose=false
            )

            for model in results_sat1.rankings[1].top_models
                @test startswith(model, "sat1")
            end
        end
    end

    @testset "max_models parameter" begin
        input = load_sample_input()
        TurbulentTransport.apply_presets!(input)

        # Test with max_models=1
        results_1 = model_selector(
            [input];
            filter_sat_rule=:sat3,
            electromagnetic=true,
            max_models=1,
            verbose=false
        )

        @test length(results_1.rankings[1].top_models) == 1

        # Test with max_models=5
        results_5 = model_selector(
            [input];
            filter_sat_rule=:sat3,
            electromagnetic=true,
            max_models=5,
            verbose=false
        )

        @test length(results_5.rankings[1].top_models) <= 5
    end

    @testset "confidence metric properties" begin
        # Create inputs with varying uncertainties
        inputs = create_regression_inputs()

        results = model_selector(
            inputs;
            filter_sat_rule=:sat3,
            electromagnetic=true,
            max_models=3,
            verbose=false
        )

        # Test that confidences are positive and reasonable
        for ranking in results.rankings
            for conf in ranking.confidences
                @test conf >= 0.0
                @test isfinite(conf)
            end

            # Confidences should be sorted (lower = more confident)
            @test issorted(ranking.confidences)
        end
    end

    @testset "basic validation" begin
        # Simple validation test - just ensure function runs and returns valid structure
        input = load_sample_input()
        TurbulentTransport.apply_presets!(input)

        results = model_selector(
            [input];
            filter_sat_rule=:sat3,
            electromagnetic=true,
            max_models=3,
            verbose=false
        )

        # Test that we got valid results
        @test !isempty(results.rankings)
        @test !isempty(results.rankings[1].top_models)
        @test all(conf -> conf >= 0 && isfinite(conf), results.rankings[1].confidences)

        # Test that at least some models succeeded
        n_success = count(x -> x[1], values(results.all_results))
        @test n_success >= 3  # Should have at least 3 successful models
    end
end
