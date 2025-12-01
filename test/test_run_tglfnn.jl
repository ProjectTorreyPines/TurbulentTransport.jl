using GACODE: FluxSolution

@testset "run_tglfnn" begin
    @testset "basic run with InputTGLF" begin
        input_tglf = load_sample_input()
        TurbulentTransport.apply_presets!(input_tglf)

        result = TurbulentTransport.run_tglfnn(
            input_tglf;
            model_filename="sat3_em_d3d_azf-1",
            warn_nn_train_bounds=false
        )

        @test result isa FluxSolution
        @test result.ENERGY_FLUX_e ≈ EXPECTED_RUN_TGLFNN_SAT3.ENERGY_FLUX_e rtol=REGRESSION_RTOL
        @test result.ENERGY_FLUX_i ≈ EXPECTED_RUN_TGLFNN_SAT3.ENERGY_FLUX_i rtol=REGRESSION_RTOL
        @test result.PARTICLE_FLUX_e ≈ EXPECTED_RUN_TGLFNN_SAT3.PARTICLE_FLUX_e rtol=REGRESSION_RTOL
        @test result.STRESS_TOR_i ≈ EXPECTED_RUN_TGLFNN_SAT3.STRESS_TOR_i rtol=REGRESSION_RTOL
    end

    @testset "run with fidelity=:GKNN" begin
        input_tglf = load_sample_input()
        TurbulentTransport.apply_presets!(input_tglf)

        result_tglfnn = TurbulentTransport.run_tglfnn(
            input_tglf;
            model_filename="sat3_em_d3d_azf-1",
            warn_nn_train_bounds=false,
            fidelity=:TGLFNN
        )

        result_gknn = TurbulentTransport.run_tglfnn(
            input_tglf;
            model_filename="sat3_em_d3d_azf-1",
            warn_nn_train_bounds=false,
            fidelity=:GKNN
        )

        # Both should be valid FluxSolutions
        @test result_tglfnn isa FluxSolution
        @test result_gknn isa FluxSolution

        # TGLFNN should match expected values
        @test result_tglfnn.ENERGY_FLUX_e ≈ EXPECTED_RUN_TGLFNN_SAT3.ENERGY_FLUX_e rtol=REGRESSION_RTOL
        @test result_tglfnn.ENERGY_FLUX_i ≈ EXPECTED_RUN_TGLFNN_SAT3.ENERGY_FLUX_i rtol=REGRESSION_RTOL
        @test result_tglfnn.PARTICLE_FLUX_e ≈ EXPECTED_RUN_TGLFNN_SAT3.PARTICLE_FLUX_e rtol=REGRESSION_RTOL
        @test result_tglfnn.STRESS_TOR_i ≈ EXPECTED_RUN_TGLFNN_SAT3.STRESS_TOR_i rtol=REGRESSION_RTOL

        # GKNN applies correction factors, should match GKNN expected values
        @test result_gknn.ENERGY_FLUX_e ≈ EXPECTED_RUN_TGLFNN_SAT3_GKNN.ENERGY_FLUX_e rtol=REGRESSION_RTOL
        @test result_gknn.ENERGY_FLUX_i ≈ EXPECTED_RUN_TGLFNN_SAT3_GKNN.ENERGY_FLUX_i rtol=REGRESSION_RTOL
        @test result_gknn.PARTICLE_FLUX_e ≈ EXPECTED_RUN_TGLFNN_SAT3_GKNN.PARTICLE_FLUX_e rtol=REGRESSION_RTOL
        @test result_gknn.STRESS_TOR_i ≈ EXPECTED_RUN_TGLFNN_SAT3_GKNN.STRESS_TOR_i rtol=REGRESSION_RTOL
    end

    @testset "run with uncertain=true" begin
        input_tglf = load_sample_input()
        TurbulentTransport.apply_presets!(input_tglf)

        result = TurbulentTransport.run_tglfnn(
            input_tglf;
            model_filename="sat3_em_d3d_azf-1",
            warn_nn_train_bounds=false,
            uncertain=true
        )

        @test result isa FluxSolution
    end

    @testset "batch run with Vector{InputTGLF}" begin
        input_tglf1 = load_sample_input()
        TurbulentTransport.apply_presets!(input_tglf1)

        input_tglf2 = load_sample_input()
        TurbulentTransport.apply_presets!(input_tglf2)
        input_tglf2.BETAE = input_tglf2.BETAE * 1.1

        input_tglfs = [input_tglf1, input_tglf2]

        results = TurbulentTransport.run_tglfnn(
            input_tglfs;
            model_filename="sat3_em_d3d_azf-1",
            warn_nn_train_bounds=false
        )

        @test results isa Vector
        @test length(results) == 2
        @test all(r -> r isa FluxSolution, results)
    end

    @testset "run with Dict input" begin
        input_tglf = load_sample_input()
        TurbulentTransport.apply_presets!(input_tglf)

        # Get model to know required xnames
        model = loadmodel("sat3_em_d3d_azf-1")
        xnames = [replace(name, "_log10" => "") for name in model.xnames]

        # Build dict with array values
        data = Dict{String,Vector{Float64}}()
        for name in xnames
            field_sym = Symbol(name)
            val = getfield(input_tglf, field_sym)
            data[name] = [val, val * 1.01]  # Two samples
        end

        result = TurbulentTransport.run_tglfnn(
            data;
            model_filename="sat3_em_d3d_azf-1",
            warn_nn_train_bounds=false
        )

        @test result isa Dict
        # Output keys use short names like "G_elec", "Q_elec", "P_ions", etc.
        @test haskey(result, "G_elec") || haskey(result, "PARTICLE_FLUX_e")
    end

    @testset "deterministic results" begin
        input_tglf = load_sample_input()
        TurbulentTransport.apply_presets!(input_tglf)

        result1 = TurbulentTransport.run_tglfnn(
            input_tglf;
            model_filename="sat3_em_d3d_azf-1",
            warn_nn_train_bounds=false
        )

        result2 = TurbulentTransport.run_tglfnn(
            input_tglf;
            model_filename="sat3_em_d3d_azf-1",
            warn_nn_train_bounds=false
        )

        @test result1.ENERGY_FLUX_e ≈ result2.ENERGY_FLUX_e
        @test result1.ENERGY_FLUX_i ≈ result2.ENERGY_FLUX_i
        @test result1.PARTICLE_FLUX_e ≈ result2.PARTICLE_FLUX_e
        @test result1.STRESS_TOR_i ≈ result2.STRESS_TOR_i
    end

    @testset "different models" begin
        input_tglf = load_sample_input()
        TurbulentTransport.apply_presets!(input_tglf)

        # Test sat3_em_d3d_azf-1
        result_sat3 = TurbulentTransport.run_tglfnn(
            input_tglf;
            model_filename="sat3_em_d3d_azf-1",
            warn_nn_train_bounds=false
        )
        @test result_sat3 isa FluxSolution
        @test result_sat3.ENERGY_FLUX_e ≈ EXPECTED_RUN_TGLFNN_SAT3.ENERGY_FLUX_e rtol=REGRESSION_RTOL
        @test result_sat3.ENERGY_FLUX_i ≈ EXPECTED_RUN_TGLFNN_SAT3.ENERGY_FLUX_i rtol=REGRESSION_RTOL
        @test result_sat3.PARTICLE_FLUX_e ≈ EXPECTED_RUN_TGLFNN_SAT3.PARTICLE_FLUX_e rtol=REGRESSION_RTOL
        @test result_sat3.STRESS_TOR_i ≈ EXPECTED_RUN_TGLFNN_SAT3.STRESS_TOR_i rtol=REGRESSION_RTOL

        # Test sat2_em_d3d_azf-1
        result_sat2 = TurbulentTransport.run_tglfnn(
            input_tglf;
            model_filename="sat2_em_d3d_azf-1",
            warn_nn_train_bounds=false
        )
        @test result_sat2 isa FluxSolution
        @test result_sat2.ENERGY_FLUX_e ≈ EXPECTED_RUN_TGLFNN_SAT2.ENERGY_FLUX_e rtol=REGRESSION_RTOL
        @test result_sat2.ENERGY_FLUX_i ≈ EXPECTED_RUN_TGLFNN_SAT2.ENERGY_FLUX_i rtol=REGRESSION_RTOL
        @test result_sat2.PARTICLE_FLUX_e ≈ EXPECTED_RUN_TGLFNN_SAT2.PARTICLE_FLUX_e rtol=REGRESSION_RTOL
        @test result_sat2.STRESS_TOR_i ≈ EXPECTED_RUN_TGLFNN_SAT2.STRESS_TOR_i rtol=REGRESSION_RTOL
    end
end
