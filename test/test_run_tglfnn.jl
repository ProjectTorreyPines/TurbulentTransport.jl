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

"""
Regression tests for run_tglfnn across models and fidelity modes.

Tests single vs vector version equivalence and regression against baseline values.
Expected values defined in fixtures.jl (REGRESSION_EXPECTED_VALUES).
"""
# Helper: test flux values against expected
function test_regression_flux_values(sol::FluxSolution, expected::NamedTuple; rtol=REGRESSION_RTOL)
    @test isapprox(sol.ENERGY_FLUX_e, expected.ENERGY_FLUX_e; rtol=rtol)
    @test isapprox(sol.ENERGY_FLUX_i, expected.ENERGY_FLUX_i; rtol=rtol)
    @test isapprox(sol.PARTICLE_FLUX_e, expected.PARTICLE_FLUX_e; rtol=rtol)
    @test isapprox(sol.STRESS_TOR_i, expected.STRESS_TOR_i; rtol=rtol)
end

@testset "run_tglfnn Regression" begin
    inputs = create_regression_inputs()

    for (model_filename, fidelity, description) in REGRESSION_MODEL_CONFIGS
        @testset "$description ($model_filename)" begin

            @testset "Single version values" begin
                for (i, input) in enumerate(inputs)
                    sol = TurbulentTransport.run_tglfnn(
                        input;
                        model_filename=model_filename,
                        warn_nn_train_bounds=false,
                        fidelity=fidelity
                    )
                    expected = REGRESSION_EXPECTED_VALUES[(model_filename, fidelity, i)]
                    test_regression_flux_values(sol, expected)
                end
            end

            @testset "Vector version values" begin
                sols = TurbulentTransport.run_tglfnn(
                    inputs;
                    model_filename=model_filename,
                    warn_nn_train_bounds=false,
                    fidelity=fidelity
                )
                for (i, sol) in enumerate(sols)
                    expected = REGRESSION_EXPECTED_VALUES[(model_filename, fidelity, i)]
                    test_regression_flux_values(sol, expected)
                end
            end

            @testset "Single vs Vector equivalence" begin
                single_results = [
                    TurbulentTransport.run_tglfnn(
                        input;
                        model_filename=model_filename,
                        warn_nn_train_bounds=false,
                        fidelity=fidelity
                    ) for input in inputs
                ]
                vector_results = TurbulentTransport.run_tglfnn(
                    inputs;
                    model_filename=model_filename,
                    warn_nn_train_bounds=false,
                    fidelity=fidelity
                )

                for (s, v) in zip(single_results, vector_results)
                    @test isapprox(s.ENERGY_FLUX_e, v.ENERGY_FLUX_e; rtol=REGRESSION_RTOL)
                    @test isapprox(s.ENERGY_FLUX_i, v.ENERGY_FLUX_i; rtol=REGRESSION_RTOL)
                    @test isapprox(s.PARTICLE_FLUX_e, v.PARTICLE_FLUX_e; rtol=REGRESSION_RTOL)
                    @test isapprox(s.STRESS_TOR_i, v.STRESS_TOR_i; rtol=REGRESSION_RTOL)
                end
            end
        end
    end
end

@testset "run_tglfnn edge cases" begin
    @testset "Empty input vector" begin
        # run_tglfnn with empty input should return empty results (not error)
        empty_inputs = InputTGLF{Float64}[]

        result = TurbulentTransport.run_tglfnn(
            empty_inputs;
            model_filename="sat3_em_d3d_azf-1",
            warn_nn_train_bounds=false,
            fidelity=:TGLFNN
        )

        @test result isa Vector
        @test isempty(result)
    end
end

@testset "run_tglfnn radial-dependent models" begin
    # Test models with radial-dependent blending (d3d, d3dnearedge, d3dedge variants)
    radial_models = (
        "sat0quench_em_d3d_azf+1_withnegD",
        "sat1_em_d3d_azf-1_withnegD",
        "sat2_em_d3d_azf-1_withnegD",
        "sat3_em_d3d_azf-1_withnegD"
    )

    @testset "Model: $model_name" for model_name in radial_models
        # Create test inputs with different RMIN_LOC values covering all radial regions
        input_core = load_sample_input()  # rmin ~ 0.573 (< 0.68, core region)
        TurbulentTransport.apply_presets!(input_core)
        
        input_overlap1 = deepcopy(input_core)
        input_overlap1.RMIN_LOC = 0.75  # 0.68 <= rmin < 0.881 (overlap: all three models)
        
        input_overlap2 = deepcopy(input_core)
        input_overlap2.RMIN_LOC = 0.92  # 0.881 <= rmin < 0.975 (overlap: nearedge & edge)
        
        input_edge = deepcopy(input_core)
        input_edge.RMIN_LOC = 0.98  # rmin >= 0.975 (edge only)

        @testset "Single input - core region (rmin < 0.68)" begin
            result = TurbulentTransport.run_tglfnn(
                input_core;
                model_filename=model_name,
                warn_nn_train_bounds=false
            )
            @test result isa FluxSolution
            @test isfinite(result.ENERGY_FLUX_e)
            @test isfinite(result.ENERGY_FLUX_i)
            @test isfinite(result.PARTICLE_FLUX_e)
            @test isfinite(result.STRESS_TOR_i)
        end

        @testset "Single input - overlap region 1 (0.68 <= rmin < 0.881)" begin
            result = TurbulentTransport.run_tglfnn(
                input_overlap1;
                model_filename=model_name,
                warn_nn_train_bounds=false
            )
            @test result isa FluxSolution
            @test isfinite(result.ENERGY_FLUX_e)
            @test isfinite(result.ENERGY_FLUX_i)
            @test isfinite(result.PARTICLE_FLUX_e)
            @test isfinite(result.STRESS_TOR_i)
        end

        @testset "Single input - overlap region 2 (0.881 <= rmin < 0.975)" begin
            result = TurbulentTransport.run_tglfnn(
                input_overlap2;
                model_filename=model_name,
                warn_nn_train_bounds=false
            )
            @test result isa FluxSolution
            @test isfinite(result.ENERGY_FLUX_e)
            @test isfinite(result.ENERGY_FLUX_i)
            @test isfinite(result.PARTICLE_FLUX_e)
            @test isfinite(result.STRESS_TOR_i)
        end

        @testset "Single input - edge region (rmin >= 0.975)" begin
            result = TurbulentTransport.run_tglfnn(
                input_edge;
                model_filename=model_name,
                warn_nn_train_bounds=false
            )
            @test result isa FluxSolution
            @test isfinite(result.ENERGY_FLUX_e)
            @test isfinite(result.ENERGY_FLUX_i)
            @test isfinite(result.PARTICLE_FLUX_e)
            @test isfinite(result.STRESS_TOR_i)
        end

        @testset "Batch input - multiple radial regions" begin
            inputs = [input_core, input_overlap1, input_overlap2, input_edge]
            results = TurbulentTransport.run_tglfnn(
                inputs;
                model_filename=model_name,
                warn_nn_train_bounds=false
            )
            
            @test results isa Vector
            @test length(results) == 4
            @test all(r -> r isa FluxSolution, results)
            
            # Check that all flux values are finite
            for result in results
                @test isfinite(result.ENERGY_FLUX_e)
                @test isfinite(result.ENERGY_FLUX_i)
                @test isfinite(result.PARTICLE_FLUX_e)
                @test isfinite(result.STRESS_TOR_i)
            end
        end

        @testset "Batch vs single equivalence" begin
            inputs = [input_core, input_overlap1, input_overlap2, input_edge]
            
            # Run individually
            single_results = [
                TurbulentTransport.run_tglfnn(
                    input;
                    model_filename=model_name,
                    warn_nn_train_bounds=false
                ) for input in inputs
            ]
            
            # Run as batch
            batch_results = TurbulentTransport.run_tglfnn(
                inputs;
                model_filename=model_name,
                warn_nn_train_bounds=false
            )
            
            # Compare results (should be identical)
            for (single, batch) in zip(single_results, batch_results)
                @test single.ENERGY_FLUX_e ≈ batch.ENERGY_FLUX_e rtol=1e-14
                @test single.ENERGY_FLUX_i ≈ batch.ENERGY_FLUX_i rtol=1e-14
                @test single.PARTICLE_FLUX_e ≈ batch.PARTICLE_FLUX_e rtol=1e-14
                @test single.STRESS_TOR_i ≈ batch.STRESS_TOR_i rtol=1e-14
            end
        end
    end
end
