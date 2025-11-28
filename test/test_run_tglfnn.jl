using GACODE: FluxSolution

@testset "run_tglfnn" begin
    @testset "basic run with InputTGLF" begin
        filepath, tmpdir = create_temp_input_tglf()
        try
            input_tglf = TurbulentTransport.load(InputTGLF(), filepath)
            TurbulentTransport.apply_presets!(input_tglf)

            result = TurbulentTransport.run_tglfnn(
                input_tglf;
                model_filename="sat3_em_d3d_azf-1",
                warn_nn_train_bounds=false
            )

            @test result isa FluxSolution
            @test isfinite(result.ENERGY_FLUX_e)
            @test isfinite(result.ENERGY_FLUX_i)
            @test isfinite(result.PARTICLE_FLUX_e)
            @test isfinite(result.STRESS_TOR_i)
        finally
            cleanup_temp_dir(tmpdir)
        end
    end

    @testset "run with fidelity=:GKNN" begin
        filepath, tmpdir = create_temp_input_tglf()
        try
            input_tglf = TurbulentTransport.load(InputTGLF(), filepath)
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

            # GKNN applies correction factors, so results should differ
            @test result_tglfnn.ENERGY_FLUX_e != result_gknn.ENERGY_FLUX_e
        finally
            cleanup_temp_dir(tmpdir)
        end
    end

    @testset "run with uncertain=true" begin
        filepath, tmpdir = create_temp_input_tglf()
        try
            input_tglf = TurbulentTransport.load(InputTGLF(), filepath)
            TurbulentTransport.apply_presets!(input_tglf)

            result = TurbulentTransport.run_tglfnn(
                input_tglf;
                model_filename="sat3_em_d3d_azf-1",
                warn_nn_train_bounds=false,
                uncertain=true
            )

            @test result isa FluxSolution
        finally
            cleanup_temp_dir(tmpdir)
        end
    end

    @testset "batch run with Vector{InputTGLF}" begin
        filepath, tmpdir = create_temp_input_tglf()
        try
            input_tglf1 = TurbulentTransport.load(InputTGLF(), filepath)
            TurbulentTransport.apply_presets!(input_tglf1)

            input_tglf2 = TurbulentTransport.load(InputTGLF(), filepath)
            TurbulentTransport.apply_presets!(input_tglf2)
            # Modify slightly
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
        finally
            cleanup_temp_dir(tmpdir)
        end
    end

    @testset "run with Dict input" begin
        filepath, tmpdir = create_temp_input_tglf()
        try
            input_tglf = TurbulentTransport.load(InputTGLF(), filepath)
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
        finally
            cleanup_temp_dir(tmpdir)
        end
    end

    @testset "deterministic results" begin
        filepath, tmpdir = create_temp_input_tglf()
        try
            input_tglf = TurbulentTransport.load(InputTGLF(), filepath)
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
        finally
            cleanup_temp_dir(tmpdir)
        end
    end

    @testset "different models" begin
        filepath, tmpdir = create_temp_input_tglf()
        try
            input_tglf = TurbulentTransport.load(InputTGLF(), filepath)
            TurbulentTransport.apply_presets!(input_tglf)

            # Test a few different model types
            for model_name in ["sat3_em_d3d_azf-1", "sat2_em_d3d_azf-1"]
                result = TurbulentTransport.run_tglfnn(
                    input_tglf;
                    model_filename=model_name,
                    warn_nn_train_bounds=false
                )
                @test result isa FluxSolution
                @test all(isfinite, [result.ENERGY_FLUX_e, result.ENERGY_FLUX_i,
                                     result.PARTICLE_FLUX_e, result.STRESS_TOR_i])
            end
        finally
            cleanup_temp_dir(tmpdir)
        end
    end
end
