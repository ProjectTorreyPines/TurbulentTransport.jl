@testset "apply_presets!" begin
    @testset "SAT_RULE = 3" begin
        input_tglf = InputTGLF{Float64}()
        input_tglf.SAT_RULE = 3
        input_tglf.UNITS = "GYRO"
        input_tglf.USE_BPER = false

        TurbulentTransport.apply_presets!(input_tglf)

        @test input_tglf.XNU_MODEL == 3
        @test input_tglf.WDIA_TRAPPED == 1.0
        @test input_tglf.UNITS == "CGYRO"
    end

    @testset "SAT_RULE = 2" begin
        input_tglf = InputTGLF{Float64}()
        input_tglf.SAT_RULE = 2
        input_tglf.UNITS = "GYRO"
        input_tglf.USE_BPER = false

        TurbulentTransport.apply_presets!(input_tglf)

        @test input_tglf.XNU_MODEL == 3
        @test input_tglf.WDIA_TRAPPED == 1.0
        @test input_tglf.UNITS == "CGYRO"
    end

    @testset "SAT_RULE = 1" begin
        input_tglf = InputTGLF{Float64}()
        input_tglf.SAT_RULE = 1
        input_tglf.USE_BPER = false

        TurbulentTransport.apply_presets!(input_tglf)

        @test input_tglf.XNU_MODEL == 2
    end

    @testset "SAT_RULE = 0" begin
        input_tglf = InputTGLF{Float64}()
        input_tglf.SAT_RULE = 0
        input_tglf.NMODES = 6
        input_tglf.USE_BPER = false

        TurbulentTransport.apply_presets!(input_tglf)

        @test input_tglf.UNITS == "GYRO"
        @test input_tglf.XNU_MODEL == 2
        @test input_tglf.NMODES == 4  # Should be adjusted when > 2
    end

    @testset "USE_BPER sets ALPHA_MACH" begin
        input_tglf = InputTGLF{Float64}()
        input_tglf.SAT_RULE = 1
        input_tglf.USE_BPER = true
        input_tglf.ALPHA_MACH = 1.0

        TurbulentTransport.apply_presets!(input_tglf)

        @test input_tglf.ALPHA_MACH == 0.0
    end

    @testset "full workflow: load then apply_presets!" begin
        filepath, tmpdir = create_temp_input_tglf()
        try
            input_tglf = TurbulentTransport.load(InputTGLF(), filepath)
            TurbulentTransport.apply_presets!(input_tglf)

            # SAT_RULE=3 with USE_BPER=true
            @test input_tglf.XNU_MODEL == 3
            @test input_tglf.WDIA_TRAPPED == 1.0
            @test input_tglf.UNITS == "CGYRO"
            @test input_tglf.ALPHA_MACH == 0.0
        finally
            cleanup_temp_dir(tmpdir)
        end
    end
end
