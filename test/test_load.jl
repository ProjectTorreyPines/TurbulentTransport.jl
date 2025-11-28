@testset "InputTGLF Loading" begin
    @testset "load from file" begin
        filepath, tmpdir = create_temp_input_tglf()
        try
            input_tglf = TurbulentTransport.load(InputTGLF(), filepath)

            # Test all expected field values
            for (field, expected) in pairs(EXPECTED_LOAD_VALUES)
                actual = getproperty(input_tglf, field)
                @test actual == expected
            end
        finally
            cleanup_temp_dir(tmpdir)
        end
    end

    @testset "load handles boolean formats" begin
        # Test .true./.false. format
        content = """
        NS = 2
        USE_BPER = .true.
        USE_BPAR = .false.
        SAT_RULE = 0
        """
        filepath, tmpdir = create_temp_input_tglf(content)
        try
            input_tglf = TurbulentTransport.load(InputTGLF(), filepath)
            @test input_tglf.USE_BPER == true
            @test input_tglf.USE_BPAR == false
        finally
            cleanup_temp_dir(tmpdir)
        end
    end

    @testset "load handles integer/float conversion" begin
        content = """
        NS = 3
        SAT_RULE = 2
        BETAE = 1.5e-3
        Q_LOC = 2.0
        """
        filepath, tmpdir = create_temp_input_tglf(content)
        try
            input_tglf = TurbulentTransport.load(InputTGLF(), filepath)
            @test input_tglf.NS == 3
            @test input_tglf.SAT_RULE == 2
            @test input_tglf.BETAE == 0.0015
            @test input_tglf.Q_LOC == 2.0
        finally
            cleanup_temp_dir(tmpdir)
        end
    end
end
