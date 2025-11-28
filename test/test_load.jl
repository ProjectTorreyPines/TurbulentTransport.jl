@testset "InputTGLF Loading" begin
    @testset "load from sample file" begin
        input_tglf = load_sample_input()

        # Test all expected field values
        for (field, expected) in pairs(EXPECTED_LOAD_VALUES)
            actual = getproperty(input_tglf, field)
            @test actual == expected
        end
    end

    @testset "load handles boolean formats" begin
        # Create temp file for specific format testing
        content = """
        NS = 2
        USE_BPER = .true.
        USE_BPAR = .false.
        SAT_RULE = 0
        """
        tmpdir = mktempdir()
        filepath = joinpath(tmpdir, "input.tglf")
        try
            write(filepath, content)
            input_tglf = TurbulentTransport.load(InputTGLF(), filepath)
            @test input_tglf.USE_BPER == true
            @test input_tglf.USE_BPAR == false
        finally
            rm(tmpdir; force=true, recursive=true)
        end
    end

    @testset "load handles integer/float conversion" begin
        content = """
        NS = 3
        SAT_RULE = 2
        BETAE = 1.5e-3
        Q_LOC = 2.0
        """
        tmpdir = mktempdir()
        filepath = joinpath(tmpdir, "input.tglf")
        try
            write(filepath, content)
            input_tglf = TurbulentTransport.load(InputTGLF(), filepath)
            @test input_tglf.NS == 3
            @test input_tglf.SAT_RULE == 2
            @test input_tglf.BETAE == 0.0015
            @test input_tglf.Q_LOC == 2.0
        finally
            rm(tmpdir; force=true, recursive=true)
        end
    end
end
