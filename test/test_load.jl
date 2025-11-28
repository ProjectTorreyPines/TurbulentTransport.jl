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

    @testset "load handles input.tglf.gen format (space-separated)" begin
        # input.tglf.gen uses "value  fieldname" format with double spaces
        content = """
        3  NS
        2.0  Q_LOC
        0.003  BETAE
        .true.  USE_BPER
        """
        tmpdir = mktempdir()
        filepath = joinpath(tmpdir, "input.tglf.gen")
        try
            write(filepath, content)
            input_tglf = TurbulentTransport.load(InputTGLF(), filepath)
            @test input_tglf.NS == 3
            @test input_tglf.Q_LOC == 2.0
            @test input_tglf.BETAE == 0.003
            @test input_tglf.USE_BPER == true
        finally
            rm(tmpdir; force=true, recursive=true)
        end
    end

    @testset "load error on invalid format" begin
        # Invalid format: neither "=" nor "  " separated
        content = """
        NS:3
        Q_LOC:2.0
        """
        tmpdir = mktempdir()
        filepath = joinpath(tmpdir, "invalid.tglf")
        try
            write(filepath, content)
            @test_throws ErrorException TurbulentTransport.load(InputTGLF(), filepath)
        finally
            rm(tmpdir; force=true, recursive=true)
        end
    end

    @testset "load handles T/F boolean shorthand" begin
        content = """
        NS = 2
        USE_BPER = T
        USE_BPAR = F
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

    @testset "load handles UNITS string field" begin
        content = """
        NS = 2
        UNITS = CGYRO
        """
        tmpdir = mktempdir()
        filepath = joinpath(tmpdir, "input.tglf")
        try
            write(filepath, content)
            input_tglf = TurbulentTransport.load(InputTGLF(), filepath)
            @test input_tglf.UNITS == "CGYRO"
        finally
            rm(tmpdir; force=true, recursive=true)
        end
    end

    @testset "load ignores unknown fields gracefully" begin
        content = """
        NS = 2
        SOME_UNKNOWN_FIELD = 999
        Q_LOC = 1.5
        """
        tmpdir = mktempdir()
        filepath = joinpath(tmpdir, "input.tglf")
        try
            write(filepath, content)
            # Should not error; unknown fields are silently skipped
            input_tglf = TurbulentTransport.load(InputTGLF(), filepath)
            @test input_tglf.NS == 2
            @test input_tglf.Q_LOC == 1.5
        finally
            rm(tmpdir; force=true, recursive=true)
        end
    end
end
