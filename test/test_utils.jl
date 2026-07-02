using GACODE

@testset "Utility Functions" begin
    @testset "flux_solution constructor (4 args)" begin
        # 4-argument version has special ordering for backward compatibility
        sol = flux_solution(1.0, 2.0, 3.0, 4.0)

        @test sol isa GACODE.FluxSolution
        # For 4 args: Qe=3, Qi=4, Γe=1, Πi=2
        @test sol.ENERGY_FLUX_e == 3.0
        @test sol.ENERGY_FLUX_i == 4.0
        @test sol.PARTICLE_FLUX_e == 1.0
        @test sol.STRESS_TOR_i == 2.0
        @test isempty(sol.PARTICLE_FLUX_i)
    end

    @testset "flux_solution constructor (>4 args)" begin
        # 6-argument version: Γe, Γi_1, Γi_2, Πi, Qe, Qi
        sol = flux_solution(1.0, 2.0, 3.0, 4.0, 5.0, 6.0)

        @test sol isa GACODE.FluxSolution
        @test sol.ENERGY_FLUX_e == 5.0
        @test sol.ENERGY_FLUX_i == 6.0
        @test sol.PARTICLE_FLUX_e == 1.0
        @test sol.STRESS_TOR_i == 4.0
        @test sol.PARTICLE_FLUX_i == [2.0, 3.0]
    end

    @testset "save and reload" begin
        input_tglf = load_sample_input()

        # Save to temp file
        tmpdir = mktempdir()
        try
            save_path = joinpath(tmpdir, "saved_input.tglf")
            TurbulentTransport.save(input_tglf, save_path)

            # Reload
            reloaded = TurbulentTransport.load(InputTGLF(), save_path)

            # Check key fields match
            @test reloaded.NS == input_tglf.NS
            @test reloaded.SAT_RULE == input_tglf.SAT_RULE
            @test reloaded.BETAE ≈ input_tglf.BETAE
            @test reloaded.Q_LOC ≈ input_tglf.Q_LOC
        finally
            rm(tmpdir; force=true, recursive=true)
        end
    end

    @testset "load parse errors" begin
        mktempdir() do tmp
            # Unparseable integer field (NS is an Int field).
            p_int = joinpath(tmp, "bad_int.tglf")
            write(p_int, "NS=notanumber\n")
            @test_throws ErrorException TurbulentTransport.load(InputTGLF(), p_int)

            # Unparseable real field (BETAE is a Float64 field).
            p_real = joinpath(tmp, "bad_real.tglf")
            write(p_real, "BETAE=xyz\n")
            @test_throws ErrorException TurbulentTransport.load(InputTGLF(), p_real)

            # Neither key=value nor gen-style "value  key" -> invalid file error.
            p_fmt = joinpath(tmp, "bad_fmt.tglf")
            write(p_fmt, "just some text\n")
            @test_throws ErrorException TurbulentTransport.load(InputTGLF(), p_fmt)
        end
    end

    @testset "load gen-style (value  key) format" begin
        mktempdir() do tmp
            # gen-style: "<value><double-space><KEY>" for every line.
            p = joinpath(tmp, "input.tglf.gen")
            write(p, "3  NS\n3  SAT_RULE\n0.00362972  BETAE\n")
            loaded = TurbulentTransport.load(InputTGLF(), p)
            @test loaded.NS == 3
            @test loaded.SAT_RULE == 3
            @test loaded.BETAE ≈ 0.00362972
        end
    end

    @testset "parse_out_tglf_gbflux" begin
        # Sample output format from TGLF
        # Format: species values for [Gam, Q, Pi, S] × [elec, ion1, ion2, ...]
        sample_output = """
        1.5
        2.0
        3.0
        4.5
        5.0
        6.0
        7.5
        8.0
        9.0
        10.5
        11.0
        12.0
        """

        result = TurbulentTransport.parse_out_tglf_gbflux(sample_output)

        @test result isa Dict
        @test haskey(result, "Gam/Gam_GB_elec")
        @test haskey(result, "Q/Q_GB_elec")
        @test haskey(result, "Pi/Pi_GB_elec")
        @test haskey(result, "S/S_GB_elec")
    end

    @testset "compare_two_input_tglfs" begin
        input1 = load_sample_input()
        input2 = load_sample_input()

        # Modify input2 slightly
        input2.BETAE = input1.BETAE * 1.1

        diff_result = TurbulentTransport.compare_two_input_tglfs(input1, input2)

        @test diff_result isa InputTGLF
        # BETAE difference should be 10% of original
        @test abs(diff_result.BETAE) ≈ abs(input1.BETAE * 0.1) rtol=0.01
    end

    @testset "diff function" begin
        input1 = load_sample_input()
        input2 = load_sample_input()

        # Identical inputs should have no differences
        differences = TurbulentTransport.diff(input1, input2)
        @test isempty(differences)

        # Modify some fields
        input2.BETAE = 0.999
        input2.Q_LOC = 5.0
        input2.NS = 5

        differences = TurbulentTransport.diff(input1, input2)

        @test :BETAE in differences
        @test :Q_LOC in differences
        @test :NS in differences
        @test length(differences) == 3
    end

    @testset "scan function" begin
        input = load_sample_input()

        # Test with single parameter scan
        inputs = TurbulentTransport.scan(input; BETAE=[0.001, 0.002, 0.003])

        @test length(inputs) == 3
        @test all(x -> x isa InputTGLF, inputs)
        @test inputs[1].BETAE == 0.001
        @test inputs[2].BETAE == 0.002
        @test inputs[3].BETAE == 0.003

        # Other fields should remain unchanged
        @test all(x -> x.Q_LOC == input.Q_LOC, inputs)
        @test all(x -> x.NS == input.NS, inputs)
    end

    @testset "scan function - multiple parameters" begin
        input = load_sample_input()

        # Test with two parameters (cartesian product)
        inputs = TurbulentTransport.scan(input; BETAE=[0.001, 0.002], Q_LOC=[1.0, 2.0])

        # Should produce 2 × 2 = 4 combinations
        @test length(inputs) == 4

        # Check all combinations exist
        betae_values = [x.BETAE for x in inputs]
        q_values = [x.Q_LOC for x in inputs]

        @test 0.001 in betae_values
        @test 0.002 in betae_values
        @test 1.0 in q_values
        @test 2.0 in q_values
    end

    @testset "scan function - empty keywords" begin
        input = load_sample_input()

        # No keywords returns array with single copy
        inputs = TurbulentTransport.scan(input)

        @test length(inputs) == 1
        @test inputs[1].BETAE == input.BETAE
        # Should be a deepcopy, not same object
        @test inputs[1] !== input
    end

    @testset "save InputCGYRO" begin
        input_tglf = load_sample_input()
        ic = TurbulentTransport.tglf_to_cgyro(input_tglf)
        mktempdir() do tmpdir
            path = joinpath(tmpdir, "input.cgyro")
            TurbulentTransport.save(ic, path)
            @test isfile(path)
            content = read(path, String)
            @test contains(content, "N_SPECIES")
            @test contains(content, "N_FIELD")
        end
    end

    @testset "save InputQLGYRO" begin
        iq = TurbulentTransport.InputQLGYRO(NKY=16, KYGRID_MODEL=1)
        mktempdir() do tmpdir
            path = joinpath(tmpdir, "input.qlgyro")
            TurbulentTransport.save(iq, path)
            @test isfile(path)
            content = read(path, String)
            @test contains(content, "NKY=16")
            @test contains(content, "KYGRID_MODEL=1")
        end
    end
end
