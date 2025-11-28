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

        diff = TurbulentTransport.compare_two_input_tglfs(input1, input2)

        @test diff isa InputTGLF
        # BETAE difference should be 10% of original
        @test abs(diff.BETAE) ≈ abs(input1.BETAE * 0.1) rtol=0.01
    end
end
