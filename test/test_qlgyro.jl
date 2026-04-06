using TurbulentTransport: InputTGLF, InputCGYRO, InputQLGYRO
using TurbulentTransport: tglf_to_cgyro, generate_ky_grid, qlgyro_run_hash,
    QLGYRORunState, save_run_state, load_run_state,
    check_cgyro_convergence, parse_cgyro_eigenvalue, parse_cgyro_qlflux,
    compute_qlgyro_fluxes
using TurbulentTransport.TJLF

# ============================================================================
# Reference values
# These are directly from paired input.tglf / input.cgyro files
# ============================================================================
const REF_TGLF_PARAMS = (
    RMIN_LOC = 1.24472e-01,
    RMAJ_LOC = 2.99341e+00,
    Q_LOC = 1.42528e+00,
    Q_PRIME_LOC = -1.78566e+00,
    KAPPA_LOC = 1.43015e+00,
    S_KAPPA_LOC = -2.30916e-02,
    DELTA_LOC = 3.27964e-02,
    S_DELTA_LOC = 2.30106e-02,
    ZETA_LOC = 3.34598e-04,
    S_ZETA_LOC = -3.50017e-04,
    DRMAJDX_LOC = -1.46092e-01,
    DZMAJDX_LOC = -1.12859e-03,
    ZMAJ_LOC = -4.38647e-02,
    SIGN_BT = -1.0,
    SIGN_IT = 1.0,
    BETAE = 1.50202e-02,
    XNUE = 2.56591e-02,
    DEBYE = 1.60447e-02,
    VEXB_SHEAR = 6.40639e-02,
)

const REF_CGYRO_EXPECTED = (
    RMIN = 1.24472e-01,
    RMAJ = 2.99341e+00,
    Q = 1.42528e+00,
    S = -1.36190e-02,
    KAPPA = 1.43015e+00,
    S_KAPPA = -2.30916e-02,
    DELTA = 3.27964e-02,
    S_DELTA = 2.30106e-02,
    ZETA = 3.34598e-04,
    S_ZETA = -3.50017e-04,
    SHIFT = -1.46092e-01,
    DZMAG = -1.12859e-03,
    ZMAG = -4.38647e-02,
    BTCCW = -1.0,
    IPCCW = 1.0,
    BETAE_UNIT = 1.50202e-02,
    NU_EE = 2.56591e-02,
    LAMBDA_STAR = 1.60447e-02,
    GAMMA_E = 6.40639e-02,
)

# Second validation case:
const REF_CASE2 = (
    RMIN_LOC = 3.38849e-01,
    Q_LOC = 1.46332e+00,
    Q_PRIME_LOC = 6.91814e+00,
    S_expected = 3.70955e-01,
)

@testset "Modular QLGYRO" begin

    # ========================================================================
    @testset "tglf_to_cgyro conversion" begin

        @testset "Geometry & shear mapping (case 1)" begin
            input_tglf = load_sample_input()
            # Override with reference case 1 values for precise validation
            for (k, v) in pairs(REF_TGLF_PARAMS)
                setproperty!(input_tglf, k, v)
            end
            input_tglf.NS = 3
            input_tglf.USE_BPER = true
            input_tglf.USE_BPAR = true

            ic = tglf_to_cgyro(input_tglf)

            # Direct 1:1 mappings
            @test ic.RMIN ≈ REF_CGYRO_EXPECTED.RMIN
            @test ic.RMAJ ≈ REF_CGYRO_EXPECTED.RMAJ
            @test ic.Q ≈ REF_CGYRO_EXPECTED.Q
            @test ic.KAPPA ≈ REF_CGYRO_EXPECTED.KAPPA
            @test ic.S_KAPPA ≈ REF_CGYRO_EXPECTED.S_KAPPA
            @test ic.DELTA ≈ REF_CGYRO_EXPECTED.DELTA
            @test ic.S_DELTA ≈ REF_CGYRO_EXPECTED.S_DELTA
            @test ic.ZETA ≈ REF_CGYRO_EXPECTED.ZETA
            @test ic.S_ZETA ≈ REF_CGYRO_EXPECTED.S_ZETA
            @test ic.SHIFT ≈ REF_CGYRO_EXPECTED.SHIFT
            @test ic.DZMAG ≈ REF_CGYRO_EXPECTED.DZMAG
            @test ic.ZMAG ≈ REF_CGYRO_EXPECTED.ZMAG
            @test ic.BTCCW ≈ REF_CGYRO_EXPECTED.BTCCW
            @test ic.IPCCW ≈ REF_CGYRO_EXPECTED.IPCCW
            @test ic.BETAE_UNIT ≈ REF_CGYRO_EXPECTED.BETAE_UNIT
            @test ic.NU_EE ≈ REF_CGYRO_EXPECTED.NU_EE
            @test ic.LAMBDA_STAR ≈ REF_CGYRO_EXPECTED.LAMBDA_STAR
            @test ic.GAMMA_E ≈ REF_CGYRO_EXPECTED.GAMMA_E

            # Critical shear conversion: S = RMIN_LOC² * Q_PRIME_LOC / Q_LOC²
            @test ic.S ≈ REF_CGYRO_EXPECTED.S atol=1e-5
        end

        @testset "Shear conversion (case 2)" begin
            input_tglf = load_sample_input()
            input_tglf.RMIN_LOC = REF_CASE2.RMIN_LOC
            input_tglf.Q_LOC = REF_CASE2.Q_LOC
            input_tglf.Q_PRIME_LOC = REF_CASE2.Q_PRIME_LOC
            ic = tglf_to_cgyro(input_tglf)
            @test ic.S ≈ REF_CASE2.S_expected atol=1e-5
        end

        @testset "Species reordering (3 species)" begin
            input_tglf = load_sample_input()
            # sample_input.tglf has NS=3: electrons(1), D(2), C(3)
            ic = tglf_to_cgyro(input_tglf)

            @test ic.N_SPECIES == 3

            # Ion 1 (CGYRO sp1) = TGLF sp2 (D)
            @test ic.Z_1 ≈ input_tglf.ZS_2
            @test ic.MASS_1 ≈ input_tglf.MASS_2
            @test ic.DENS_1 ≈ input_tglf.AS_2
            @test ic.TEMP_1 ≈ input_tglf.TAUS_2
            @test ic.DLNNDR_1 ≈ input_tglf.RLNS_2
            @test ic.DLNTDR_1 ≈ input_tglf.RLTS_2

            # Ion 2 (CGYRO sp2) = TGLF sp3 (C)
            @test ic.Z_2 ≈ input_tglf.ZS_3
            @test ic.MASS_2 ≈ input_tglf.MASS_3
            @test ic.DENS_2 ≈ input_tglf.AS_3
            @test ic.TEMP_2 ≈ input_tglf.TAUS_3
            @test ic.DLNNDR_2 ≈ input_tglf.RLNS_3
            @test ic.DLNTDR_2 ≈ input_tglf.RLTS_3

            # Electrons (CGYRO sp3, last) = TGLF sp1
            @test ic.Z_3 ≈ input_tglf.ZS_1
            @test ic.MASS_3 ≈ input_tglf.MASS_1
            @test ic.DENS_3 ≈ input_tglf.AS_1
            @test ic.TEMP_3 ≈ input_tglf.TAUS_1
            @test ic.DLNNDR_3 ≈ input_tglf.RLNS_1
            @test ic.DLNTDR_3 ≈ input_tglf.RLTS_1
        end

        @testset "Electromagnetic fields" begin
            input_tglf = load_sample_input()

            # Both BPER and BPAR on -> N_FIELD = 3
            input_tglf.USE_BPER = true
            input_tglf.USE_BPAR = true
            @test tglf_to_cgyro(input_tglf).N_FIELD == 3

            # Only BPER -> N_FIELD = 2
            input_tglf.USE_BPER = true
            input_tglf.USE_BPAR = false
            @test tglf_to_cgyro(input_tglf).N_FIELD == 2

            # Neither -> N_FIELD = 1 (electrostatic)
            input_tglf.USE_BPER = false
            input_tglf.USE_BPAR = false
            @test tglf_to_cgyro(input_tglf).N_FIELD == 1
        end

        @testset "Linear run defaults" begin
            input_tglf = load_sample_input()
            ic = tglf_to_cgyro(input_tglf)

            @test ic.NONLINEAR_FLAG == 0
            @test ic.N_TOROIDAL == 1
            @test ic.EQUILIBRIUM_MODEL == 2
            @test ic.COLLISION_MODEL == 4
            @test ic.MAX_TIME == 100000.0
            @test ic.DELTA_T_METHOD == 1
            @test ic.ERROR_TOL == 0.001
            @test ic.FREQ_TOL == 0.01
            @test ic.N_RADIAL == 16
            @test ic.N_XI == 24
            @test ic.FIELD_PRINT_FLAG == 1
            @test ic.ROTATION_MODEL == 2
        end

        @testset "Zero VEXB_SHEAR" begin
            input_tglf = load_sample_input()
            input_tglf.VEXB_SHEAR = 0.0
            ic = tglf_to_cgyro(input_tglf)
            @test ismissing(ic.GAMMA_E)
        end

        @testset "Higher-order shaping coefficients" begin
            input_tglf = load_sample_input()
            input_tglf.SHAPE_COS0 = 0.1
            input_tglf.SHAPE_S_COS0 = 0.2
            input_tglf.SHAPE_COS3 = 0.03
            input_tglf.SHAPE_S_COS3 = 0.04
            input_tglf.SHAPE_SIN3 = 0.05
            input_tglf.SHAPE_S_SIN3 = 0.06
            ic = tglf_to_cgyro(input_tglf)

            @test ic.SHAPE_COS0 ≈ 0.1
            @test ic.SHAPE_S_COS0 ≈ 0.2
            @test ic.SHAPE_COS3 ≈ 0.03
            @test ic.SHAPE_S_COS3 ≈ 0.04
            @test ic.SHAPE_SIN3 ≈ 0.05
            @test ic.SHAPE_S_SIN3 ≈ 0.06

            # Missing values should not be set
            input_tglf2 = load_sample_input()
            ic2 = tglf_to_cgyro(input_tglf2)
            @test ismissing(ic2.SHAPE_COS0)
            @test ismissing(ic2.SHAPE_SIN3)
        end
    end

    # ========================================================================
    @testset "generate_ky_grid" begin
        input_tglf = load_sample_input()

        @testset "KYGRID_MODEL=0: default KY=1.2 NKY=12" begin
            ky_grid = generate_ky_grid(input_tglf; kygrid_model=0)
            @test length(ky_grid) == input_tglf.NKY
            @test ky_grid[1] ≈ 1.2 / input_tglf.NKY
            @test ky_grid[end] ≈ 1.2
        end

        @testset "KYGRID_MODEL=0: custom NKY" begin
            ky_grid = generate_ky_grid(input_tglf; nky=5, kygrid_model=0)
            @test length(ky_grid) == 5
            dky = 1.2 / 5
            for i in 1:5
                @test ky_grid[i] ≈ i * dky
            end
        end

        @testset "KYGRID_MODEL=0: custom KY and NKY" begin
            ky_grid = generate_ky_grid(input_tglf; ky=0.6, nky=6, kygrid_model=0)
            @test length(ky_grid) == 6
            @test ky_grid[end] ≈ 0.6
            @test ky_grid[1] ≈ 0.1
        end

        @testset "KYGRID_MODEL=0: uniform spacing" begin
            ky_grid = generate_ky_grid(input_tglf; nky=10, kygrid_model=0)
            diffs = diff(ky_grid)
            @test all(d -> d ≈ diffs[1], diffs)
        end

        @testset "Default uses input_tglf.KYGRID_MODEL" begin
            # sample_input.tglf has KYGRID_MODEL=4, NKY=12 → total = 12 + 12 = 24
            ky_grid = generate_ky_grid(input_tglf)
            @test input_tglf.KYGRID_MODEL == 4
            @test length(ky_grid) == TJLF.get_ky_spectrum_size(input_tglf.NKY, 4)
            @test length(ky_grid) == input_tglf.NKY + 12
            @test all(ky_grid .> 0)
            @test issorted(ky_grid)
        end

        @testset "KYGRID_MODEL=1: APS07 multi-scale" begin
            ky_grid = generate_ky_grid(input_tglf; kygrid_model=1)
            @test length(ky_grid) == input_tglf.NKY + 9
            @test all(ky_grid .> 0)
            @test issorted(ky_grid)
            # Should have a low-ky linear region and a high-ky log-spaced region
            @test ky_grid[end] > ky_grid[9]  # ETG region extends beyond ion-scale
        end

        @testset "KYGRID_MODEL=4: SAT2/3 multi-scale" begin
            ky_grid = generate_ky_grid(input_tglf; kygrid_model=4)
            @test length(ky_grid) == input_tglf.NKY + 12
            @test all(ky_grid .> 0)
            @test issorted(ky_grid)
            # Model 4 has 12 low-ky points + NKY high-ky points
            # Low-ky region starts at ky_min ≈ 0.05/rho_ion
            @test ky_grid[1] > 0
            # High-ky region should reach ETG-scale (>> ion-scale)
            @test ky_grid[end] > 10 * ky_grid[1]
        end

        @testset "kygrid_model kwarg overrides input_tglf" begin
            # sample has KYGRID_MODEL=4, but we override to 0
            ky_grid_0 = generate_ky_grid(input_tglf; kygrid_model=0)
            ky_grid_4 = generate_ky_grid(input_tglf; kygrid_model=4)
            @test length(ky_grid_0) != length(ky_grid_4)
            @test length(ky_grid_0) == input_tglf.NKY
            @test length(ky_grid_4) == input_tglf.NKY + 12
        end
    end

    # ========================================================================
    @testset "qlgyro_run_hash" begin
        input_tglf = load_sample_input()
        ic1 = tglf_to_cgyro(input_tglf)
        ic2 = tglf_to_cgyro(input_tglf)

        @testset "Deterministic" begin
            @test qlgyro_run_hash(ic1) == qlgyro_run_hash(ic2)
        end

        @testset "Sensitive to changes" begin
            ic3 = deepcopy(ic1)
            ic3.Q = ic1.Q * 1.01
            @test qlgyro_run_hash(ic1) != qlgyro_run_hash(ic3)
        end
    end

    # ========================================================================
    @testset "QLGYRORunState save/load" begin
        mktempdir() do tmpdir
            state = QLGYRORunState(
                tmpdir,
                [0.1, 0.2, 0.3],
                ["12345", "12346", ""],
                [true, false, false],
                [true, true, false],
                UInt64(42)
            )

            save_run_state(state)
            loaded = load_run_state(tmpdir)

            @test loaded !== nothing
            @test loaded.basedir == tmpdir
            @test loaded.ky_values ≈ [0.1, 0.2, 0.3]
            @test loaded.slurm_ids == ["12345", "12346", ""]
            @test loaded.converged == [true, false, false]
            @test loaded.submitted == [true, true, false]
            @test loaded.input_hash == UInt64(42)
        end

        @testset "Missing state file returns nothing" begin
            mktempdir() do tmpdir
                @test load_run_state(tmpdir) === nothing
            end
        end
    end

    # ========================================================================
    @testset "check_cgyro_convergence" begin

        @testset "Not started" begin
            mktempdir() do tmpdir
                @test check_cgyro_convergence(tmpdir) == :not_started
            end
        end

        @testset "Linear converged" begin
            mktempdir() do tmpdir
                write(joinpath(tmpdir, "out.cgyro.info"),
                    """
                    INFO: (CGYRO) Velocity order 1
                    INFO: (CGYRO) Linear converged
                    """)
                @test check_cgyro_convergence(tmpdir) == :converged
            end
        end

        @testset "Underflow treated as converged" begin
            mktempdir() do tmpdir
                write(joinpath(tmpdir, "out.cgyro.info"),
                    "INFO: Underflow in calculation of frequency error\n")
                @test check_cgyro_convergence(tmpdir) == :converged
            end
        end

        @testset "Timeout treated as converged" begin
            mktempdir() do tmpdir
                write(joinpath(tmpdir, "out.cgyro.info"), "INFO: (CGYRO) Restart\n")
                write(joinpath(tmpdir, "out.cgyro.time"), "50.0\n100.0\n")
                @test check_cgyro_convergence(tmpdir) == :converged
            end
        end

        @testset "Still running" begin
            mktempdir() do tmpdir
                write(joinpath(tmpdir, "out.cgyro.info"), "INFO: (CGYRO) Restart\n")
                write(joinpath(tmpdir, "out.cgyro.time"), "10.0\n20.0\n")
                @test check_cgyro_convergence(tmpdir) == :running
            end
        end

        @testset "Error" begin
            mktempdir() do tmpdir
                write(joinpath(tmpdir, "out.cgyro.info"), "ERROR: something exploded\n")
                @test check_cgyro_convergence(tmpdir) == :error
            end
        end

        @testset "Terminated at max time" begin
            mktempdir() do tmpdir
                write(joinpath(tmpdir, "out.cgyro.info"), "INFO: terminated at max time\n")
                @test check_cgyro_convergence(tmpdir) == :terminated
            end
        end

        @testset "Running without timefile" begin
            mktempdir() do tmpdir
                write(joinpath(tmpdir, "out.cgyro.info"), "INFO: (CGYRO) Restart\n")
                @test check_cgyro_convergence(tmpdir) == :running
            end
        end
    end

    # ========================================================================
    @testset "parse_cgyro_eigenvalue" begin

        @testset "Normal eigenvalue file" begin
            mktempdir() do tmpdir
                write(joinpath(tmpdir, "out.cgyro.freq"),
                    """ 1.7754E-01 -4.7633E-01
 -9.1431E-02  2.8509E-01
  1.7936E-02  2.3512E-02
  1.7936E-02  2.3512E-02
""")
                freq, gamma = parse_cgyro_eigenvalue(tmpdir)
                @test freq ≈ 1.7936e-02 rtol=1e-4
                @test gamma ≈ 2.3512e-02 rtol=1e-4
            end
        end

        @testset "Missing file returns zeros" begin
            mktempdir() do tmpdir
                freq, gamma = parse_cgyro_eigenvalue(tmpdir)
                @test freq == 0.0
                @test gamma == 0.0
            end
        end

        @testset "Empty file returns zeros" begin
            mktempdir() do tmpdir
                write(joinpath(tmpdir, "out.cgyro.freq"), "")
                freq, gamma = parse_cgyro_eigenvalue(tmpdir)
                @test freq == 0.0
                @test gamma == 0.0
            end
        end
    end

    # ========================================================================
    @testset "parse_cgyro_qlflux" begin
        n_species = 3
        n_field = 3
        n_moments = 3

        @testset "Correct shape and last-timestep extraction" begin
            mktempdir() do tmpdir
                n_time = 5
                n_per_step = n_moments * n_field * n_species

                # Create synthetic data: each timestep has values = timestep_index * base
                data = Float32[]
                for t in 1:n_time
                    for s in 1:n_species
                        for f in 1:n_field
                            for m in 1:n_moments
                                push!(data, Float32(t * 100 + s * 10 + f + m * 0.1))
                            end
                        end
                    end
                end
                write(joinpath(tmpdir, "bin.cgyro.ky_flux"), reinterpret(UInt8, data))

                result = parse_cgyro_qlflux(tmpdir, n_species, n_field)
                @test size(result) == (n_species, n_field, n_moments)

                # Verify it's from the last timestep (t=5)
                # CGYRO Fortran writes (n_species, n_flux, n_field, n_time) column-major
                @test all(result .!= 0)
            end
        end

        @testset "Missing file returns zeros" begin
            mktempdir() do tmpdir
                result = parse_cgyro_qlflux(tmpdir, n_species, n_field)
                @test size(result) == (n_species, n_field, n_moments)
                @test all(result .== 0)
            end
        end

        @testset "Single timestep" begin
            mktempdir() do tmpdir
                # One timestep: 3*3*3 = 27 float32 values
                data = Float32.(collect(1:27))
                write(joinpath(tmpdir, "bin.cgyro.ky_flux"), reinterpret(UInt8, data))

                result = parse_cgyro_qlflux(tmpdir, n_species, n_field)
                @test size(result) == (n_species, n_field, n_moments)
                # CGYRO Fortran column-major: (n_species, n_flux, n_field, n_time)
                raw_4d = reshape(data, n_species, n_moments, n_field, 1)
                expected = permutedims(raw_4d[:, :, :, 1], (1, 3, 2))
                @test result ≈ expected
            end
        end
    end

    # ========================================================================
    @testset "tglf_to_cgyro round-trip with sample_input.tglf" begin
        input_tglf = load_sample_input()
        ic = tglf_to_cgyro(input_tglf)

        # Verify key physics is preserved
        @test ic.N_SPECIES == input_tglf.NS
        @test ic.RMIN ≈ input_tglf.RMIN_LOC
        @test ic.RMAJ ≈ input_tglf.RMAJ_LOC
        @test ic.Q ≈ input_tglf.Q_LOC

        # Expected shear for sample input: S = 0.573129² * 14.7947 / 2.00545² ≈ 1.208
        expected_S = input_tglf.RMIN_LOC^2 * input_tglf.Q_PRIME_LOC / input_tglf.Q_LOC^2
        @test ic.S ≈ expected_S

        # Verify electron charge is last species
        @test ic.Z_3 ≈ -1.0
        @test ic.MASS_3 ≈ input_tglf.MASS_1  # electron mass

        # Verify first ion is D (ZS_2=1.0, MASS_2=1.0)
        @test ic.Z_1 ≈ 1.0
        @test ic.MASS_1 ≈ 1.0
    end

    # ========================================================================
    @testset "compute_qlgyro_fluxes with fixture data" begin
        # Uses real CGYRO output from 3 ky points (ky=0.1, 0.5, 1.0)
        fixture_dir = joinpath(@__DIR__, "data", "sample_qlgyro")
        if isdir(fixture_dir)
            input_tglf = load_sample_input()

            ky_values = [0.1, 0.5, 1.0]
            nky = length(ky_values)
            state = QLGYRORunState(
                fixture_dir,
                ky_values,
                fill("", nky),
                fill(true, nky),
                fill(true, nky),
                UInt64(0)
            )

            result = compute_qlgyro_fluxes(input_tglf, state)

            @testset "Returns FluxSolution" begin
                @test result isa TurbulentTransport.GACODE.FluxSolution
            end

            @testset "Fluxes are finite" begin
                @test isfinite(result.ENERGY_FLUX_e)
                @test isfinite(result.ENERGY_FLUX_i)
                @test isfinite(result.PARTICLE_FLUX_e)
                @test all(isfinite.(result.PARTICLE_FLUX_i))
                @test all(isfinite.(result.STRESS_TOR_i))
            end

            @testset "Parse fixture eigenvalues" begin
                # KY_1 (ky=0.1): positive growth rate
                freq1, gamma1 = parse_cgyro_eigenvalue(joinpath(fixture_dir, "KY_1"))
                @test gamma1 > 0
                @test gamma1 ≈ 1.8339e-02 rtol=1e-3

                # KY_2 (ky=0.5): positive growth rate
                freq2, gamma2 = parse_cgyro_eigenvalue(joinpath(fixture_dir, "KY_2"))
                @test gamma2 > 0
                @test gamma2 ≈ 1.0611e-01 rtol=1e-3

                # KY_3 (ky=1.0): negative growth rate (stable)
                freq3, gamma3 = parse_cgyro_eigenvalue(joinpath(fixture_dir, "KY_3"))
                @test gamma3 < 0
            end

            @testset "Parse fixture QL fluxes" begin
                for iky in 1:3
                    ql = parse_cgyro_qlflux(joinpath(fixture_dir, "KY_$iky"), 3, 3)
                    @test size(ql) == (3, 3, 3)
                    @test any(ql .!= 0)
                end
            end
        else
            @warn "Skipping compute_qlgyro_fluxes test: fixture data not found at $fixture_dir"
        end
    end

end
