# Tests for the `dd`-based input constructors that translate an IMAS `dd`
# (equilibrium + core_profiles) into TGLF / CGYRO / TGLF-EP inputs. These
# exercise the large evaluation functions in tglf.jl, cgyro.jl, and tglf_ep.jl
# using the committed `sample_dd.json` fixture, so no FUSE/scheduler is needed.

using TurbulentTransport: InputTGLF, InputCGYRO, InputTGLFEP, InputTGLFs

@testset "dd-based inputs" begin
    dd = load_sample_dd()
    rho = [0.3, 0.5, 0.7]

    @testset "InputTGLF(dd, rho)" begin
        its = InputTGLF(dd, rho, :sat0, false, true)
        @test its isa InputTGLFs
        @test length(its.tglfs) == length(rho)
        @test its[1] isa InputTGLF

        # Every radius gets a species count and finite core geometry/beta.
        for k in eachindex(rho)
            it = its[k]
            @test it.NS >= 2
            @test it.SAT_RULE == 0
            @test isfinite(it.Q_LOC)
            @test isfinite(it.RMIN_LOC)
            @test isfinite(it.BETAE)
            @test it.BETAE > 0
        end
    end

    @testset "InputTGLF(dd, rho) sat/em/lump variants" begin
        # SAT_RULE and electromagnetic flags propagate.
        its1 = InputTGLF(dd, rho, :sat1, true, true)
        @test all(==(1), its1.SAT_RULE)
        @test all(its1[k].USE_BPER == true for k in eachindex(rho))

        # lump_ions=false keeps individual ion species (>= lumped count).
        its_lump = InputTGLF(dd, rho, :sat0, false, true)
        its_full = InputTGLF(dd, rho, :sat0, false, false)
        @test its_full[1].NS >= its_lump[1].NS
    end

    @testset "InputTGLF(dd, gridpoint_cp)" begin
        gp = [3, 6, 9]
        its = InputTGLF(dd, gp, :sat0, false, true)
        @test its isa InputTGLFs
        @test length(its.tglfs) == length(gp)
    end

    @testset "InputCGYRO(dd, gridpoint)" begin
        ic = InputCGYRO(dd, 5, true)
        @test ic isa InputCGYRO
        @test !ismissing(ic.N_SPECIES)
        @test ic.N_SPECIES >= 2
        for f in (:RMIN, :RMAJ, :Q, :S, :KAPPA, :BETAE_UNIT, :NU_EE, :GAMMA_E, :GAMMA_P, :MACH)
            v = getproperty(ic, f)
            @test !ismissing(v)
            @test isfinite(v)
        end
        # Sign conventions are +/- 1.
        @test abs(ic.BTCCW) == 1
        @test abs(ic.IPCCW) == 1
    end

    @testset "InputCGYRO variants (lump/fast/MXH)" begin
        # Un-lumped, no fast ions.
        ic_nf = InputCGYRO(dd, 5, false; fast_ions=false)
        @test ic_nf.N_SPECIES >= 2

        # MXH multi-mode geometry path (fills SHAPE_* coefficients).
        ic_mxh = InputCGYRO(dd, 5, true; MXH_modes=3)
        @test !ismissing(ic_mxh.KAPPA)
        @test !ismissing(ic_mxh.SHAPE_COS0)
    end

    @testset "InputTGLFEP(dd, rho)" begin
        res = InputTGLFEP(dd, rho, :sat0, false, false)
        @test res isa Tuple
        its, extraEP = res
        @test its isa InputTGLFs
        @test extraEP isa Dict
        @test length(its.tglfs) == length(rho)
        # is_ep appends the fast energetic-particle species, so NS grows.
        @test its[1].NS >= 3

        # gridpoint variant returns the same shape.
        res_gp = InputTGLFEP(dd, [3, 6, 9], :sat0, false, false)
        @test res_gp[1] isa InputTGLFs
    end

    @testset "InputTGLFEP units regression (BETAE/XNUE/DEBYE match InputTGLF)" begin
        # Regression for a units bug where InputTGLF_EP fed its working-unit
        # variables (Te [keV], ne [1e19 m^-3], a [m]) into the cgs formulas
        # copied from tglf.jl (which expect eV, cm^-3, cm). BETAE came out
        # ~1e16 too small, silently disabling all electromagnetic physics in
        # the TJLFEP IMAS path (every Alfvenic scan returned "always stable").
        # The electron-scale scalars must agree with the standard path exactly.
        its_std = InputTGLF(dd, rho, :sat0, false, false)
        its_ep, _ = InputTGLFEP(dd, rho, :sat0, false, false)
        for k in eachindex(rho)
            @test its_ep[k].BETAE ≈ its_std[k].BETAE rtol = 1e-10
            @test its_ep[k].XNUE ≈ its_std[k].XNUE rtol = 1e-10
            @test its_ep[k].DEBYE ≈ its_std[k].DEBYE rtol = 1e-10
            # sanity: physically plausible magnitudes (the bug gave ~1e-19)
            @test 1e-5 < its_ep[k].BETAE < 1e-1
            @test its_ep[k].XNUE > 1e-4
        end
    end
end
