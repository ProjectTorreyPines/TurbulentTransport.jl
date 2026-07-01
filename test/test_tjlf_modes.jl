# Tests for src/tjlf.jl: turbulence-mode classification, the SAT Lorentzian
# spectrum builders, and radial-correlation-length reconstruction.
#
# The private classification / spectrum helpers are pure and are tested directly.
# `identify_modes` and `radial_correlation_length` are exercised end-to-end using
# TJLF (already a dependency) on the shared sample input.

using TurbulentTransport: ITG, TEM, KBM, ETG, MTM, TurbulenceMode
using TurbulentTransport: TJLFModeIdentification, TJLFFluctuationSpectra
using TurbulentTransport: identify_modes, radial_correlation_length
using TurbulentTransport: MODE_COLORS, MODE_LABELS

@testset "tjlf mode identification & spectra" begin
    @testset "TurbulenceMode enum metadata" begin
        modes = collect(instances(TurbulenceMode))
        @test Set(modes) == Set([ITG, TEM, KBM, ETG, MTM])
        for m in modes
            @test haskey(MODE_COLORS, m)
            @test haskey(MODE_LABELS, m)
        end
        @test MODE_LABELS[ITG] == "ITG"
        @test MODE_COLORS[MTM] == :red
    end

    @testset "_classify_mode_with_sign" begin
        # em_threshold=0.5, ion_electron_threshold=0.5, ky_etg=2.0, signetg=+1
        f = TurbulentTransport._classify_mode_with_sign
        @test f(1.0, 1.0, 1.0, 1.0, 1.0, 0.5, 0.5, 2.0) == MTM   # EM, freq·sign ≥ 0
        @test f(1.0, -1.0, 1.0, 1.0, 1.0, 0.5, 0.5, 2.0) == KBM  # EM, freq·sign < 0
        @test f(0.1, 1.0, 1.0, 1.0, 1.0, 0.5, 0.5, 2.0) == TEM   # ES, e-dir, high ie_ratio
        @test f(0.1, 1.0, 3.0, 0.1, 1.0, 0.5, 0.5, 2.0) == ETG   # ES, e-dir, low ie, high ky
        @test f(0.1, 1.0, 1.0, 0.1, 1.0, 0.5, 0.5, 2.0) == TEM   # ES, e-dir, low ie, low ky
        @test f(0.1, -1.0, 1.0, 0.1, 1.0, 0.5, 0.5, 2.0) == ITG  # ES, ion-dir
    end

    @testset "_classify_mode_no_sign" begin
        g = TurbulentTransport._classify_mode_no_sign
        @test g(1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 2.0) == KBM   # EM, emi_ratio > 0.25
        @test g(1.0, 1.0, 1.0, 0.1, 0.5, 0.5, 2.0) == MTM   # EM, emi_ratio ≤ 0.25
        @test g(0.1, 1.0, 2.0, 0.0, 0.5, 0.5, 2.0) == ITG   # ES, ie_ratio ≥ 1.5
        @test g(0.1, 1.0, 1.0, 0.0, 0.5, 0.5, 2.0) == TEM   # ES, 0.5 ≤ ie < 1.5
        @test g(0.1, 3.0, 0.1, 0.0, 0.5, 0.5, 2.0) == ETG   # ES, low ie, high ky
        @test g(0.1, 1.0, 0.1, 0.0, 0.5, 0.5, 2.0) == TEM   # ES, low ie, low ky
    end

    @testset "_sat_lorentz_shape" begin
        h = TurbulentTransport._sat_lorentz_shape
        @test h(0.0, 1.0, 1.0, 4) == 1.0        # peak at u = 0
        @test h(1.0, 1.0, 1.0, 4) < 1.0         # decays away from peak
        @test h(2.0, 1.0, 1.0, 4) < h(1.0, 1.0, 1.0, 4)
        @test h(-1.0, 1.0, 1.0, 4) ≈ h(1.0, 1.0, 1.0, 4)  # even in u
    end

    @testset "_build_phi2_from_sat_params" begin
        ky = [0.5]
        phinorm = reshape([2.0], 1, 1)     # (nky, nmodes)
        kx_width = [0.5]
        kx0_e = [0.0]

        # Explicit kx grid: with kx0_e = 0, phi2 at kx = 0 equals phinorm.
        kx_in = [-1.0, 0.0, 1.0]
        kxv, phi2 = TurbulentTransport._build_phi2_from_sat_params(
            phinorm, kx_width, kx0_e, 1.0, 1.0, 4, ky; kx=kx_in)
        @test kxv == kx_in
        @test size(phi2) == (3, 1, 1)
        @test phi2[2, 1, 1] ≈ 2.0            # peak equals phinorm at kx=0
        @test phi2[1, 1, 1] < phi2[2, 1, 1]  # decays off-peak
        @test all(isfinite, phi2)

        # Auto-sized kx grid honors n_kx.
        kxv2, phi2b = TurbulentTransport._build_phi2_from_sat_params(
            phinorm, kx_width, kx0_e, 1.0, 1.0, 4, ky; n_kx=41, kx_max_sigma=3.0)
        @test length(kxv2) == 41
        @test size(phi2b) == (41, 1, 1)

        # Dimension-mismatch assertions.
        @test_throws AssertionError TurbulentTransport._build_phi2_from_sat_params(
            reshape([1.0, 2.0], 2, 1), kx_width, kx0_e, 1.0, 1.0, 4, ky; kx=kx_in)
    end

    @testset "_density_spectrum_from_weights" begin
        phi2 = reshape(Float64[1, 2, 3, 4], 2, 2, 1)  # (nkx=2, nky=2, nmodes=1)
        @test TurbulentTransport._density_spectrum_from_weights(nothing, phi2, 1, 2, Float64) === nothing

        # dw shape (ns, nmodes, nky) = (2, 1, 2)
        dw = reshape(Float64[10, 20, 30, 40], 2, 1, 2)
        d2 = TurbulentTransport._density_spectrum_from_weights(dw, phi2, 1, 2, Float64)
        @test size(d2) == (2, 2, 1, 2)  # (nkx, nky, nmodes, ns)
        @test d2[1, 1, 1, 1] ≈ dw[1, 1, 1] * phi2[1, 1, 1]
        @test d2[2, 2, 1, 2] ≈ dw[2, 1, 2] * phi2[2, 2, 1]

        # Wrong ndims / sizes error.
        @test_throws AssertionError TurbulentTransport._density_spectrum_from_weights(
            reshape(Float64[1, 2], 2, 1), phi2, 1, 2, Float64)
    end

    @testset "TJLFFluctuationSpectra show + radial_correlation_length" begin
        ky = [0.3, 0.6]
        phinorm = reshape([1.0, 0.8], 2, 1)
        kx_width = [0.5, 0.4]
        kx0_e = [0.0, 0.1]
        kxv, phi2 = TurbulentTransport._build_phi2_from_sat_params(
            phinorm, kx_width, kx0_e, 1.0, 1.0, 4, ky; n_kx=128, kx_max_sigma=6.0)

        fs = TJLFFluctuationSpectra{Float64}(
            kxv, ky, phi2, copy(phinorm), kx_width, kx0_e, 1.0, 1.0, 4, 2, nothing)

        io = IOBuffer()
        show(io, MIME"text/plain"(), fs)
        s = String(take!(io))
        @test occursin("TJLFFluctuationSpectra", s)
        @test occursin("SAT_RULE=2", s)
        @test occursin("not computed", s)  # density2 === nothing branch

        # With density2 populated, show reports the 4-D array.
        dw = ones(Float64, 2, 1, length(ky))  # (ns, nmodes, nky)
        d2 = TurbulentTransport._density_spectrum_from_weights(dw, phi2, 1, length(ky), Float64)
        fs2 = TJLFFluctuationSpectra{Float64}(
            kxv, ky, phi2, copy(phinorm), kx_width, kx0_e, 1.0, 1.0, 4, 2, d2)
        io2 = IOBuffer()
        show(io2, MIME"text/plain"(), fs2)
        @test occursin("density2: Array", String(take!(io2)))

        # radial_correlation_length: all three width methods.
        for method in (:hwhm, :efold, :integral)
            res = radial_correlation_length(fs; method=method)
            @test length(res.Lr) == length(ky)
            @test all(isfinite, res.Lr)
            @test all(>=(0), res.Lr)
            @test isfinite(res.L_avg)
        end

        # mode slice path + error paths.
        res_modes = radial_correlation_length(fs; modes=[1])
        @test length(res_modes.Lr) == length(ky)
        @test_throws ArgumentError radial_correlation_length(fs; dr_max=-1.0)
    end

    @testset "identify_modes end-to-end (TJLF)" begin
        input_tglf = load_sample_input()
        mid = identify_modes(input_tglf)

        @test mid isa TJLFModeIdentification
        @test mid.dominant_mode isa TurbulenceMode
        @test 0.0 <= mid.dominant_mode_fraction <= 1.0
        @test length(mid.mode_per_ky) == length(mid.ky_spectrum)
        @test all(m -> m isa TurbulenceMode, mid.mode_per_ky)
        # energy_flux_per_mode has an entry for every mode.
        for m in instances(TurbulenceMode)
            @test haskey(mid.energy_flux_per_mode, m)
        end
        @test mid.flux_solution isa GACODE.FluxSolution

        io = IOBuffer()
        show(io, MIME"text/plain"(), mid)
        s = String(take!(io))
        @test occursin("TJLFModeIdentification", s)
        @test occursin("Dominant mode", s)
    end
end
