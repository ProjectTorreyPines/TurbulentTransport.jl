# Smoke tests for `qlnn_fluctuation_spectra`.
#
# What we cover:
#   1. Shape contract: `fs.kx`, `fs.ky`, `fs.phi2` have the documented sizes
#      (so downstream synth-diag code that consumes `bands_meta[k].fs_*`
#      keeps working).
#   2. Finiteness of every entry — NaN / Inf would silently propagate
#      through the kx-Lorentzian build into the rasterized turbulence frames.
#   3. Comparison with the reference TJLF path: `qlnn_fluctuation_spectra`
#      reuses the same `intensity_sat` saturation rule, but with QLNN-
#      predicted γ/ω/QL_weights upstream. We do NOT expect numerical agreement
#      (the NN is approximating γ/ω/QL), only that both paths return a
#      `TJLFFluctuationSpectra{Float64}` with positive `phi2` and the same
#      ky grid (TJLF's `get_ky_spectrum` is the same on both paths).

import TJLF

const QLNN_FS_BUNDLE_NAME = "QLNN"
const QLNN_FS_BUNDLE_DIR  = joinpath(dirname(@__DIR__), "models", QLNN_FS_BUNDLE_NAME)

if !isdir(QLNN_FS_BUNDLE_DIR)
    @info "Skipping qlnn_fluctuation_spectra tests; bundle directory not found: $QLNN_FS_BUNDLE_DIR"
else
    @testset "qlnn_fluctuation_spectra smoke tests" begin
        # Load the same fixture the rest of the QLNN tests use. The TJLF
        # reference path requires SAT_RULE ∈ {1,2,3} and ALPHA_QUENCH=0.
        input_tglf = load_sample_input()
        input_tjlf_qlnn = InputTJLF{Float64}(input_tglf)
        input_tjlf_tjlf = InputTJLF{Float64}(input_tglf)
        for it in (input_tjlf_qlnn, input_tjlf_tjlf)
            it.ALPHA_QUENCH = 0.0
            if !(it.SAT_RULE in (1, 2, 3))
                it.SAT_RULE = 2
            end
        end

        # --- QLNN path -----------------------------------------------------
        out_q = TurbulentTransport.qlnn_fluctuation_spectra(input_tjlf_qlnn;
                                                            bundle_name=QLNN_FS_BUNDLE_NAME,
                                                            n_kx=33, kx_max_sigma=2.5)
        # Single-input convenience wrapper returns a NamedTuple.
        @test out_q isa NamedTuple
        @test haskey(out_q, :fs)
        @test haskey(out_q, :eigenvalue)

        fs_q = out_q.fs
        eig_q = out_q.eigenvalue
        nky_q = length(fs_q.ky)
        nkx_q = length(fs_q.kx)

        @testset "shape contract" begin
            @test nky_q > 0
            @test nkx_q == 33
            # NMODES=1 in the QLNN path.
            @test size(fs_q.phi2) == (nkx_q, nky_q, 1)
            @test size(eig_q)     == (1, nky_q, 2)
            @test length(fs_q.kx_width) == nky_q
            @test length(fs_q.kx0_e)    == nky_q
            @test size(fs_q.phi2_ky)    == (nky_q, 1)
        end

        @testset "finiteness" begin
            @test all(isfinite, fs_q.kx)
            @test all(isfinite, fs_q.ky)
            @test all(isfinite, fs_q.phi2)
            @test all(>=(0), fs_q.phi2)        # |φ|² ≥ 0
            @test all(isfinite, eig_q)
            @test all(isfinite, fs_q.kx_width)
            @test all(>(0), fs_q.kx_width)     # widths must be positive
            @test all(isfinite, fs_q.phi2_ky)
        end

        @testset "ky grid matches the TJLF path" begin
            # `get_ky_spectrum` is the same on both paths, so the ky grid
            # should be byte-identical (same SaturationParameters, same
            # KYGRID_MODEL).
            fs_t = TurbulentTransport.fluctuation_spectra(input_tjlf_tjlf;
                                                          n_kx=33, kx_max_sigma=2.5)
            @test length(fs_t.ky) == nky_q
            @test fs_t.ky ≈ fs_q.ky atol=1e-10
            # kx/ky dimensions match; mode count may differ (QLNN=1, TJLF≥1).
            @test size(fs_t.phi2)[1:2] == size(fs_q.phi2)[1:2]
            # SAT_RULE flows through the saturation rule; QLNN uses the same
            # `intensity_sat`, so the recorded sat_rule must match.
            @test fs_t.sat_rule == fs_q.sat_rule
        end

        @testset "Vector wrapper returns one entry per radial point" begin
            its = [InputTJLF{Float64}(input_tglf) for _ in 1:3]
            for it in its
                it.ALPHA_QUENCH = 0.0
                it.SAT_RULE = 2
            end
            outs = TurbulentTransport.qlnn_fluctuation_spectra(its;
                                                               bundle_name=QLNN_FS_BUNDLE_NAME,
                                                               n_kx=21, kx_max_sigma=2.0)
            @test length(outs) == 3
            for o in outs
                @test o.fs isa TurbulentTransport.TJLFFluctuationSpectra
                @test all(isfinite, o.fs.phi2)
                @test size(o.eigenvalue, 1) == 1
                @test size(o.eigenvalue, 3) == 2
            end
        end

        @testset "WIDTH_SPECTRUM is populated by the helper" begin
            # qlnn_fluctuation_spectra mutates input_tjlf.WIDTH_SPECTRUM in
            # place before calling intensity_sat — either via the chained
            # width regressor (when `width_regressor.bson` is present in the
            # bundle dir) or via the constant-WIDTH fallback.
            it = InputTJLF{Float64}(input_tglf)
            it.ALPHA_QUENCH = 0.0
            it.SAT_RULE = 2
            _ = TurbulentTransport.qlnn_fluctuation_spectra(it;
                                                            bundle_name=QLNN_FS_BUNDLE_NAME,
                                                            n_kx=21, kx_max_sigma=2.0)
            @test !(it.WIDTH_SPECTRUM === missing)
            @test length(it.WIDTH_SPECTRUM) > 0
            @test all(isfinite, it.WIDTH_SPECTRUM)
            @test all(>(0), it.WIDTH_SPECTRUM)
        end
    end
end
