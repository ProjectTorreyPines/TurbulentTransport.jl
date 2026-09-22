# Smoke tests for the :QLNN flux matcher.
#
# What we cover:
#   1. The bundle directory `models/QLNN/` is discoverable via
#      `available_qlnn_bundles()` and `loadqlnnbundle("QLNN")`.
#   2. `run_qlnn` produces a `GACODE.FluxSolution` with finite, real-valued
#      fluxes from `sample_input.tglf`.
#   3. `SAT_RULE` and `ALPHA_ZF` actually flow through `TJLF.sum_ky_spectrum`
#      — toggling either one changes the integrated flux.
#
# These tests intentionally check *behavior* (sat-rule sensitivity), not
# numerical regression. They run only if `models/QLNN/` is present so they
# don't fail in CI environments that don't ship the QLNN bundle.

import TJLF
import ForwardDiff
import Logging

const QLNN_BUNDLE_NAME = "QLNN"
const QLNN_BUNDLE_DIR = joinpath(dirname(@__DIR__), "models", QLNN_BUNDLE_NAME)

# Skip the whole testset if the bundle isn't present (keeps CI green when
# running on a worker without the QLNN models checked in).
if !isdir(QLNN_BUNDLE_DIR)
    @info "Skipping QLNN smoke tests; bundle directory not found: $QLNN_BUNDLE_DIR"
else
    @testset "QLNN smoke tests" begin
        @testset "Discovery + bundle load" begin
            bundles = TurbulentTransport.available_qlnn_bundles()
            @test QLNN_BUNDLE_NAME in bundles

            bundle = TurbulentTransport.loadqlnnbundle(QLNN_BUNDLE_NAME)
            @test bundle isa TurbulentTransport.QLNNbundle
            # Note: bundles can be either `QLNNmodel` or `QLNNensemble`.
            # Use property access (not `getfield`) so the ensemble's
            # `Base.getproperty` forwarding to its first member is exercised.
            for fname in (:energy, :particle, :momentum, :eigenvalue)
                m = getfield(bundle, fname)
                @test m isa TurbulentTransport.AbstractQLNNmodel
                @test m.target === fname
                @test !isempty(m.xnames)
                @test !isempty(m.ynames)
            end
            # The QLNN bundle that ships in this repo contains the
            # stability classifier; if a future bundle drops it, the loader
            # gracefully sets `bundle.stability = nothing`.
            if isfile(joinpath(QLNN_BUNDLE_DIR, "stability_classifier.bson"))
                @test bundle.stability !== nothing
                @test bundle.stability.target === :stability
            else
                @test bundle.stability === nothing
            end
        end

        @testset "run_qlnn produces a finite FluxSolution" begin
            input_tglf = load_sample_input_cgyro()
            input_tjlf = InputTJLF{Float64}(input_tglf)

            sol = TurbulentTransport.run_qlnn(input_tjlf;
                                              bundle_name=QLNN_BUNDLE_NAME,
                                              warn_nn_train_bounds=false)

            @test sol isa TurbulentTransport.GACODE.FluxSolution
            @test isfinite(sol.ENERGY_FLUX_e)
            @test isfinite(sol.ENERGY_FLUX_i)
            @test isfinite(sol.PARTICLE_FLUX_e)
            @test isfinite(sol.STRESS_TOR_i)
            for v in sol.PARTICLE_FLUX_i
                @test isfinite(v)
            end
        end

        # Helper: integrated electron energy flux for a given SAT_RULE / ALPHA_ZF
        # override. Each call gets a fresh InputTJLF (so width memory and
        # KY_SPECTRUM resets don't bleed across runs).
        function _qlnn_qe(; sat_rule::Int, alpha_zf::Real)
            input_tglf = load_sample_input_cgyro()
            input_tjlf = InputTJLF{Float64}(input_tglf)
            input_tjlf.SAT_RULE = sat_rule
            input_tjlf.ALPHA_ZF = Float64(alpha_zf)
            sol = TurbulentTransport.run_qlnn(input_tjlf;
                                              bundle_name=QLNN_BUNDLE_NAME,
                                              warn_nn_train_bounds=false)
            return sol.ENERGY_FLUX_e
        end

        @testset "SAT_RULE flows into the integrated flux" begin
            # Pick two saturation rules that share the same QL packing
            # (sat1 vs sat2) so the only thing that changes is the rule.
            qe_sat1 = _qlnn_qe(; sat_rule=1, alpha_zf=-1.0)
            qe_sat2 = _qlnn_qe(; sat_rule=2, alpha_zf=-1.0)
            @test isfinite(qe_sat1)
            @test isfinite(qe_sat2)
            @test !isapprox(qe_sat1, qe_sat2; rtol=1e-4)
        end

        @testset "ALPHA_ZF flows into the integrated flux" begin
            # SAT_RULE=2 honors ALPHA_ZF via `czf = abs(alpha_zf)` in
            # `intensity_sat`, so changing the magnitude must change the
            # integrated electron heat flux. (Don't test the sign flip: the
            # sign only gates a low-k kymin cutoff in the zonal-mixing peak
            # search, which is a no-op whenever the NN-predicted gamma/ky
            # spectrum peaks above the cutoff — model- and input-dependent.)
            qe_zf_1 = _qlnn_qe(; sat_rule=2, alpha_zf=-1.0)
            qe_zf_h = _qlnn_qe(; sat_rule=2, alpha_zf=-0.5)
            @test isfinite(qe_zf_1)
            @test isfinite(qe_zf_h)
            @test !isapprox(qe_zf_1, qe_zf_h; rtol=1e-4)
        end

        @testset "ALPHA_ZF sign flows into the integrated flux" begin
            # The sign branch (kymin cutoff) is only active when the predicted
            # gamma/ky spectrum peaks below the cutoff, so it needs an input
            # known to exercise it for the current default bundle:
            # input_zf_sign.tglf (TJLF regression case tglf13).
            function _qe_sign(alpha_zf::Real)
                it = TJLF.readInput(joinpath(@__DIR__, "data", "input_zf_sign.tglf"))
                it.SAT_RULE = 2
                it.UNITS = "CGYRO"
                it.ALPHA_ZF = Float64(alpha_zf)
                sol = TurbulentTransport.run_qlnn(it;
                                                  bundle_name=QLNN_BUNDLE_NAME,
                                                  warn_nn_train_bounds=false)
                return sol.ENERGY_FLUX_e
            end
            qe_neg = _qe_sign(-1.0)
            qe_pos = _qe_sign(+1.0)
            @test isfinite(qe_neg)
            @test isfinite(qe_pos)
            @test !isapprox(qe_neg, qe_pos; rtol=1e-4)
        end

        @testset "warn_nn_train_bounds emits @warn for out-of-range inputs" begin
            # Build a baseline `InputTJLF` that's well inside the training
            # distribution, then push a single feature (RMIN_LOC) far past
            # `bundle.energy.xbounds[i, 2]` and confirm `run_qlnn` issues a
            # warning when `warn_nn_train_bounds=true` and stays silent when
            # `warn_nn_train_bounds=false`. We can't predict exactly which
            # feature index will be RMIN_LOC, so we set every feature in turn
            # to the upper bound + 10·xσ — guaranteed to trigger.
            input_tglf = load_sample_input_cgyro()
            bundle = TurbulentTransport.loadqlnnbundle(QLNN_BUNDLE_NAME)
            xnames = bundle.energy.xnames
            xbounds = bundle.energy.xbounds
            # Skip the test if the loaded bundle has ±Inf bounds (legacy BSON
            # without `:xbounds`) — there'd be nothing to warn about.
            has_finite_bounds = any(isfinite, xbounds)
            if !has_finite_bounds
                @info "Skipping warn_nn_train_bounds test: bundle has no :xbounds field."
            else
                # Find a plain scalar `InputTJLF` field with a finite upper
                # bound so we can deterministically push it past the training
                # range. Skip `ky` (per-column, not a struct field) and skip
                # species-suffixed names (e.g. `AS_2`, `VPAR_SHEAR_3`) which
                # require indexing into a vector field.
                species_rx = r"^(.+)_(\d+)$"
                plain_idx = findfirst(eachindex(xnames)) do i
                    nm = xnames[i]
                    nm == "ky" && return false
                    isfinite(xbounds[i, 2]) || return false
                    base = endswith(nm, "_log10") ? nm[1:end-6] : nm
                    match(species_rx, base) === nothing &&
                        hasfield(typeof(InputTJLF{Float64}(input_tglf)), Symbol(base))
                end
                if plain_idx === nothing
                    @info "Skipping warn_nn_train_bounds value-injection: no plain scalar feature with finite bounds."
                else
                    nm = xnames[plain_idx]
                    base = endswith(nm, "_log10") ? nm[1:end-6] : nm
                    fname = Symbol(base)
                    # Set the value well above the upper bound. For `_log10`
                    # features we want post-log10(value) > xbounds[i, 2], so we
                    # raise the linear value to 10^(xbounds[i, 2] + 2).
                    over_value = if endswith(nm, "_log10")
                        10.0 ^ (Float64(xbounds[plain_idx, 2]) + 2.0)
                    else
                        Float64(xbounds[plain_idx, 2]) + 1.0e6
                    end
                    input_tjlf = InputTJLF{Float64}(input_tglf)
                    setfield!(input_tjlf, fname, over_value)
                    @test_logs (:warn,) match_mode=:any TurbulentTransport.run_qlnn(
                        input_tjlf; bundle_name=QLNN_BUNDLE_NAME, warn_nn_train_bounds=true)
                    # And: no warning when the flag is off.
                    input_tjlf2 = InputTJLF{Float64}(input_tglf)
                    setfield!(input_tjlf2, fname, over_value)
                    @test_logs min_level=Logging.Warn TurbulentTransport.run_qlnn(
                        input_tjlf2; bundle_name=QLNN_BUNDLE_NAME, warn_nn_train_bounds=false)
                end
            end
        end

        @testset "ForwardDiff: predict preserves Dual eltype" begin
            # Narrow AD smoke test for the new code in qlnn.jl: feed a
            # `Matrix{Dual}` straight into `predict` (regressor + classifier)
            # and check the eltype + partials survive the Flux.Chain forward
            # pass. End-to-end Dual propagation through `run_qlnn` →
            # `TJLF.sum_ky_spectrum` is exercised by FUSE's existing `:TJLF`
            # `:forward_ad` integration test, which we extended to cover
            # `:QLNN` in the same code path (`ad_flux_match_errors!`).
            #
            # If this test passes, the regressors and classifier are
            # Dual-compatible; the rest of the pipeline (QL packing + sum_ky)
            # already supports Duals via the `T<:Real` parameterization that
            # the `:TJLF` AD path relies on.
            bundle = TurbulentTransport.loadqlnnbundle(QLNN_BUNDLE_NAME)
            nf = length(bundle.energy.xnames)
            nky = 5

            # Build a Float64 input matrix that's strictly positive so
            # `_qlnn_apply_log10!` can run on `_log10`-suffixed feature rows
            # without hitting `log10(negative) → DomainError`. In production
            # the log10-tagged features (BETAE, DEBYE, XNUE, ...) are always
            # positive; `1.0 + 0.01·N(0,1)` keeps every entry well above
            # zero (1 - 4σ = 0.96 > 0) and stays in a numerically tame range.
            xs_f64 = 1.0 .+ 0.01 .* randn(nf, nky)

            # Promote to Matrix{Dual{Tag,Float64,1}} with the partial seeded
            # on the first feature row. After a successful forward pass the
            # output should carry non-zero partials in the same row.
            Tag = ForwardDiff.Tag{:qlnn_predict_test, Float64}
            DualT = ForwardDiff.Dual{Tag, Float64, 1}
            xs_dual = Matrix{DualT}(undef, nf, nky)
            for i in 1:nf, j in 1:nky
                partial = (i == 1 ? 1.0 : 0.0)
                xs_dual[i, j] = ForwardDiff.Dual{Tag}(xs_f64[i, j], partial)
            end

            # Regressor predict: output must be Dual-typed with finite values
            # and finite partials.
            y_energy = TurbulentTransport.predict(bundle.energy, xs_dual)
            @test eltype(y_energy) <: ForwardDiff.Dual
            @test all(isfinite, ForwardDiff.value.(y_energy))
            @test all(isfinite, ForwardDiff.partials.(y_energy, 1))
            # Sanity: at least one partial is non-zero — otherwise the chain
            # silently down-converted to Float64 and stripped derivatives.
            @test any(p -> p != 0.0, ForwardDiff.partials.(y_energy, 1))

            # Classifier predict: same contract for the stability head (when
            # the bundle ships one). `predict_unstable_prob` runs σ on top of
            # `predict`, so it also has to preserve the Dual eltype.
            if bundle.stability !== nothing
                p_un = TurbulentTransport.predict_unstable_prob(bundle.stability, xs_dual)
                @test eltype(p_un) <: ForwardDiff.Dual
                @test all(0.0 .<= ForwardDiff.value.(p_un) .<= 1.0)
                @test all(isfinite, ForwardDiff.partials.(p_un, 1))
            end
        end
    end
end

# DT-lumped bundle (QLNN_ukstep26): sidecars, automatic D+T lumping, and the NS=3 guard.
# Skipped when the bundle directory is absent (e.g. a checkout without LFS content).
const QLNN_DT_BUNDLE = "QLNN_ukstep26"

# (e, D, T, C) NS=4 variant of the sample input: split the D density 70/30 into D and T,
# move C to slot 4. Mirrors the corpus D-T split (`_apply_stfpp_transform!`).
function _sample_input_unbundled(; he_ash::Bool=false)
    g = load_sample_input_cgyro()
    for p in (:AS, :ZS, :MASS, :RLNS, :RLTS, :TAUS, :VPAR, :VPAR_SHEAR)
        setproperty!(g, Symbol(p, "_4"), getproperty(g, Symbol(p, "_3")))
    end
    as2 = g.AS_2
    g.AS_2 = 0.7 * as2; g.AS_3 = 0.3 * as2
    g.ZS_3 = 1.0; g.MASS_3 = 1.49760170089
    g.RLNS_3 = g.RLNS_2 + 0.2; g.RLTS_3 = g.RLTS_2 - 0.1; g.TAUS_3 = g.TAUS_2 * 1.1
    g.NS = 4
    if he_ash
        for p in (:AS, :ZS, :MASS, :RLNS, :RLTS, :TAUS, :VPAR, :VPAR_SHEAR)
            setproperty!(g, Symbol(p, "_5"), getproperty(g, Symbol(p, "_4")))
        end
        g.ZS_5 = 2.0; g.MASS_5 = 2.0; g.AS_5 = 0.02; g.RLNS_5 = 0.1
        g.NS = 5
    end
    return g
end

if isdir(joinpath(dirname(@__DIR__), "models", QLNN_DT_BUNDLE))
    @testset "DT-lumped QLNN bundle ($QLNN_DT_BUNDLE)" begin
        bundle = TurbulentTransport.loadqlnnbundle(QLNN_DT_BUNDLE)
        @test bundle.vexb_convention === :legacy
        @test bundle.momentum_sign == 1.0
        @test bundle.stability !== nothing
        info = TurbulentTransport._qlnn_parse_qlweight_ynames(bundle.energy.ynames)
        @test info.species_set == ["e", "DT", "imp"]
        @test info.ns == 3
        @test TurbulentTransport._qlnn_is_dt_lumped(bundle)
        @test !TurbulentTransport._qlnn_is_dt_lumped(TurbulentTransport.loadqlnnbundle(QLNN_BUNDLE_NAME))
        @test isfile(joinpath(bundle.dir, "input.tglf.template"))

        # NS=3 input (sample input is e, D, C): used as is, finite fluxes
        it3 = InputTJLF{Float64}(load_sample_input_cgyro())
        @test it3.NS == 3
        @test TurbulentTransport.qlnn_lump_dt(it3) === it3
        sol3 = TurbulentTransport.run_qlnn(it3; bundle_name=QLNN_DT_BUNDLE, warn_nn_train_bounds=false)
        @test isfinite(sol3.ENERGY_FLUX_e) && isfinite(sol3.ENERGY_FLUX_i)

        @testset "qlnn_lump_dt NS=4" begin
            g = _sample_input_unbundled()
            it4 = InputTJLF{Float64}(g)
            @test it4.NS == 4
            l = TurbulentTransport.qlnn_lump_dt(it4)
            @test l !== it4 && it4.NS == 4          # copy; caller untouched
            @test l.NS == 3
            @test length(l.AS) == 3 && length(l.MASS) == 3
            n = g.AS_2 + g.AS_3; wd = g.AS_2 / n; wt = g.AS_3 / n
            @test l.AS[2] ≈ n
            @test l.ZS[2] == 1
            @test l.MASS[2] ≈ wd * 1.0 + wt * 1.49760170089
            @test l.TAUS[2] ≈ wd * g.TAUS_2 + wt * g.TAUS_3
            @test l.RLNS[2] ≈ wd * g.RLNS_2 + wt * g.RLNS_3
            @test l.RLTS[2] ≈ wd * g.RLTS_2 + wt * g.RLTS_3
            # impurity: old slot 4 (= the sample's C) moves to slot 3
            @test l.ZS[3] == g.ZS_4 && l.MASS[3] == g.MASS_4 && l.AS[3] == g.AS_4
            @test l.RLNS[3] == g.RLNS_4 && l.TAUS[3] == g.TAUS_4
            # electrons and scalars untouched
            @test l.AS[1] == it4.AS[1] && l.RLTS[1] == it4.RLTS[1]
            @test l.BETAE == it4.BETAE && l.Q_LOC == it4.Q_LOC && isequal(l.KY_SPECTRUM, it4.KY_SPECTRUM)
        end

        @testset "qlnn_lump_dt NS=5 (He ash into the impurity)" begin
            g = _sample_input_unbundled(he_ash=true)
            it5 = InputTJLF{Float64}(g)
            l = TurbulentTransport.qlnn_lump_dt(it5)
            @test l.NS == 3
            q  = g.AS_4 * g.ZS_4 + g.AS_5 * g.ZS_5
            z2 = g.AS_4 * g.ZS_4^2 + g.AS_5 * g.ZS_5^2
            @test l.AS[3] * l.ZS[3] ≈ q            # charge density conserved
            @test l.AS[3] * l.ZS[3]^2 ≈ z2         # Zeff contribution conserved
            @test l.MASS[3] ≈ (g.AS_4 * g.MASS_4 + g.AS_5 * g.MASS_5) / (g.AS_4 + g.AS_5)
            @test l.RLNS[3] ≈ (g.AS_4 * g.ZS_4 * g.RLNS_4 + g.AS_5 * g.ZS_5 * g.RLNS_5) / q
        end

        @testset "run_qlnn lumps automatically" begin
            it4 = InputTJLF{Float64}(_sample_input_unbundled())
            sol4 = TurbulentTransport.run_qlnn(it4; bundle_name=QLNN_DT_BUNDLE, warn_nn_train_bounds=false)
            @test it4.NS == 4                        # caller's input not modified
            itl = TurbulentTransport.qlnn_lump_dt(InputTJLF{Float64}(_sample_input_unbundled()))
            soll = TurbulentTransport.run_qlnn(itl; bundle_name=QLNN_DT_BUNDLE, warn_nn_train_bounds=false)
            @test sol4.ENERGY_FLUX_e == soll.ENERGY_FLUX_e
            @test sol4.ENERGY_FLUX_i == soll.ENERGY_FLUX_i
            @test sol4.PARTICLE_FLUX_e == soll.PARTICLE_FLUX_e
            @test sol4.STRESS_TOR_i == soll.STRESS_TOR_i
            # spectra path takes the same route
            sp = TurbulentTransport.qlnn_fluctuation_spectra(InputTJLF{Float64}(_sample_input_unbundled());
                                                             bundle_name=QLNN_DT_BUNDLE)
            @test sp !== nothing
        end

        @testset "non-hydrogenic slot 3 is refused" begin
            g = _sample_input_unbundled()
            g.ZS_3 = 6.0; g.MASS_3 = 6.0             # (e, D, C, C): not a D-T pair
            it = InputTJLF{Float64}(g)
            @test_throws ErrorException TurbulentTransport.qlnn_lump_dt(it)
            @test_throws ErrorException TurbulentTransport.run_qlnn(it; bundle_name=QLNN_DT_BUNDLE, warn_nn_train_bounds=false)
            # direct predictor call with NS=4 hits the safety net
            @test_throws ErrorException TurbulentTransport._run_qlnn_predict([InputTJLF{Float64}(_sample_input_unbundled())], bundle)
        end
    end
end
