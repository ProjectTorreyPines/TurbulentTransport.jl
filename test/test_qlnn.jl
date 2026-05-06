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
            input_tglf = load_sample_input()
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
            input_tglf = load_sample_input()
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
            # SAT_RULE=2 honors ALPHA_ZF in `intensity_sat`; flipping the sign
            # must change the integrated electron heat flux.
            qe_zf_neg = _qlnn_qe(; sat_rule=2, alpha_zf=-1.0)
            qe_zf_pos = _qlnn_qe(; sat_rule=2, alpha_zf=+1.0)
            @test isfinite(qe_zf_neg)
            @test isfinite(qe_zf_pos)
            @test !isapprox(qe_zf_neg, qe_zf_pos; rtol=1e-4)
        end

        @testset "warn_nn_train_bounds emits @warn for out-of-range inputs" begin
            # Build a baseline `InputTJLF` that's well inside the training
            # distribution, then push a single feature (RMIN_LOC) far past
            # `bundle.energy.xbounds[i, 2]` and confirm `run_qlnn` issues a
            # warning when `warn_nn_train_bounds=true` and stays silent when
            # `warn_nn_train_bounds=false`. We can't predict exactly which
            # feature index will be RMIN_LOC, so we set every feature in turn
            # to the upper bound + 10·xσ — guaranteed to trigger.
            input_tglf = load_sample_input()
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
