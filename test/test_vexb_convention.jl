# VEXB_SHEAR sign convention: every NN receives the sign it was trained on.
#
# Internal (canonical) convention is TGYRO: VEXB_SHEAR = -SIGN_BT*gamma_e*a/c_s, same sign
# as VPAR_SHEAR_1 (what `InputTGLF(dd, ...)` writes and what plain input.tglf files carry).
# Models tagged `:legacy` were trained on SIGN_BT * that value and are fed the flipped feature.

using TurbulentTransport: _default_vexb_convention, _vexb_row, vexb_sign

# Feature vector of `it` in the TGYRO convention, as the model's xnames order (log10 applied by flux_array)
function _features(model, it)
    return Float64[getfield(it, Symbol(replace(n, "_log10" => ""))) for n in model.xnames]
end
_solution(model, x) = TurbulentTransport.flux_solution(TurbulentTransport.flux_array(model, x; warn_nn_train_bounds=false))

_fluxes(s) = (s.ENERGY_FLUX_e, s.ENERGY_FLUX_i, s.PARTICLE_FLUX_e, s.STRESS_TOR_i)
# vector and matrix inference paths differ at the last ulp
_same(a, b) = all(isapprox.(_fluxes(a), _fluxes(b); rtol=REGRESSION_RTOL))

@testset "VEXB_SHEAR sign convention" begin
    it = load_sample_input()
    @test it.SIGN_BT == -1           # the sample input exercises the flip
    @test vexb_sign(it) == -1

    @testset "convention table" begin
        @test _default_vexb_convention("sat3_em_d3d_azf-1") === :legacy
        @test _default_vexb_convention("sat3_em_iter4_azf-1") === :legacy
        @test _default_vexb_convention("sat3_em_d3d_azf-1_gknne24") === :tgyro
        @test _default_vexb_convention("sat3_em_d3d_azf-1_withnegD_gknn37") === :tgyro
        @test _default_vexb_convention("sat3_em_mastuedge+nstxedge_azf-1_withnegD") === :tgyro
        @test _default_vexb_convention("sat3_em_mastu+nstx_azf-1_withnegD") === :legacy
        @test _default_vexb_convention("modeid_qlgyro_sat3_azf-1") === :tgyro
        @test _default_vexb_convention("some_future_model") === :legacy
        # every listed name is a shipped model
        for name in TurbulentTransport._VEXB_TGYRO_MODELS
            @test isfile(joinpath(pkgdir(TurbulentTransport), "models", name * ".bson"))
        end
        @test_throws ErrorException TurbulentTransport._validate_vexb_convention(:bogus, "test")
    end

    @testset "tags survive loading (bson without the key -> table default)" begin
        @test loadmodel("sat3_em_d3d_azf-1").vexb_convention === :legacy
        @test loadmodel("sat3_em_d3d_azf-1_gknne24").vexb_convention === :tgyro
        @test loadmodel("sat3_em_mastuedge+nstxedge_azf-1_withnegD").vexb_convention === :tgyro
        @test TurbulentTransport.load_modeid_model(TEST_MODEID_MODEL).vexb_convention === :tgyro
        # explicit key wins over the table
        d = Dict{Any,Any}(k => v for (k, v) in TurbulentTransport.mod2dict(loadmodel("sat3_em_d3d_azf-1").models[1]))
        d[:vexb_convention] = :tgyro
        @test TurbulentTransport.dict2mod(d; model_basename="sat3_em_d3d_azf-1").vexb_convention === :tgyro
        d[:vexb_convention] = "legacy"     # trainers may write a String
        @test TurbulentTransport.dict2mod(d; model_basename="sat3_em_d3d_azf-1_gknne24").vexb_convention === :legacy
        d[:vexb_convention] = :bogus
        @test_throws ErrorException TurbulentTransport.dict2mod(d; model_basename="x")
    end

    @testset "legacy model sees SIGN_BT * VEXB_SHEAR (bit-identical to the pre-1.4 feed)" begin
        model = loadmodel("sat3_em_d3d_azf-1")
        sol = TurbulentTransport.run_tglfnn(it; model_filename="sat3_em_d3d_azf-1", warn_nn_train_bounds=false)
        x = _features(model, it)
        x[_vexb_row(model.xnames)] *= it.SIGN_BT
        @test _same(sol, _solution(model, x))
        # SIGN_BT=+1 with the same VEXB_SHEAR reproduces what the model was fed before 1.4
        it_plus = deepcopy(it); it_plus.SIGN_BT = 1
        sol_plus = TurbulentTransport.run_tglfnn(it_plus; model_filename="sat3_em_d3d_azf-1", warn_nn_train_bounds=false)
        for (a, b) in zip(_fluxes(sol_plus), EXPECTED_RUN_TGLFNN_SAT3_SIGNBT_PLUS)
            @test isapprox(a, b; rtol=REGRESSION_RTOL)
        end
        # (SIGN_BT=-1, v) and (SIGN_BT=+1, -v) are the same legacy feed
        it_neg = deepcopy(it); it_neg.SIGN_BT = 1; it_neg.VEXB_SHEAR = -it.VEXB_SHEAR
        sol_neg = TurbulentTransport.run_tglfnn(it_neg; model_filename="sat3_em_d3d_azf-1", warn_nn_train_bounds=false)
        @test _fluxes(sol) == _fluxes(sol_neg)
        # Dict path with SIGN_BT present agrees with the InputTGLF path
        data = Dict{String,Any}(n => [getfield(it, Symbol(n))] for n in replace.(model.xnames, "_log10" => ""))
        data["SIGN_BT"] = [it.SIGN_BT]
        y = TurbulentTransport.run_tglfnn(data; model_filename="sat3_em_d3d_azf-1", warn_nn_train_bounds=false)
        @test y["Q_elec"][1] == sol.ENERGY_FLUX_e
        # ... and without SIGN_BT the dictionary is taken as already in the model's convention
        delete!(data, "SIGN_BT")
        y2 = TurbulentTransport.run_tglfnn(data; model_filename="sat3_em_d3d_azf-1", warn_nn_train_bounds=false)
        @test y2["Q_elec"][1] == sol_plus.ENERGY_FLUX_e
    end

    @testset "tgyro nets are not flipped: GKNN correction distinguishes the two feeds" begin
        it_neg = deepcopy(it); it_neg.SIGN_BT = 1; it_neg.VEXB_SHEAR = -it.VEXB_SHEAR
        # the legacy base net sees identical inputs (checked above) ...
        a = TurbulentTransport.run_tglfnn(it; model_filename="sat3_em_d3d_azf-1", warn_nn_train_bounds=false)
        b = TurbulentTransport.run_tglfnn(it_neg; model_filename="sat3_em_d3d_azf-1", warn_nn_train_bounds=false)
        @test _fluxes(a) == _fluxes(b)
        # ... while the tgyro GKNN nets see v vs -v, so the corrected fluxes must differ
        ga = TurbulentTransport.run_tglfnn(it; model_filename="sat3_em_d3d_azf-1", warn_nn_train_bounds=false, fidelity=:GKNN)
        gb = TurbulentTransport.run_tglfnn(it_neg; model_filename="sat3_em_d3d_azf-1", warn_nn_train_bounds=false, fidelity=:GKNN)
        @test _fluxes(ga) != _fluxes(gb)
        # direct check of the composition for the electron heat channel
        base = loadmodel("sat3_em_d3d_azf-1_tglfnn24")
        gknne = loadmodel("sat3_em_d3d_azf-1_gknne24")
        xb = _features(base, it); xb[_vexb_row(base.xnames)] *= it.SIGN_BT
        yb = TurbulentTransport.flux_array(base, xb; warn_nn_train_bounds=false)
        xg = vcat(_features(base, it), yb[3])           # tgyro feed, unflipped, + base Q_elec
        err = TurbulentTransport.flux_array(gknne, xg; warn_nn_train_bounds=false, fidelity=:GKNN)[1]
        @test ga.ENERGY_FLUX_e ≈ yb[3] * err rtol=REGRESSION_RTOL
    end

    @testset "mixed conventions inside one call" begin
        # ST family: core net is :legacy, its radial variants are :tgyro
        core_name = "sat3_em_mastu+nstx_azf-1_withnegD"
        nearedge_name, edge_name, r_ne, r_e = TurbulentTransport.radial_blend_variants(core_name)
        core = loadmodel(core_name); edge = loadmodel(edge_name)
        @test core.vexb_convention === :legacy
        @test edge.vexb_convention === :tgyro
        it_core = deepcopy(it); it_core.RMIN_LOC = 0.5
        it_edge = deepcopy(it); it_edge.RMIN_LOC = 0.99
        sols = TurbulentTransport.run_tglfnn([it_core, it_edge]; model_filename=core_name, warn_nn_train_bounds=false)
        xc = _features(core, it_core); xc[_vexb_row(core.xnames)] *= it.SIGN_BT
        @test _same(sols[1], _solution(core, xc))
        @test _same(sols[2], _solution(edge, _features(edge, it_edge)))   # no flip

        # mixed SIGN_BT columns in one GKNN batch == singles
        it_plus = deepcopy(it); it_plus.SIGN_BT = 1
        batch = TurbulentTransport.run_tglfnn([it, it_plus, it]; model_filename="sat3_em_d3d_azf-1_withnegD", warn_nn_train_bounds=false, fidelity=:GKNN)
        singles = [TurbulentTransport.run_tglfnn(x; model_filename="sat3_em_d3d_azf-1_withnegD", warn_nn_train_bounds=false, fidelity=:GKNN) for x in (it, it_plus, it)]
        for (b, s) in zip(batch, singles)
            @test all(isapprox.(_fluxes(b), _fluxes(s); rtol=REGRESSION_RTOL))
        end
        @test _fluxes(batch[1]) != _fluxes(batch[2])
    end

    @testset "InputTGLF(dd) writes the TGYRO sign: VEXB_SHEAR == VPAR_SHEAR_1 * r/(|q| R)" begin
        dd = load_sample_dd()
        its = InputTGLF(dd, [0.3, 0.5, 0.7], :sat3, true, true)
        for k in eachindex(its.tglfs)
            t = its.tglfs[k]
            expected = t.VPAR_SHEAR_1 * t.RMIN_LOC / (t.Q_LOC * t.RMAJ_LOC)
            @test isapprox(t.VEXB_SHEAR, expected; rtol=1e-10, atol=1e-14)
        end
    end

    @testset "QLNN bundles" begin
        qlnn_dir = joinpath(pkgdir(TurbulentTransport), "models")
        if isdir(joinpath(qlnn_dir, "QLNN_d3d_1")) && isdir(joinpath(qlnn_dir, "QLNN"))
            @test TurbulentTransport._qlnn_read_vexb_convention(joinpath(qlnn_dir, "QLNN_d3d_1")) === :legacy
            @test TurbulentTransport._qlnn_read_vexb_convention(joinpath(qlnn_dir, "QLNN")) === :tgyro
            @test TurbulentTransport.loadqlnnbundle("QLNN_d3d_1").vexb_convention === :legacy
            bundle = TurbulentTransport.loadqlnnbundle("QLNN")
            @test bundle.vexb_convention === :tgyro
            itj = load_sample_input_tjlf()
            @test itj.SIGN_BT == -1
            xnames = bundle.energy.xnames
            ks = [0.1, 0.3, 1.0]
            xs_t = TurbulentTransport._qlnn_build_xs(itj, ks, xnames; vexb_convention=:tgyro)
            xs_l = TurbulentTransport._qlnn_build_xs(itj, ks, xnames; vexb_convention=:legacy)
            k = findfirst(==("VEXB_SHEAR"), xnames)
            if k !== nothing
                @test xs_l[k, :] == -xs_t[k, :]
                rows = setdiff(1:length(xnames), k)
                @test xs_l[rows, :] == xs_t[rows, :]
            end
        else
            @info "Skipping QLNN vexb_convention tests; bundles not found"
        end
    end

    @testset "table agrees with the stored input means" begin
        # For SIGN_BT=-1 training devices, sign(xm[VEXB_SHEAR]) == sign(xm[VPAR_SHEAR_1]) iff the
        # corpus was in the TGYRO convention. Models with negligible rotation are skipped.
        names = if get(ENV, "TT_TEST_ALL_MODELS", "false") == "true"
            filter(n -> isfile(joinpath(pkgdir(TurbulentTransport), "models", n * ".bson")), available_models())
        else
            ["sat3_em_d3d_azf-1", "sat3_em_iter4_azf-1", "sat3_em_nstx_azf-1", "sat3_em_mastu+nstx_azf-1_withnegD",
             "sat3_em_d3dedge_azf-1_withnegD", "sat3_em_mastuedge+nstxedge_azf-1_withnegD",
             "sat3_em_mastunearedge+nstxnearedge_azf-1_withnegD", "sat3_em_d3d_azf-1_gknne24",
             "sat3_em_d3d_azf-1_withnegD_gknn37", "sat3_em_d3d+mastu+nstx_azf-1_gknn31", "sat3_em_d3d_azf-1_qlnn",
             "sat3_em_d3d_azf-1_gkdb_gknn31_cgyro"]
        end
        for name in names
            startswith(name, "finn") && continue
            m = loadmodel(name)
            @test m.vexb_convention in TurbulentTransport._VEXB_CONVENTIONS
            kv = findfirst(==("VEXB_SHEAR"), m.xnames); kp = findfirst(==("VPAR_SHEAR_1"), m.xnames)
            (kv === nothing || kp === nothing) && continue
            (abs(m.xm[kv]) > 0.01 && abs(m.xm[kp]) > 0.05) || continue
            same = sign(m.xm[kv]) == sign(m.xm[kp])
            @test (m.vexb_convention === :tgyro) == same
        end
    end
end
