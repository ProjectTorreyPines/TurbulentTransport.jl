@testset "FINN model loading" begin
    @testset "load_finn_model returns FINNmodel or FINNensemble" begin
        model = TurbulentTransport.load_finn_model(TEST_FINN_MODEL)
        @test model isa Union{TurbulentTransport.FINNmodel, TurbulentTransport.FINNensemble}
    end

    @testset "FINNmodel fields are well-formed" begin
        model = TurbulentTransport.load_finn_model(TEST_FINN_MODEL)
        @test !isempty(model.name)
        @test length(model.xnames) == 17
        @test length(model.ynames) == 4
        @test length(model.xm) == length(model.xnames)
        @test length(model.xσ) == length(model.xnames)
        @test size(model.xbounds) == (length(model.xnames), 2)
        @test length(model.ym) == length(model.ynames)
        @test length(model.yσ) == length(model.ynames)
        @test size(model.ybounds) == (length(model.ynames), 2)
        @test all(model.xσ .> 0)
        @test all(model.yσ .> 0)
    end

    @testset "ensemble has multiple models" begin
        model = TurbulentTransport.load_finn_model(TEST_FINN_MODEL)
        if model isa TurbulentTransport.FINNensemble
            @test length(model.models) > 1
            for m in model.models
                @test m isa TurbulentTransport.FINNmodel
                @test m.xnames == model.models[1].xnames
                @test m.ynames == model.models[1].ynames
            end
        end
    end

    @testset "xnames include expected geometry and source inputs" begin
        model = TurbulentTransport.load_finn_model(TEST_FINN_MODEL)
        @test "Q_LOC"       ∈ model.xnames
        @test "KAPPA_LOC"   ∈ model.xnames
        @test "DRMAJDX_LOC" ∈ model.xnames
        @test "Qe"          ∈ model.xnames
        @test "Qi"          ∈ model.xnames
        @test "Ge"          ∈ model.xnames
        @test "Pi"          ∈ model.xnames
        @test "rho"         ∈ model.xnames
    end

    @testset "ynames include expected gradient outputs" begin
        model = TurbulentTransport.load_finn_model(TEST_FINN_MODEL)
        ynames_clean = [replace(yn, "OUT_" => "") for yn in model.ynames]
        @test "RLTS_1"     ∈ ynames_clean
        @test "RLTS_2"     ∈ ynames_clean
        @test "RLNS_1"     ∈ ynames_clean
        @test "VEXB_SHEAR" ∈ ynames_clean
    end

    @testset "load_finn_model_once is memoized (same object)" begin
        m1 = TurbulentTransport.load_finn_model_once(TEST_FINN_MODEL)
        m2 = TurbulentTransport.load_finn_model_once(TEST_FINN_MODEL)
        @test m1 === m2
    end

    @testset "show methods do not error" begin
        model = TurbulentTransport.load_finn_model(TEST_FINN_MODEL)
        @test !isempty(sprint(show, MIME"text/plain"(), model))
    end
end

@testset "FINN prediction" begin
    model = TurbulentTransport.load_finn_model(TEST_FINN_MODEL)

    x_mid = (model.xbounds[:, 1] .+ model.xbounds[:, 2]) ./ 2

    @testset "predict_finn vector input returns correct shape" begin
        y = TurbulentTransport.predict_finn(model, x_mid)
        @test y isa Vector{Float64}
        @test length(y) == length(model.ynames)
    end

    @testset "predict_finn matrix input returns correct shape" begin
        X = hcat(x_mid, x_mid .* 1.05, x_mid .* 0.95)
        Y = TurbulentTransport.predict_finn(model, X)
        @test Y isa Matrix{Float64}
        @test size(Y) == (length(model.ynames), 3)
    end

    @testset "predict_finn outputs are all finite" begin
        y = TurbulentTransport.predict_finn(model, x_mid)
        @test all(isfinite, y)

        X = hcat(x_mid, x_mid .* 1.05, x_mid .* 0.95)
        Y = TurbulentTransport.predict_finn(model, X)
        @test all(isfinite, Y)
    end

    @testset "vector and matrix predictions are consistent" begin
        y_vec = TurbulentTransport.predict_finn(model, x_mid)
        X = hcat(x_mid, x_mid .* 1.05)
        Y_mat = TurbulentTransport.predict_finn(model, X)
        @test Y_mat[:, 1] ≈ y_vec rtol=REGRESSION_RTOL
    end

    @testset "predictions are deterministic" begin
        y1 = TurbulentTransport.predict_finn(model, x_mid)
        y2 = TurbulentTransport.predict_finn(model, x_mid)
        @test y1 ≈ y2 rtol=REGRESSION_RTOL
    end

    @testset "different inputs give different outputs" begin
        x1 = x_mid
        x2 = x_mid .* 1.1
        y1 = TurbulentTransport.predict_finn(model, x1)
        y2 = TurbulentTransport.predict_finn(model, x2)
        @test !(y1 ≈ y2)
    end

    @testset "gradient outputs have physically reasonable signs and magnitudes" begin
        y = TurbulentTransport.predict_finn(model, x_mid)
        ynames_clean = [replace(yn, "OUT_" => "") for yn in model.ynames]
        result = Dict(name => y[k] for (k, name) in enumerate(ynames_clean))

        @test result["RLTS_1"] > 0
        @test result["RLTS_2"] > 0
        @test result["RLNS_1"] > 0
        @test abs(result["RLTS_1"]) < 50
        @test abs(result["RLTS_2"]) < 50
        @test abs(result["RLNS_1"]) < 50
        @test abs(result["VEXB_SHEAR"]) < 5
    end
end
