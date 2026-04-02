@testset "FINN model loading" begin
    @testset "load_finn_model returns FINNmodel" begin
        model = TurbulentTransport.load_finn_model(TEST_FINN_MODEL)
        @test model isa TurbulentTransport.FINNmodel
    end

    @testset "FINNmodel fields are well-formed" begin
        model = TurbulentTransport.load_finn_model(TEST_FINN_MODEL)
        @test !isempty(model.name)
        @test length(model.xnames) == 16
        @test length(model.ynames) >= 4
        @test length(model.xm) == length(model.xnames)
        @test length(model.xσ) == length(model.xnames)
        @test size(model.xbounds) == (length(model.xnames), 2)
        @test length(model.ym) == length(model.ynames)
        @test length(model.yσ) == length(model.ynames)
        @test size(model.ybounds) == (length(model.ynames), 2)
        @test all(model.xσ .> 0)
        @test all(model.yσ .> 0)
    end

    @testset "xnames include expected geometry and source inputs" begin
        model = TurbulentTransport.load_finn_model(TEST_FINN_MODEL)
        @test "Q_LOC"      ∈ model.xnames
        @test "KAPPA_LOC"  ∈ model.xnames
        @test "DRMINDX_LOC" ∈ model.xnames
        @test "DRMINDX_LOC" ∈ model.xnames
        @test "Qe"         ∈ model.xnames
        @test "Qi"         ∈ model.xnames
        @test "Ge"         ∈ model.xnames
        @test "Pi"         ∈ model.xnames
        @test "rho"        ∈ model.xnames
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

    # Midpoint of training bounds as a representative in-distribution input
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

    @testset "midpoint regression values" begin
        ynames_clean = [replace(yn, "OUT_" => "") for yn in model.ynames]
        y = TurbulentTransport.predict_finn(model, x_mid)
        result = Dict(name => y[k] for (k, name) in enumerate(ynames_clean))

        @test result["RLTS_1"]     ≈ EXPECTED_FINN_MIDPOINT.RLTS_1     rtol=REGRESSION_RTOL
        @test result["RLTS_2"]     ≈ EXPECTED_FINN_MIDPOINT.RLTS_2     rtol=REGRESSION_RTOL
        @test result["RLNS_1"]     ≈ EXPECTED_FINN_MIDPOINT.RLNS_1     rtol=REGRESSION_RTOL
        @test result["VEXB_SHEAR"] ≈ EXPECTED_FINN_MIDPOINT.VEXB_SHEAR rtol=REGRESSION_RTOL
    end

    @testset "matrix regression values (col 1 == vector)" begin
        ynames_clean = [replace(yn, "OUT_" => "") for yn in model.ynames]
        X = hcat(x_mid, x_mid .* 1.05, x_mid .* 0.95)
        Y = TurbulentTransport.predict_finn(model, X)
        result = Dict(name => Y[k, 1] for (k, name) in enumerate(ynames_clean))

        @test result["RLTS_1"]     ≈ EXPECTED_FINN_MATRIX_COL1.RLTS_1     rtol=REGRESSION_RTOL
        @test result["RLTS_2"]     ≈ EXPECTED_FINN_MATRIX_COL1.RLTS_2     rtol=REGRESSION_RTOL
        @test result["RLNS_1"]     ≈ EXPECTED_FINN_MATRIX_COL1.RLNS_1     rtol=REGRESSION_RTOL
        @test result["VEXB_SHEAR"] ≈ EXPECTED_FINN_MATRIX_COL1.VEXB_SHEAR rtol=REGRESSION_RTOL
    end

    @testset "different inputs give different outputs" begin
        x1 = x_mid
        x2 = x_mid .* 1.1
        y1 = TurbulentTransport.predict_finn(model, x1)
        y2 = TurbulentTransport.predict_finn(model, x2)
        @test !(y1 ≈ y2)
    end
end
