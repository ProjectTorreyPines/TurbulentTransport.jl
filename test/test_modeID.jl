@testset "ModeID model loading" begin
    @testset "load_modeid_model returns ModeIDmodel or ModeIDensemble" begin
        model = TurbulentTransport.load_modeid_model(TEST_MODEID_MODEL)
        @test model isa Union{TurbulentTransport.ModeIDmodel, TurbulentTransport.ModeIDensemble}
    end

    @testset "ModeIDmodel fields are well-formed" begin
        model = TurbulentTransport.load_modeid_model(TEST_MODEID_MODEL)
        @test !isempty(model.name)
        @test length(model.xnames) == MODEID_N_INPUTS
        @test length(model.ynames) == MODEID_N_CLASSES
        @test length(model.xm) == MODEID_N_INPUTS
        @test length(model.xσ) == MODEID_N_INPUTS
        @test size(model.xbounds) == (MODEID_N_INPUTS, 2)
        @test all(model.xσ .> 0)
    end

    @testset "ynames are exactly the expected mode labels" begin
        model = TurbulentTransport.load_modeid_model(TEST_MODEID_MODEL)
        @test sort(model.ynames) == sort(MODEID_YNAMES)
        @test length(model.ynames) == length(MODEID_YNAMES)
    end

    @testset "ensemble has multiple models with consistent structure" begin
        model = TurbulentTransport.load_modeid_model(TEST_MODEID_MODEL)
        if model isa TurbulentTransport.ModeIDensemble
            @test length(model.models) > 1
            for m in model.models
                @test m isa TurbulentTransport.ModeIDmodel
                @test m.xnames == model.models[1].xnames
                @test m.ynames == model.models[1].ynames
            end
        end
    end

    @testset "getproperty delegation on ensemble" begin
        model = TurbulentTransport.load_modeid_model(TEST_MODEID_MODEL)
        if model isa TurbulentTransport.ModeIDensemble
            @test length(model.xnames) == MODEID_N_INPUTS
            @test length(model.ynames) == MODEID_N_CLASSES
        end
    end

    @testset "load_modeid_model_once is memoized (same object)" begin
        m1 = TurbulentTransport.load_modeid_model_once(TEST_MODEID_MODEL)
        m2 = TurbulentTransport.load_modeid_model_once(TEST_MODEID_MODEL)
        @test m1 === m2
    end

    @testset "show methods do not error" begin
        model = TurbulentTransport.load_modeid_model(TEST_MODEID_MODEL)
        @test !isempty(sprint(show, MIME"text/plain"(), model))
        if model isa TurbulentTransport.ModeIDensemble
            @test !isempty(sprint(show, MIME"text/plain"(), model.models[1]))
        end
    end
end

@testset "ModeID prediction" begin
    model = TurbulentTransport.load_modeid_model(TEST_MODEID_MODEL)
    x_mid = Float64.((model.xbounds[:, 1] .+ model.xbounds[:, 2]) ./ 2)

    @testset "predict_modeid vector input returns shape (N_classes,)" begin
        y = TurbulentTransport.predict_modeid(model, x_mid)
        @test length(y) == MODEID_N_CLASSES
    end

    @testset "predict_modeid matrix input returns shape (N_classes, N_samples)" begin
        X = hcat(x_mid, x_mid .* 1.05, x_mid .* 0.95)
        Y = TurbulentTransport.predict_modeid(model, X)
        @test size(Y) == (MODEID_N_CLASSES, 3)
    end

    @testset "predict_modeid outputs are all finite" begin
        y = TurbulentTransport.predict_modeid(model, x_mid)
        @test all(isfinite, y)

        X = hcat(x_mid, x_mid .* 1.05, x_mid .* 0.95)
        Y = TurbulentTransport.predict_modeid(model, X)
        @test all(isfinite, Y)
    end

    @testset "outputs are valid probabilities (non-negative, sum to 1)" begin
        y = TurbulentTransport.predict_modeid(model, x_mid)
        @test all(y .>= 0)
        @test sum(y) ≈ 1.0 atol=1e-5

        X = hcat(x_mid, x_mid .* 1.05, x_mid .* 0.95)
        Y = TurbulentTransport.predict_modeid(model, X)
        @test all(Y .>= 0)
        for j in axes(Y, 2)
            @test sum(Y[:, j]) ≈ 1.0 atol=1e-5
        end
    end

    @testset "vector and matrix predictions are consistent" begin
        y_vec = TurbulentTransport.predict_modeid(model, x_mid)
        X = hcat(x_mid, x_mid .* 1.05)
        Y_mat = TurbulentTransport.predict_modeid(model, X)
        @test Y_mat[:, 1] ≈ y_vec rtol=sqrt(eps(Float32))
    end

    @testset "predictions are deterministic" begin
        y1 = TurbulentTransport.predict_modeid(model, x_mid)
        y2 = TurbulentTransport.predict_modeid(model, x_mid)
        @test y1 == y2
    end

    @testset "different inputs give different outputs" begin
        y1 = TurbulentTransport.predict_modeid(model, x_mid)
        y2 = TurbulentTransport.predict_modeid(model, x_mid .* 1.5)
        @test !(y1 ≈ y2)
    end

    @testset "dominant mode is valid TurbulenceMode" begin
        y = TurbulentTransport.predict_modeid(model, x_mid)
        ynames = model.ynames
        dominant_idx = argmax(y)
        dominant_name = ynames[dominant_idx]
        @test dominant_name ∈ MODEID_YNAMES
    end
end

@testset "NNModeIdentification" begin
    model = TurbulentTransport.load_modeid_model(TEST_MODEID_MODEL)
    x_mid = Float64.((model.xbounds[:, 1] .+ model.xbounds[:, 2]) ./ 2)
    probs_vec = TurbulentTransport.predict_modeid(model, x_mid)

    @testset "construct NNModeIdentification manually" begin
        output_modes = [TurbulentTransport._YNAME_TO_MODE[yn] for yn in model.ynames]
        probs_dict = Dict{TurbulentTransport.TurbulenceMode, Float64}(
            mode => Float64(probs_vec[k]) for (k, mode) in enumerate(output_modes)
        )
        dominant = output_modes[argmax(probs_vec)]
        dominant_prob = maximum(probs_vec)

        mid = TurbulentTransport.NNModeIdentification{Float64}(probs_dict, dominant, dominant_prob)

        @test mid isa TurbulentTransport.NNModeIdentification
        @test mid.dominant_mode isa TurbulentTransport.TurbulenceMode
        @test 0.0 <= mid.dominant_mode_fraction <= 1.0
        @test length(mid.probabilities) == MODEID_N_CLASSES
        @test sum(values(mid.probabilities)) ≈ 1.0 atol=1e-5
    end

    @testset "show method does not error" begin
        output_modes = [TurbulentTransport._YNAME_TO_MODE[yn] for yn in model.ynames]
        probs_dict = Dict{TurbulentTransport.TurbulenceMode, Float64}(
            mode => Float64(probs_vec[k]) for (k, mode) in enumerate(output_modes)
        )
        dominant = output_modes[argmax(probs_vec)]
        dominant_prob = maximum(Float64.(probs_vec))

        mid = TurbulentTransport.NNModeIdentification{Float64}(probs_dict, dominant, dominant_prob)
        @test !isempty(sprint(show, MIME"text/plain"(), mid))
    end
end
