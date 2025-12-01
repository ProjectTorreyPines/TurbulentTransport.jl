using Dates
import Flux

@testset "Model Loading" begin
    @testset "available_models" begin
        models = available_models()

        @test !isempty(models)
        @test models isa Vector{String}

        # Check for known models
        @test TEST_MODEL_ENSEMBLE in models
        @test "sat3_em_d3d_azf-1" in models
        @test "sat2_em_d3d_azf-1" in models
    end

    @testset "loadmodel ensemble structure" begin
        # All models in the repository are ensembles
        model = loadmodel("sat0_em_d3d")

        @test model isa TGLFNNensemble
        @test length(model.models) > 1

        # Test that ensemble getproperty delegates to first model
        @test !isempty(model.xnames)
        @test !isempty(model.ynames)
        @test length(model.xm) == length(model.xnames)
        @test length(model.xσ) == length(model.xnames)
        @test length(model.ym) == length(model.ynames)
        @test length(model.yσ) == length(model.ynames)
        @test size(model.xbounds, 1) == length(model.xnames)
        @test size(model.xbounds, 2) == 2  # min and max bounds
    end

    @testset "loadmodel ensemble" begin
        # sat3_em_d3d_azf-1 should be an ensemble
        model = loadmodel(TEST_MODEL_ENSEMBLE)

        @test model isa TGLFNNensemble
        @test length(model.models) > 1

        # Check that all models in ensemble have consistent structure
        first_model = model.models[1]
        for m in model.models
            @test length(m.xnames) == length(first_model.xnames)
            @test length(m.ynames) == length(first_model.ynames)
        end
    end

    @testset "TGLFNNmodel structure (from ensemble)" begin
        # Access individual model from ensemble
        ensemble = loadmodel("sat0_em_d3d")
        model = ensemble.models[1]

        @test model isa TGLFNNmodel
        @test model.fluxmodel isa Flux.Chain
        @test model.name isa String
        @test model.date isa Dates.DateTime
        @test model.nions isa Int
        @test model.nions >= 1
    end

    @testset "TGLFNNensemble getproperty" begin
        ensemble = loadmodel(TEST_MODEL_ENSEMBLE)

        # getproperty should delegate to first model for most fields
        @test ensemble.xnames == ensemble.models[1].xnames
        @test ensemble.ynames == ensemble.models[1].ynames
        @test ensemble.nions == ensemble.models[1].nions

        # fluxmodel should error
        @test_throws ErrorException ensemble.fluxmodel
    end

    @testset "loadmodel with .bson extension" begin
        model1 = loadmodel("sat0_em_d3d")
        model2 = loadmodel("sat0_em_d3d.bson")

        @test model1.xnames == model2.xnames
        @test model1.ynames == model2.ynames
    end

    @testset "loadmodel error for nonexistent model" begin
        @test_throws ErrorException loadmodel("nonexistent_model_xyz")
    end

    @testset "TGLFNNmodel show method" begin
        ensemble = loadmodel("sat0_em_d3d")
        model = ensemble.models[1]

        # Capture show output
        io = IOBuffer()
        show(io, MIME"text/plain"(), model)
        output = String(take!(io))

        # Verify key components are displayed
        @test contains(output, "TGLFNNmodel")
        @test contains(output, "date:")
        @test contains(output, "nions:")
        @test contains(output, "xnames")
        @test contains(output, "ynames")
        @test contains(output, string(length(model.xnames)))
        @test contains(output, string(length(model.ynames)))
    end

    @testset "TGLFNNensemble show method" begin
        ensemble = loadmodel(TEST_MODEL_ENSEMBLE)

        # Capture show output
        io = IOBuffer()
        show(io, MIME"text/plain"(), ensemble)
        output = String(take!(io))

        # Verify ensemble header
        @test contains(output, "TGLFNNensemble")
        @test contains(output, "n models:")
        @test contains(output, string(length(ensemble.models)))

        # Verify first model info is also shown
        @test contains(output, "TGLFNNmodel")
        @test contains(output, "xnames")
        @test contains(output, "ynames")
    end

    @testset "display works without error" begin
        ensemble = loadmodel("sat0_em_d3d")
        model = ensemble.models[1]

        # display() via show to IOBuffer should not throw
        io = IOBuffer()
        @test (show(io, MIME"text/plain"(), ensemble); true)
        @test (show(io, MIME"text/plain"(), model); true)

        # repr should also work
        @test !isempty(repr(MIME"text/plain"(), ensemble))
        @test !isempty(repr(MIME"text/plain"(), model))
    end
end

