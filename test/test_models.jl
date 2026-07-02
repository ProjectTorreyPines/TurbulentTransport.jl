using Dates
import Flux
import SHA

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

    @testset "register_model_path!" begin
        saved = copy(TurbulentTransport._MODEL_SEARCH_PATHS)
        try
            mktempdir() do tmpdir
                TurbulentTransport.register_model_path!(tmpdir)
                @test first(TurbulentTransport._MODEL_SEARCH_PATHS) == tmpdir

                # prepend=false appends to end
                tmpdir2 = mktempdir()
                TurbulentTransport.register_model_path!(tmpdir2; prepend=false)
                @test last(TurbulentTransport._MODEL_SEARCH_PATHS) == tmpdir2
            end
        finally
            empty!(TurbulentTransport._MODEL_SEARCH_PATHS)
            append!(TurbulentTransport._MODEL_SEARCH_PATHS, saved)
        end
    end

    @testset "resolve_model_path" begin
        # Absolute path to existing file
        path = TurbulentTransport.resolve_model_path("sat0_em_d3d")
        @test isfile(path)
        @test endswith(path, ".bson")

        # Direct file path resolves immediately
        @test TurbulentTransport.resolve_model_path(path) == path

        # Nonexistent model throws
        @test_throws ErrorException TurbulentTransport.resolve_model_path("nonexistent_zzz")
    end

    @testset "Git LFS pointer detection" begin
        mktempdir() do tmpdir
            pointer = joinpath(tmpdir, "pointer.bson")
            write(pointer, "version https://git-lfs.github.com/spec/v1\noid sha256:abc\nsize 1\n")
            real = joinpath(tmpdir, "real.bson")
            write(real, 0x01)

            @test TurbulentTransport.is_lfs_pointer(pointer)
            @test !TurbulentTransport.is_lfs_pointer(real)
        end
    end

    @testset "LFS pointer parsing" begin
        mktempdir() do tmpdir
            good = joinpath(tmpdir, "good.bson")
            write(good,
                "version https://git-lfs.github.com/spec/v1\n" *
                "oid sha256:0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef\n" *
                "size 42\n")
            info = TurbulentTransport._lfs_pointer_info(good)
            @test info !== nothing
            @test info.oid == "0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef"
            @test info.size == 42

            bad = joinpath(tmpdir, "bad.bson")
            write(bad, "version https://git-lfs.github.com/spec/v1\nbogus\n")
            @test TurbulentTransport._lfs_pointer_info(bad) === nothing
        end
    end

    @testset "_sha256_of_file matches SHA.sha256" begin
        mktempdir() do tmpdir
            blob = rand(UInt8, 1024)
            f = joinpath(tmpdir, "blob.bin")
            write(f, blob)
            @test TurbulentTransport._sha256_of_file(f) == bytes2hex(SHA.sha256(blob))
        end
    end

    @testset "_lfs_media_url default (GitHub) vs override" begin
        old = TurbulentTransport._LFS_URL_OVERRIDE[]
        try
            # With no override, builds the media.githubusercontent.com LFS URL.
            TurbulentTransport._LFS_URL_OVERRIDE[] = nothing
            url = TurbulentTransport._lfs_media_url("abc123", "models/foo/bar.pt")
            @test occursin("media.githubusercontent.com", url)
            @test occursin("ProjectTorreyPines/TurbulentTransport.jl", url)
            @test endswith(url, "abc123/models/foo/bar.pt")

            # With an override installed, it is used verbatim instead.
            TurbulentTransport._LFS_URL_OVERRIDE[] = (ref, rel) -> "file:///$(ref)/$(rel)"
            @test TurbulentTransport._lfs_media_url("r", "p") == "file:///r/p"
        finally
            TurbulentTransport._LFS_URL_OVERRIDE[] = old
        end
    end

    @testset "ensure_model_file! materializes via SHA-verified override" begin
        # Use a `file://` URL override so the test is hermetic (no network).
        # Two refs are exposed: "bad_ref" serves wrong bytes, "good_ref" serves
        # the real payload. The cascade tries both; only the SHA-matching one
        # may replace the pointer in place.
        mktempdir() do tmpdir
            payload = rand(UInt8, 257)
            oid = bytes2hex(SHA.sha256(payload))
            wrong = copy(payload)
            wrong[end] ⊻= 0xFF

            target = joinpath(tmpdir, "foo.bson")
            write(target,
                "version https://git-lfs.github.com/spec/v1\n" *
                "oid sha256:$oid\n" *
                "size $(length(payload))\n")
            @test TurbulentTransport.is_lfs_pointer(target)

            good_blob = joinpath(tmpdir, "good.blob")
            bad_blob = joinpath(tmpdir, "bad.blob")
            write(good_blob, payload)
            write(bad_blob, wrong)

            url_map = Dict("good_ref" => "file://" * good_blob,
                           "bad_ref"  => "file://" * bad_blob)
            TurbulentTransport._LFS_URL_OVERRIDE[] = (ref, _rel) ->
                get(url_map, String(ref), "file:///does/not/exist")

            try
                # Only bad ref available → SHA mismatch on every candidate → error.
                ENV["TURBULENTTRANSPORT_MODELS_REF"] = "bad_ref"
                @test_throws ErrorException TurbulentTransport.ensure_model_file!(target)
                @test TurbulentTransport.is_lfs_pointer(target)  # untouched

                # Good ref available → materializes, file replaced in place.
                ENV["TURBULENTTRANSPORT_MODELS_REF"] = "good_ref"
                returned = TurbulentTransport.ensure_model_file!(target)
                @test returned == target
                @test !TurbulentTransport.is_lfs_pointer(target)
                @test read(target) == payload
                @test TurbulentTransport._sha256_of_file(target) == oid

                # Idempotent: real file → no-op.
                @test TurbulentTransport.ensure_model_file!(target) == target
            finally
                delete!(ENV, "TURBULENTTRANSPORT_MODELS_REF")
                TurbulentTransport._LFS_URL_OVERRIDE[] = nothing
            end
        end
    end
end

