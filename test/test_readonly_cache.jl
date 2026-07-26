import SHA
import TurbulentTransport
using TurbulentTransport: ensure_model_file!, is_lfs_pointer

# The read-only scenario (container SIF, shared site install) only exists on
# POSIX systems, and Windows chmod cannot make a directory non-writable to
# simulate it — so this testset is POSIX-only.
if Sys.iswindows()
    @testset "Read-only model cache" begin
        @test_skip "read-only directories cannot be simulated on Windows"
    end
else
@testset "Read-only model cache" begin
    mktempdir() do dir
        # A fake model payload and a Git LFS pointer stub for it.
        payload = joinpath(dir, "payload.bin")
        write(payload, rand(UInt8, 4096))
        oid = open(io -> bytes2hex(SHA.sha256(io)), payload)
        stubdir = joinpath(dir, "models")
        mkpath(stubdir)
        stub = joinpath(stubdir, "fake_model.bson")
        write(stub, "version https://git-lfs.github.com/spec/v1\noid sha256:$oid\nsize $(filesize(payload))\n")

        downloads = Ref(0)
        TurbulentTransport._LFS_URL_OVERRIDE[] = (ref, rel) -> begin
            downloads[] += 1
            return "file://" * payload
        end
        try
            # Simulate a read-only install (container SIF / shared site dir).
            chmod(stubdir, 0o555)
            @test !TurbulentTransport._dir_writable(stubdir)

            resolved = ensure_model_file!(stub)
            @test resolved != stub                    # redirected to the cache
            @test read(resolved) == read(payload)     # correct bytes
            @test is_lfs_pointer(stub)                # stub left untouched
            n = downloads[]
            @test n >= 1

            resolved2 = ensure_model_file!(stub)      # second call: cache hit
            @test resolved2 == resolved
            @test downloads[] == n                    # no re-download

            # Writable again: materialization happens in place, path unchanged.
            chmod(stubdir, 0o755)
            resolved3 = ensure_model_file!(stub)
            @test resolved3 == stub
            @test read(stub) == read(payload)
            @test !is_lfs_pointer(stub)
        finally
            TurbulentTransport._LFS_URL_OVERRIDE[] = nothing
            chmod(stubdir, 0o755)
        end
    end
end
end
