# Tests for the generic SLURM submit/poll helpers in src/slurm_utils.jl.
# These functions are pure (parsing / path building / (de)serialization) and do
# not require an actual SLURM scheduler, so they run anywhere.

@testset "slurm_utils" begin
    @testset "parse_slurm_jobid" begin
        @test TurbulentTransport.parse_slurm_jobid("Submitted batch job 12345678") == "12345678"
        @test TurbulentTransport.parse_slurm_jobid("Submitted batch job 987654 on cluster perlmutter") == "987654"
        # Leading text with a short number (<5 digits) should NOT match; only the 5+ digit id.
        @test TurbulentTransport.parse_slurm_jobid("job 42 -> Submitted batch job 55555") == "55555"
        @test TurbulentTransport.parse_slurm_jobid("no job id here") == ""
        @test TurbulentTransport.parse_slurm_jobid("") == ""
    end

    @testset "check_slurm_status empty id" begin
        # An empty job id short-circuits to :unknown without shelling out to sacct.
        @test TurbulentTransport.check_slurm_status("") == :unknown
    end

    # The fake `sacct` is a bash shim and SLURM is unix-only, so skip on Windows.
    Sys.iswindows() || @testset "check_slurm_status state mapping (fake sacct)" begin
        # Map every documented sacct State to the driver's status symbol using a
        # fake `sacct` on PATH, so no real scheduler is required.
        for (state, expected) in (
            ("PENDING", :pending),
            ("RUNNING", :running),
            ("CONFIGURING", :running),
            ("COMPLETING", :running),
            ("COMPLETED", :completed),
            ("FAILED", :failed),
            ("TIMEOUT", :failed),
            ("CANCELLED", :failed),
            ("NODE_FAIL", :failed),
            ("OUT_OF_MEMORY", :failed),
            ("CANCELLED by 12345", :failed),  # startswith("CANCELLED") branch
            ("SOME_NEW_STATE", :unknown),     # unrecognized -> :unknown
            ("", :unknown),                   # empty sacct output -> :unknown
        )
            with_fake_sacct(state) do
                @test TurbulentTransport.check_slurm_status("12345") == expected
            end
        end

        # A failing `sacct` (non-zero exit) is swallowed by the catch -> :unknown.
        with_fake_sacct("__FAIL__") do
            @test TurbulentTransport.check_slurm_status("12345") == :unknown
        end
    end

    @testset "default_results_dir" begin
        mktempdir() do tmp
            dir = TurbulentTransport.default_results_dir("MYSUB"; prefix=tmp)
            # <prefix>/MYSUB/<runprefix>_<user>_<timestamp>
            @test startswith(dir, joinpath(tmp, "MYSUB"))
            @test occursin("mysub_", basename(dir))  # default runprefix is lowercase(subdir)
            # The base <prefix>/<subdir> directory is created.
            @test isdir(joinpath(tmp, "MYSUB"))
        end
    end

    @testset "default_results_dir custom runprefix" begin
        mktempdir() do tmp
            dir = TurbulentTransport.default_results_dir("SUB"; runprefix="custompref", prefix=tmp)
            @test occursin("custompref_", basename(dir))
        end
    end

    @testset "default_results_dir falls back when prefix not writable" begin
        # A prefix that cannot be created (nested under a file) forces the PSCRATCH/tempdir fallback.
        mktempdir() do tmp
            blocker = joinpath(tmp, "not_a_dir")
            write(blocker, "x")  # a regular file where a directory would be needed
            withenv("PSCRATCH" => tmp, "SCRATCH" => nothing) do
                dir = TurbulentTransport.default_results_dir("FALLBACKSUB"; prefix=joinpath(blocker, "cannot"))
                @test occursin(joinpath("results", "FUSE", "FALLBACKSUB"), dir)
                @test startswith(dir, tmp)
            end
        end
    end

    @testset "save_state / load_state round trip" begin
        mktempdir() do tmp
            state = (a=1, b="hello", c=[1.0, 2.0, 3.0])
            path = joinpath(tmp, "nested", "state.jls")
            returned = TurbulentTransport.save_state(state, path)
            @test returned == path
            @test isfile(path)

            loaded = TurbulentTransport.load_state(path)
            @test loaded == state
            @test loaded.b == "hello"
            @test loaded.c == [1.0, 2.0, 3.0]
        end
    end

    @testset "load_state missing file returns nothing" begin
        mktempdir() do tmp
            @test TurbulentTransport.load_state(joinpath(tmp, "does_not_exist.jls")) === nothing
        end
    end
end
