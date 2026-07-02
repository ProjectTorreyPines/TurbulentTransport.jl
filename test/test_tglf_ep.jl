# Tests for src/tglf_ep.jl.
#
# Two groups that do not need a running TGLF/TJLFEP or a SLURM scheduler:
#   1. The pure GACODE-style numeric helpers (`expro_bound_deriv`,
#      `expro_log_gradients`) — exact 3-point Lagrange derivatives.
#   2. The `run_tjlfep` SLURM driver in dry-run mode (`submit=false`), which only
#      renders a batch script and builds/persists a `TGLFEPRunState`, plus the
#      state pollers (`refresh_tjlfep!`, `tjlfep_status`, `load_tjlfep_results`)
#      and the sysimage-freshness check (`_tjlfep_sysimage_ok`).

@testset "tglf_ep" begin
    @testset "expro_bound_deriv" begin
        r = collect(0.0:0.25:2.0)

        # Exact for linear f = 2r + 3 -> f' = 2 everywhere.
        f_lin = 2 .* r .+ 3
        df_lin = TurbulentTransport.expro_bound_deriv(f_lin, r)
        @test all(≈(2.0), df_lin)

        # Exact for quadratic f = r^2 -> f' = 2r (3-point Lagrange is exact to deg 2).
        f_quad = r .^ 2
        df_quad = TurbulentTransport.expro_bound_deriv(f_quad, r)
        @test df_quad ≈ 2 .* r atol = 1e-9

        # Integer input is accepted (converted to Float64 internally).
        @test TurbulentTransport.expro_bound_deriv([1, 3, 5, 7], [0, 1, 2, 3]) ≈ fill(2.0, 4) atol = 1e-9

        # Length mismatch errors.
        @test_throws ErrorException TurbulentTransport.expro_bound_deriv([1.0, 2.0, 3.0], [0.0, 1.0])
    end

    @testset "expro_log_gradients" begin
        r = collect(0.0:0.2:2.0)
        # ni = exp(-2r) -> ln(ni) = -2r -> dlnnidr = -d(ln n)/dr = 2
        # ti = exp(3r)  -> ln(ti) =  3r -> dlntidr = -d(ln t)/dr = -3
        ni = exp.(-2 .* r)
        ti = exp.(3 .* r)
        dlnnidr, dlntidr = TurbulentTransport.expro_log_gradients(ni, ti, r)
        @test dlnnidr ≈ fill(2.0, length(r)) atol = 1e-8
        @test dlntidr ≈ fill(-3.0, length(r)) atol = 1e-8

        @test_throws ErrorException TurbulentTransport.expro_log_gradients([1.0, 2.0], [1.0, 2.0, 3.0], [1.0, 2.0, 3.0])
    end

    @testset "_git_sha / _tjlfep_sysimage_ok" begin
        mktempdir() do tmp
            # Non-git directories: _git_sha swallows the error and returns "".
            root1 = joinpath(tmp, "root1")
            root2 = joinpath(tmp, "root2")
            mkpath(root1)
            mkpath(root2)

            # Missing/empty sysimage path -> false.
            @test TurbulentTransport._tjlfep_sysimage_ok("", root1, root2) == false
            @test TurbulentTransport._tjlfep_sysimage_ok(joinpath(tmp, "nope.so"), root1, root2) == false

            # Existing sysimage but no .sha sidecar -> false.
            sysimg = joinpath(tmp, "img.so")
            write(sysimg, "not-a-real-sysimage")
            @test TurbulentTransport._tjlfep_sysimage_ok(sysimg, root1, root2) == false

            # Sidecar matching the current SHAs -> true (compute expected via the
            # same helper so this is robust whether or not the roots are git repos).
            want = "TJLFEP=$(TurbulentTransport._git_sha(root1))\nTJLF=$(TurbulentTransport._git_sha(root2))"
            write(sysimg * ".sha", want)
            @test TurbulentTransport._tjlfep_sysimage_ok(sysimg, root1, root2) == true

            # Mismatched sidecar -> false.
            write(sysimg * ".sha", "TJLFEP=deadbeef\nTJLF=cafef00d")
            @test TurbulentTransport._tjlfep_sysimage_ok(sysimg, root1, root2) == false
        end
    end

    @testset "run_tjlfep dry run (threads template)" begin
        mktempdir() do tmp
            state = run_tjlfep(:ITER; submit=false, gpu=false, basedir=tmp, nodes=5, n_scan=20)
            @test state isa TurbulentTransport.TGLFEPRunState
            @test state.status == :rendered
            @test state.job_id == ""
            @test state.build_job_id == ""
            @test state.case == :ITER
            @test isfile(state.batchfile)

            content = read(state.batchfile, String)
            @test occursin("CASE=\"ITER\"", content)
            @test occursin("SCAN_N=\"20\"", content)
            @test occursin("#SBATCH -N 5", content)
            # threads template is NOT the 3-phase SPMD one.
            @test !occursin("phase 1: prepare", content)

            # State was serialized next to the batch script.
            @test isfile(joinpath(tmp, ".tjlfep_state.jls"))
        end
    end

    @testset "run_tjlfep dry run (mps_team SPMD template)" begin
        mktempdir() do tmp
            state = run_tjlfep(:D3D; submit=false, gpu=false, basedir=tmp, inner=:mps_team, mps_team=8)
            @test state.status == :rendered
            @test state.case == :D3D
            content = read(state.batchfile, String)
            @test occursin("CASE=\"D3D\"", content)
            # SPMD template renders the three phases.
            @test occursin("phase 1: prepare", content)
            @test occursin("phase 2", content)
            @test occursin("phase 3: merge", content)
            @test occursin("INNER=\"mps_team\"", content)
        end
    end

    @testset "state pollers on a rendered (un-submitted) run" begin
        mktempdir() do tmp
            state = run_tjlfep(:ITER; submit=false, gpu=false, basedir=tmp)

            # refresh on an empty job id keeps status :rendered (no sacct call).
            refreshed = TurbulentTransport.refresh_tjlfep!(state)
            @test refreshed === state
            @test state.status == :rendered

            status_str = TurbulentTransport.tjlfep_status(state)
            @test occursin("TJLFEP", status_str)
            @test occursin("ITER", status_str)

            # No results file yet -> nothing.
            @test TurbulentTransport.load_tjlfep_results(state) === nothing
        end
    end

    # The fake `sacct` is a bash shim and SLURM is unix-only, so skip on Windows.
    Sys.iswindows() || @testset "state pollers on a submitted run (fake sacct)" begin
        mktempdir() do tmp
            state = run_tjlfep(:ITER; submit=false, gpu=false, basedir=tmp)
            state.job_id = "12345"  # pretend it was submitted

            # A resolvable %j log file next to the run gets picked up on refresh.
            logfile = joinpath(tmp, "tjlfep_12345.out")
            write(logfile, "some log line\nTIMING_RESULT phase=total_job seconds=42\n")

            with_fake_sacct("COMPLETED") do
                refreshed = TurbulentTransport.refresh_tjlfep!(state)
                @test refreshed === state
                @test state.status == :completed
                @test state.logfile == logfile
            end
            # refresh persisted the updated state.
            @test isfile(joinpath(tmp, ".tjlfep_state.jls"))

            # A finished job with a TIMING_RESULT/total_job line appends timing.
            status_str = TurbulentTransport.tjlfep_status(state)
            @test occursin("[completed]", status_str)
            @test occursin("TIMING_RESULT", status_str)
            @test occursin("phase=total_job", status_str)

            # Results loader returns the serialized NamedTuple plus dd_json path.
            results = (rho_scan=[0.5, 0.6], SFmin=[1.0, 2.0], width=[0.3, 0.4],
                       kymark=[0.1, 0.2], n_EP=[1.0], p_EP=[2.0])
            TurbulentTransport.save_state(results, joinpath(tmp, "tjlfep_results.jls"))
            write(joinpath(tmp, "dd_out.json"), "{}")

            loaded = TurbulentTransport.load_tjlfep_results(state)
            @test loaded.rho_scan == results.rho_scan
            @test loaded.SFmin == results.SFmin
            @test loaded.dd_json == joinpath(tmp, "dd_out.json")
        end
    end

    Sys.iswindows() || @testset "refresh_tjlfep! failed job (fake sacct)" begin
        mktempdir() do tmp
            state = run_tjlfep(:D3D; submit=false, gpu=false, basedir=tmp)
            state.job_id = "99999"
            with_fake_sacct("FAILED") do
                TurbulentTransport.refresh_tjlfep!(state)
                @test state.status == :failed
            end
        end
    end
end
