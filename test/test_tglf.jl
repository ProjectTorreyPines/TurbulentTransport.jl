# Tests for the environment-only / display helpers in src/tglf.jl that do not
# require an installed TGLF binary (`run_tglf` shells out to `tglf` and is
# covered elsewhere / in integration).

@testset "tglf helpers" begin
    @testset "gacode_preamble" begin
        # Off NERSC/Perlmutter: empty preamble (GACODE assumed on PATH).
        withenv("NERSC_HOST" => "not_perlmutter") do
            @test TurbulentTransport.gacode_preamble() == ""
        end

        # On Perlmutter with GACODE_ROOT_CPU set: sources that root + platform.
        withenv("NERSC_HOST" => "perlmutter", "GACODE_ROOT_CPU" => "/opt/gacode") do
            pre = TurbulentTransport.gacode_preamble()
            @test occursin("export GACODE_ROOT=/opt/gacode", pre)
            @test occursin("PERLMUTTER_CPU", pre)
            @test occursin("gacode_setup", pre)
        end

        # On Perlmutter without GACODE_ROOT_CPU: informative error.
        withenv("NERSC_HOST" => "perlmutter", "GACODE_ROOT_CPU" => nothing) do
            @test_throws ErrorException TurbulentTransport.gacode_preamble()
        end
    end

    @testset "Base.show(InputTGLF)" begin
        input = load_sample_input()
        io = IOBuffer()
        show(io, MIME"text/plain"(), input)
        str = String(take!(io))
        @test occursin("NS = 3", str)
        @test occursin("SAT_RULE = 3", str)
        # Species fields with index > NS (=3) are suppressed.
        @test !occursin("MASS_4", str)
        # A present low-index species field is shown.
        @test occursin("MASS_1", str)
    end
end
