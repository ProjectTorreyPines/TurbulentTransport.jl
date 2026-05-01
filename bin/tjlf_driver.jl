#!/usr/bin/env julia
# Driver script invoked by bin/tjlf. Loads input.tglf (or input.tjlf) from the
# simdir, runs TJLF end-to-end, and writes out.tglf.gbflux.

using TurbulentTransport
const TT = TurbulentTransport
import TJLF

function parse_args(args)
    simdir = pwd()
    for arg in args
        if startswith(arg, "--simdir=")
            simdir = arg[length("--simdir=")+1:end]
        else
            error("tjlf_driver: unknown argument: $arg")
        end
    end
    return simdir
end

simdir = parse_args(ARGS)

tglf_path = joinpath(simdir, "input.tglf")
tjlf_path = joinpath(simdir, "input.tjlf")

input_file = if isfile(tglf_path)
    tglf_path
elseif isfile(tjlf_path)
    tjlf_path
else
    error("missing input.tglf or input.tjlf in $simdir")
end

@info "tjlf: loading inputs" simdir input_file
input_tglf = TT.InputTGLF{Float64}()
TT.load(input_tglf, input_file)
TT.apply_presets!(input_tglf)

@info "tjlf: running TJLF" NS=input_tglf.NS NKY=input_tglf.NKY KYGRID_MODEL=input_tglf.KYGRID_MODEL SAT_RULE=input_tglf.SAT_RULE
result = TJLF.run(input_tglf)

gbflux_path = joinpath(simdir, "out.tglf.gbflux")
TT.write_gbflux(result.QL_flux_out, gbflux_path)
@info "tjlf: done" gbflux_path

println("Q_e = ", TJLF.Qe(result.QL_flux_out))
println("Q_i = ", TJLF.Qi(result.QL_flux_out))
println("Γ_e = ", TJLF.Γe(result.QL_flux_out))
println("Γ_i = ", TJLF.Γi(result.QL_flux_out))
println("Π_i = ", TJLF.Πi(result.QL_flux_out))
