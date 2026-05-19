#!/usr/bin/env julia
# Driver script invoked by bin/qlgyro. Loads input.cgyro + input.qlgyro from
# the simdir, runs the modular Julia QLGYRO (linear CGYRO + TJLF saturation
# rules), and writes out.qlgyro.gbflux.

using TurbulentTransport
const TT = TurbulentTransport

function parse_args(args)
    simdir = pwd()
    n_mpi = 1
    n_omp = 1
    local_mode = nothing  # nothing = auto-detect
    for arg in args
        if startswith(arg, "--simdir=")
            simdir = arg[length("--simdir=")+1:end]
        elseif startswith(arg, "--n=")
            n_mpi = parse(Int, arg[length("--n=")+1:end])
        elseif startswith(arg, "--nomp=")
            n_omp = parse(Int, arg[length("--nomp=")+1:end])
        elseif arg == "--local"
            local_mode = true
        elseif arg == "--cluster"
            local_mode = false
        else
            error("qlgyro_driver: unknown argument: $arg")
        end
    end
    return simdir, n_mpi, n_omp, local_mode
end

simdir, n_mpi, n_omp, local_mode = parse_args(ARGS)

cgyro_path = joinpath(simdir, "input.cgyro")
qlgyro_path = joinpath(simdir, "input.qlgyro")
isfile(cgyro_path) || error("missing input.cgyro in $simdir")
isfile(qlgyro_path) || error("missing input.qlgyro in $simdir")

input_cgyro = TT.InputCGYRO()
TT.load(input_cgyro, cgyro_path)

input_qlgyro = TT.InputQLGYRO()
TT.load(input_qlgyro, qlgyro_path)

@info "qlgyro: running modular QLGYRO" simdir n_mpi n_omp NKY=input_qlgyro.NKY KYGRID_MODEL=input_qlgyro.KYGRID_MODEL SAT_RULE=input_qlgyro.SAT_RULE

sol = TT.run_qlgyro(input_cgyro, input_qlgyro;
                    basedir=simdir,
                    n_mpi=n_mpi,
                    n_omp=n_omp,
                    local_mode=local_mode)

@info "qlgyro: done" gbflux_path=joinpath(simdir, "out.qlgyro.gbflux")
println("Q_e = ", sol.ENERGY_FLUX_e)
println("Q_i = ", sol.ENERGY_FLUX_i)
println("Γ_e = ", sol.PARTICLE_FLUX_e)
println("Γ_i = ", sol.PARTICLE_FLUX_i)
println("Π_i = ", sol.STRESS_TOR_i)
