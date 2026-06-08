module TurbulentTransport

using IMAS
using IMASutils: argmin_abs
using AdaptiveArrayPools
import TJLF
import TJLF: InputTGLF, InputTJLF
import GACODE

include("models.jl")

include("tglf.jl")

include("input_tglfs.jl")

include("tjlf.jl")

include("pooled_layers.jl") 

include("tglf_nn.jl")

include("qlnn.jl")

include("finn.jl")

include("modeID.jl")

include("cgyro.jl")

include("slurm_utils.jl")

include("qlgyro.jl")

include("utils.jl")

include("tglf_ep.jl")

export InputTGLF, InputTJLF, available_models, available_qlnn_bundles, model_selector
export InputTGLFEP, run_tjlfep
export run_qlnn, qlnn_fluctuation_spectra, loadqlnnbundle, loadqlnnmodel
export QLNNmodel, QLNNensemble, QLNNbundle
export qlnn_to_gpu, qlnn_fluctuation_spectra_gpu

const document = Dict()
document[Symbol(@__MODULE__)] = [name for name in Base.names(@__MODULE__; all=false, imported=false) if name != Symbol(@__MODULE__)]

# Runtime initialization. Anything that depends on `Threads.nthreads()` must
# happen here, not at top-level / precompile time, since precompilation runs
# single-threaded and would otherwise bake in a 1-thread fallback.
function __init__()
    init_qlnn_thread_pools!()
    return nothing
end

end # module
