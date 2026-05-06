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

include("qlgyro.jl")

include("utils.jl")

export InputTGLF, InputTJLF, available_models, available_qlnn_bundles, model_selector
export run_qlnn, loadqlnnbundle, loadqlnnmodel, QLNNmodel, QLNNensemble, QLNNbundle

const document = Dict()
document[Symbol(@__MODULE__)] = [name for name in Base.names(@__MODULE__; all=false, imported=false) if name != Symbol(@__MODULE__)]

end # module
