module TokamakEnergyTransportExt

using TurbulentTransport
using TokamakEnergyTransport

function __init__()
    # Register the model path when extension loads (prepend=true by default,
    # so provider models take precedence over built-ins on name collision)
    TurbulentTransport.register_model_path!(TokamakEnergyTransport.models_path())
    @debug "TokamakEnergyTransport models registered" path=TokamakEnergyTransport.models_path()
end

end # module
