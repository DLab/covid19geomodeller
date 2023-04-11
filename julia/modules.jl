module SIRAbstractSpace

include("./SIR.jl")
include("./AbstractSpace.jl")

end #module

module SIRGraphSpace

include("./SIR.jl")
include("./GraphSpace.jl")

end #module

module SEIRAbstractSpace

include("./SEIR.jl")
include("./AbstractSpace.jl")

end #module

module SEIRGraphSpace

include("./SEIR.jl")
include("./GraphSpace.jl")

end #module


module SEIRHVDAbstractSpace

include("./SEIRHVD.jl")
include("./AbstractSpace.jl")

end #module

module SEIRHVDGraphSpace

include("./SEIRHVD.jl")
include("./GraphSpace.jl")

end #module
