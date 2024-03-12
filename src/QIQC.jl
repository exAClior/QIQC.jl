module QIQC

using Yao
using MLStyle

export sdp_Choi_rep
include("channel.jl")

export variational_distance, sdp_measurement
include("measurement.jl")

include("registers.jl")

export chsh_op_povm
include("entanglement.jl")
end
