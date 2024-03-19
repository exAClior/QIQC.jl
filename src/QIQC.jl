module QIQC

using Yao
using MLStyle

export sdp_Choi_rep, partial_T
include("channel.jl")

export variational_distance, sdp_measurement
include("measurement.jl")

include("registers.jl")

export chsh_op_povm, cabxy
include("entanglement.jl")

export nonsignaling_game
include("nonsignaling.jl")
end
