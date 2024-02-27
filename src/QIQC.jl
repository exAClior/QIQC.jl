module QIQC

using Yao, QuantumInformation

export sdp_Choi_rep
include("channel.jl")

export variational_distance
include("measurement.jl")
end
