module QIQC

using Yao

export sdp_Choi_rep
include("channel.jl")

export variational_distance, sdp_measurement
include("measurement.jl")
end
