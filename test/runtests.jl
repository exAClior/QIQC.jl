using QIQC
using Test

@testset "measurements" begin
	include("measurement.jl")
end

@testset "quantum channel" begin
	include("channel.jl")
end