using Test
using QIQC
using Yao

@testset "trace distance via measurement" begin
	ρ1 = rand_density_matrix(1)
	ρ2 = rand_density_matrix(1)

	@test isapprox(sdp_measurement(ρ1.state, ρ2.state), tracedist(ρ1, ρ2) / 2.0, atol = 0.00001)
	sdp_measurement(ρ1.state, ρ2.state)
	tracedist(ρ1, ρ2) / 2.0
	@time sdp_measurement(ρ1.state, ρ2.state)
end

# @testset "state distance" begin
# 	ρ1 = rand_density_matrix(1)
# 	ρ2 = rand_density_matrix(1)
# 	@test variational_distance(ρ1, ρ2, 0.001) ≈ tracedist(ρ1, ρ2)
# end
