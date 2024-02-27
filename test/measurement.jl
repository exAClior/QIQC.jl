using Test
using QIQC
using Yao

@testset "state distance" begin
	ρ1 = rand_density_matrix(1)
	ρ2 = rand_density_matrix(1)
	@test variational_distance(ρ1, ρ2, 0.001) ≈ tracedist(ρ1, ρ2)
end
