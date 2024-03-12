using Test
using QIQC
using Yao
using JuMP

@testset "entanglement" begin
	dit = 2
	ρ_init = density_matrix(rand_state(2; nlevel = dit)).state
	for _ in 1:10
		As = nothing
		Bs = nothing
		obj_val1, As = chsh_op_povm("As", nothing, nothing, ρ_init, dit)
		obj_val2, Bs = chsh_op_povm("As", nothing, As, ρ_init, dit)
		obj_val3, ρ = chsh_op_povm("ρ", As, Bs, nothing, dit)
		println(obj_val1, " ", obj_val2, " ", obj_val3)
	end

	@test isapprox(obj_val2, 2 * sqrt(2), atol = 1e-4)
end

@testset "entanglement" begin
	dit = 3
	ρ_init = density_matrix(rand_state(2; nlevel = dit)).state
	for _ in 1:10
		As = nothing
		Bs = nothing
		obj_val1, As = chsh_op_povm("As", nothing, nothing, ρ_init, dit)
		obj_val2, Bs = chsh_op_povm("As", nothing, As, ρ_init, dit)
		obj_val3, ρ = chsh_op_povm("ρ", As, Bs, nothing, dit)
		println(obj_val1, " ", obj_val2, " ", obj_val3)
	end

	@test isapprox(obj_val2, 2 * sqrt(2), atol = 1e-4)
end

@testset "entanglement" begin
	dit = 5
	ρ_init = density_matrix(rand_state(2; nlevel = dit)).state
	for _ in 1:10
		As = nothing
		Bs = nothing
		obj_val1, As = chsh_op_povm("As", nothing, nothing, ρ_init, dit)
		obj_val2, Bs = chsh_op_povm("As", nothing, As, ρ_init, dit)
		obj_val3, ρ = chsh_op_povm("ρ", As, Bs, nothing, dit)
		println(obj_val1, " ", obj_val2, " ", obj_val3)
	end

	using BitBasis

	@test isapprox(obj_val2, 2 * sqrt(2), atol = 1e-4)
end
