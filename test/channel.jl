using Test
using Yao, QIQC
using LinearAlgebra

@testset "SDP Choi" begin
	num_qubits = 1

	ρ1 = rand_density_matrix(num_qubits).state
	ρ2 = rand_density_matrix(num_qubits).state
	σ1 = rand_density_matrix(num_qubits).state
	σ2 = rand_density_matrix(num_qubits).state

	J1 = sdp_Choi_rep(ρ1, ρ2, ρ1, ρ2)

	@test J1 ≈ mat(igate(2^num_qubits)) 
end
