using Test
using Yao, QIQC
using LinearAlgebra

@testset "SDP Choi trivial" begin
	num_qubits = 1

	ρ1 = rand_density_matrix(num_qubits).state
	ρ2 = rand_density_matrix(num_qubits).state

	J1 = sdp_Choi_rep(ρ1, ρ2, ρ1, ρ2)

	@test isapprox(partial_tr(DensityMatrix(J1 * kron(ρ1', mat(igate(num_qubits)))), tuple((num_qubits + 1):(num_qubits + num_qubits)...)).state, ρ1, atol = 1e-6)

	@test isapprox(partial_tr(DensityMatrix(J1 * kron(ρ2', mat(igate(num_qubits)))), tuple((num_qubits + 1):(num_qubits + num_qubits)...)).state, ρ2, atol = 1e-6)
end

@testset "SDP Choi non-trivial" begin
	num_qubits = 1

	ρ1 = density_matrix(zero_state(num_qubits)).state
	ρ2 = density_matrix(ghz_state(num_qubits)).state

	σ1 = density_matrix(arrayreg(bit"1")).state
	σ2 = density_matrix((arrayreg(bit"0") - arrayreg(bit"1")) / sqrt(2)).state

	J1 = sdp_Choi_rep(ρ1, ρ2, σ1, σ2)

	@test isapprox(partial_tr(DensityMatrix(J1 * kron(ρ1', mat(igate(num_qubits)))), tuple((num_qubits + 1):(num_qubits + num_qubits)...)).state, σ1, atol = 1e-6)

	@test isapprox(partial_tr(DensityMatrix(J1 * kron(ρ2', mat(igate(num_qubits)))), tuple((num_qubits + 1):(num_qubits + num_qubits)...)).state, σ2, atol = 1e-6)
end
