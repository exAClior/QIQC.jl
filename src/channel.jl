using Yao
using JuMP
using SCS: SCS
using LinearAlgebra

import Yao.YaoAPI.partial_tr as partial_tr
import Yao.YaoArrayRegister.partial_tr! as partial_tr!

"""
Partial assuming the system is always qubits
"""
function partial_tr(dm::AbstractMatrix{T}, locs; D::Int = 2) where T <: AbstractJuMPScalar
	locs = Tuple(locs)
	nbits = Yao.log2i(Base.size(dm, 1))
	@assert_locs_safe nbits (locs...,)
	m = nbits - length(locs)
	strides = ntuple(i -> D^(i - 1), nbits)
	out_strides = ntuple(i -> D^(i - 1), m)
	remainlocs = (setdiff(1:nbits, locs)...,)
	remain_strides = map(i -> strides[i], remainlocs)
	trace_strides = ntuple(i -> strides[locs[i]], length(locs))
	state = similar(dm, D^m, D^m)   # NOTE: does not preserve adjoint
	fill!(state, zero(T))
	partial_tr!(Val{D}(), state, dm, trace_strides, out_strides, remain_strides)
	return state
end

"""
Use Semidefinite Programming to find the Choi representation of a channel.
The channel maps ρ1 to σ1 and ρ2 to σ2.

Following https://shuvomoy.github.io/blogs/posts/Solving_semidefinite_programming_problems_in_Julia/

# TODO
This is the naive version without considering c1 and c2 which
always gives a solution and the validity of the solution is 
guaranteed when only c1 = 1.0 and c2 = 0.0.
"""
function sdp_Choi_rep(ρ1::AbstractMatrix, ρ2::AbstractMatrix, σ1::AbstractMatrix, σ2::AbstractMatrix;
	optimizer = SCS.Optimizer, silent = true)
	model = Model(optimizer)
	silent && set_silent(model)

	N_A = Base.size(ρ1, 1)
	N_B = Base.size(σ1, 1)

	A_qubits = Yao.log2i(N_A)
	B_qubits = Yao.log2i(N_B)

	# this is a CP map from A'A to A'B
	J1 = @variable(model, [1:(N_A * N_B), 1:(N_A * N_B)] in HermitianPSDCone())

	@objective(model, Min, 1.0)

	# trace out the system B
	@constraint(model, partial_tr(J1, tuple(1:A_qubits...)) .== LinearAlgebra.I)

	@constraint(model, partial_tr(J1 * kron(ρ1', mat(igate(B_qubits))), tuple((A_qubits + 1):(A_qubits + B_qubits)...)) .== σ1)

	@constraint(model, partial_tr(J1 * kron(ρ2', mat(igate(B_qubits))), tuple((A_qubits + 1):(A_qubits + B_qubits)...)) .== σ2)

	optimize!(model)

	is_solved_and_feasible(model) && return JuMP.value.(J1)

	return error("SDP failed to find a solution")
end
