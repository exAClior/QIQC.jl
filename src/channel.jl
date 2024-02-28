using Yao
using JuMP
using SCS: SCS
using LinearAlgebra

import Yao.YaoAPI.partial_tr as partial_tr
import Yao.YaoArrayRegister.partial_tr! as partial_tr!

"""
Partial assuming the system is always qubits
"""
function partial_tr(dm::AbstractMatrix{T}, locs) where T <: AbstractJuMPScalar
	D = 2
	locs = Tuple(locs)
	nbits = Yao.log2i(size(dm, 1))
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

ρ1 = rand_density_matrix(1).state
ρ2 = rand_density_matrix(1).state
σ1 = rand_density_matrix(1).state
σ2 = rand_density_matrix(1).state

model = Model(SCS.Optimizer)
set_silent(model)

N_A = size(ρ1, 1)
N_B = size(σ1, 1)
A_qubits = Yao.log2i(N_A)
B_qubits = Yao.log2i(N_B)

@variable(model, c1 >= 0)
@variable(model, c2 >= 0)
J1 = @variable(model, [1:(N_A * N_B), 1:(N_A * N_B)] in HermitianPSDCone())
J2 = @variable(model, [1:(N_A * N_B), 1:(N_A * N_B)] in HermitianPSDCone())

@objective(model, Min, c1 + c2)

set_start_value(c1, 0.5)
set_start_value(c2, 0.5)

@constraint(model, partial_tr(J1, tuple(1:A_qubits...)) .== mat(igate(A_qubits)))
@constraint(model, partial_tr(J2, tuple(1:A_qubits...)) .== mat(igate(A_qubits)))

@constraint(model, partial_tr((c1 .* J1 - c2 .* J2) * kron(ρ1, mat(igate(B_qubits))), tuple((A_qubits + 1):(A_qubits + B_qubits)...)) .== σ1)
@constraint(model, partial_tr((c1 .* J1 - c2 .* J2) * kron(ρ2, mat(igate(B_qubits))), tuple((A_qubits + 1):(A_qubits + B_qubits)...)) .== σ2)

# need another optimizer or change the condition
optimize!(model)

status = JuMP.termination_status(model)
c1_sol = JuMP.value.(c1)
c2_sol = JuMP.value.(c2)
obj_value = JuMP.objective_value(model)

"""
Use Semidefinite Programming to find the Choi representation of a channel.
The channel maps ρ1 to σ1 and ρ2 to σ2.

Following https://shuvomoy.github.io/blogs/posts/Solving_semidefinite_programming_problems_in_Julia/
"""
function sdp_Choi_rep(ρ1::AbstractMatrix, ρ2::AbstractMatrix, σ1::AbstractMatrix, σ2::AbstractMatrix;
	optimizer = SCS.Optimizer, silent = true)
	model = Model(SCS.Optimizer)
	silent && set_silent(model)

	N_A = size(ρ1, 1)
	N_B = size(σ1, 1)
	@variable(model, c1 >= 0)
	@variable(model, c2 >= 0)
	@variable(model, J1[1:(N_A * N_B), 1:(N_A * N_B)], PSD)
	@variable(model, J2[1:(N_A * N_B), 1:(N_A * N_B)], PSD)

	@objective(model, Min, c1 + c2)

	set_start_value(c1, 0.5)
	set_start_value(c2, 0.5)

	for i in 1:(N_A * N_B)
		set_start_value(J1[i, i], 1.0)
		set_start_value(J2[i, i], 1.0)
	end

	@constraint(model, partial_tr(J1) == mat(igate(Yao.log2i(N_A))))
	@constraint(model, partial_tr(J2) == mat(igate(Yao.log2i(N_A))))

	@constraint(model, partial_tr((c1 .* J1 - c2 .* J2) * kron(ρ1, mat(igate(Yao.log2i(N_B))))) == σ1)
	@constraint(model, partial_tr((c1 .* J1 - c2 .* J2) * kron(ρ2, mat(igate(Yao.log2i(N_B))))) == σ2)
	optimize!(model)

	status = JuMP.termination_status(model)
	c1_sol = JuMP.value.(c1)
	c2_sol = JuMP.value.(c2)
	obj_value = JuMP.objective_value(model)

	return status, c1_sol, c2_sol, obj_value
end
