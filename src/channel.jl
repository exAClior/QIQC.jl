using Yao
using JuMP
using SCS: SCS
using LinearAlgebra

"""
Use Semidefinite Programming to find the Choi representation of a channel.
The channel maps ρ1 to σ1 and ρ2 to σ2.

Following https://shuvomoy.github.io/blogs/posts/Solving_semidefinite_programming_problems_in_Julia/
"""
function sdp_Choi_rep(ρ1::DensityMatrix, ρ2::DensityMatrix, σ1::DensityMatrix, σ2::DensityMatrix)
	N_A = size(ρ1, 1)
	N_B = size(σ1, 1)
	model = Model(SCS.Optimizer)
	set_silent(model)
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
