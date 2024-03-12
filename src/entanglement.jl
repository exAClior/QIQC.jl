using JuMP
using Yao
using LinearAlgebra
using SCS

function cabxy(a, b, x, y)
	return (-1)^(a + b + x * y)
end

function chsh_op_state(optimizer = SCS.Optimizer, silent = true)
	model = Model(optimizer)
	silent && set_silent(model)

	ρ = density_matrix(ghz_state(2)).state

	v0 = rand_state(1)
	v1 = rand_state(1)
	B00 = v0.state * v0.state'
	B10 = I - B00
	B01 = v1.state * v1.state'
	B11 = I - B01

	Bs = [B00, B01, B10, B11]

	As = [@variable(model, [1:2, 1:2] in HermitianPSDCone()) for _ in 1:4]

	@constraint(model, As[1] + As[3] == LinearAlgebra.I)
	@constraint(model, As[2] + As[4] == LinearAlgebra.I)

	@objective(model, Max, real(sum([cabxy(a, b, x, y) * tr(kron(As[(a) * 2 + x + 1], Bs[(b) * 2 + y + 1]) * ρ) for a in 0:1, b in 0:1, x in 0:1, y in 0:1])))

	optimize!(model)

	@assert is_solved_and_feasible(model)

	solution_summary(model)

	return objective_value(model)
end
