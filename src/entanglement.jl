using JuMP
using Yao
using LinearAlgebra
using SCS

function cabxy(a, b, x, y)
	return (-1)^(a + b + x * y)
end

# function cabxy(a, b, x, y)
# 	return (-1)^(a * x + b * y)
# end

# function cabxy(a, b, x, y)
# 	return (0.5)^(a + b + x * y) + (-0.5)^(x + y + a * b)
# end

# function cabxy(a, b, x, y)
# 	return exp(a)^b^x^y 
# end

function chsh_op_povm(target::String, As, Bs, ρ, n, optimizer = SCS.Optimizer, silent = true)
	model = Model(optimizer)
	silent && set_silent(model)
	if isnothing(ρ)
		if target == "ρ"
			ρ = @variable(model, [1:(n^2), 1:(n^2)] in HermitianPSDCone())
			@constraint(model, real(tr(ρ)) == 1.0)
			# enforce PPT critesion
			@constraint(model, partial_T(ρ) >=  0 , PSDCone())
		end
	end

	v0 = rand_state(1; nlevel = n)
	v1 = rand_state(1; nlevel = n)

	if isnothing(Bs)
		Bs = Matrix{ComplexF64}[]
		B00 = v0.state * v0.state'
		B10 = I - B00
		B01 = v1.state * v1.state'
		B11 = I - B01
		Bs = [B00, B01, B10, B11]
	end

	if isnothing(As)
		target == "ρ" && error("whoa")
		As = [@variable(model, [1:n, 1:n] in HermitianPSDCone()) for _ in 1:4]
		@constraint(model, As[1] + As[3] == LinearAlgebra.I)
		@constraint(model, As[2] + As[4] == LinearAlgebra.I)
	end

	@objective(model, Max, real(sum([cabxy(a, b, x, y) * tr(kron(As[(a) * 2 + x + 1], Bs[(b) * 2 + y + 1]) * ρ) for a in 0:1, b in 0:1, x in 0:1, y in 0:1])))

	optimize!(model)

	@assert is_solved_and_feasible(model)

	solution_summary(model)

	if target == "ρ"
		return objective_value(model), value.(ρ)
	else
		return objective_value(model), [value.(As[ii]) for ii in 1:4]
	end
end
