using JuMP
using Yao
using LinearAlgebra
using SCS

# function cabxy(a, b, x, y)
# 	return (-1)^(a + b + x * y)
# end

# function cabxy(a, b, x, y)
# 	return (-1)^(a * x + b * y)
# end

# function cabxy(a, b, x, y)
# 	return (0.5)^(a + b + x * y) + (-0.5)^(x + y + a * b)
# end

# function cabxy(a, b, x, y)
# 	return exp(a)^b^x^y 
# end
# using Random
# rand(16)
# res = rand(16)
# res

function cabxy(a, b, x, y)
	table = [0.7954447180190837,
		0.4068562539710824,
		0.8863966495001979,
		0.4730542413777712,
		0.5406657962237318,
		0.42004910733789425,
		0.24600566923294576,
		0.34251279156147585,
		0.6700891867459674,
		0.7300059855279646,
		0.7347221543014663,
		0.6534003223880149,
		0.5516352706310135,
		0.5992477843905545,
		0.35285383974630824,
		0.1663026323377389]
	table = reshape(table, (2, 2, 2, 2))
	return table[a + 1, b + 1, x + 1, y + 1]
end

function chsh_op_povm(target::String, As, Bs, ρ, n, meas_n, optimizer = SCS.Optimizer, silent = true)
	model = Model(optimizer)
	silent && set_silent(model)
	if isnothing(ρ)
		if target == "ρ"
			ρ = @variable(model, [1:(n^2), 1:(n^2)] in HermitianPSDCone())
			@constraint(model, real(tr(ρ)) == 1.0)
		end
		# else
		# 	ρ = density_matrix(ghz_state(2)).state
		# end
	end
	v0 = rand_state(1; nlevel = n)
	v1 = rand_state(1; nlevel = n)

	if isnothing(Bs)
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
