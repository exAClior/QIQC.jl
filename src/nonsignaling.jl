using JuMP
using SCS

function nonsignaling_game(optimizer = SCS.Optimizer, silent = true, signaling = false)
	# optimizer = SCS.Optimizer
	model = Model(optimizer)
	# silent = true
	silent && set_silent(model)

	pabxy = @variable(model, [1:2, 1:2, 1:2, 1:2], lower_bound = 0.0, upper_bound = 1.0)
	pax = @variable(model, [1:2, 1:2], lower_bound = 0.0, upper_bound = 1.0)
	pby = @variable(model, [1:2, 1:2], lower_bound = 0.0, upper_bound = 1.0)

	for given_val in 1:2
		@constraint(model, sum(pax[:, given_val]) == 1.0)
		@constraint(model, sum(pby[:, given_val]) == 1.0)
	end

	for y in 1:2
		for x in 1:2
			@constraint(model, sum(pabxy[:, :, x, y]) == 1.0)
			@constraint(model, sum(pabxy[:, :, x, y]) == 1.0)
			if !signaling
				for a in 1:2
					@constraint(model, sum(pabxy[a, :, x, y]) == pax[a, x])
				end
				for b in 1:2
					@constraint(model, sum(pabxy[:, b, x, y]) == pby[b, y])
				end
			end
		end
	end

	@objective(model, Max, sum([cabxy(a, b, x, y) * pabxy[a + 1, b + 1, x + 1, y + 1] for a in 0:1, b in 0:1, x in 0:1, y in 0:1]))

	optimize!(model)

	@assert is_solved_and_feasible(model)

	println(solution_summary(model))

    return value.(pabxy)
end
