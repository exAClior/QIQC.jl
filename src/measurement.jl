using Yao

function eval_measure(ρ::DensityMatrix, θs::AbstractVector)
end

"""
Compute the trace distance between two mixed states variationally. 

D(ρ,σ) = max |Tr(M(ρ - σ))|, where the maximum is taken over all POVM operators M with 0 ≤ M ≤ 1.
"""
function variational_distance(ρ::DensityMatrix, σ::DensityMatrix, stepsize::AbstractFloat)
	loss = 1.0
	θs = rand(3 * Yao.log2i(size(ρ, 1)))
	while loss >= 1e-6
		loss = 1.0 - trace()
		θs = θs .- stepsize .* grad
	end
	return 1.0 - loss
end
