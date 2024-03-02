using Yao
using JuMP
using SCS: SCS
using LinearAlgebra

"""
Given two mixed states ρ and σ, find the measurement that maximizes the
measurement result difference between the two states. They should be equal to
the trace distance between the two states.
"""
function sdp_measurement(
    ρ::AbstractMatrix, σ::AbstractMatrix; optimizer=SCS.Optimizer, silent=true
)
    model = Model(optimizer)
    silent && set_silent(model)

    M = @variable(model, [1:Base.size(ρ, 1), 1:Base.size(ρ, 2)] in HermitianPSDCone())

    @constraint(model, LinearAlgebra.I - M in HermitianPSDCone())

    @objective(model, Max, real(tr(M * (ρ - σ))))

    optimize!(model)

    @assert is_solved_and_feasible(model)

    solution_summary(model)

    return objective_value(model)
end

function eval_measure(ρ::DensityMatrix, θs::AbstractVector) end

"""
Compute the trace distance between two mixed states variationally. 

D(ρ,σ) = max |Tr(M(ρ - σ))|, where the maximum is taken over all POVM operators M with 0 ≤ M ≤ 1.
"""
function variational_distance(ρ::DensityMatrix, σ::DensityMatrix, stepsize::AbstractFloat)
    loss = 1.0
    θs = rand(3 * Yao.log2i(Base.size(ρ, 1)))
    while loss >= 1e-6
        loss = 1.0 - trace()
        θs = θs .- stepsize .* grad
    end
    return 1.0 - loss
end
