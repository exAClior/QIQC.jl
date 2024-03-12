using Test
using QIQC
using Yao
using JuMP

@testset "entanglement" begin
	obj_val, As = chsh_op_povm(nothing)
	obj_val2, Bs = chsh_op_povm(As)
	@test isapprox(obj_val2, 2 * sqrt(2), atol = 1e-4)
end
