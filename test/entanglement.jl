using Test
using QIQC
using Yao

@testset "entanglement" begin
	@time chsh_op_povm()
end
