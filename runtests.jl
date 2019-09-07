using Test
include("decomposition.jl")

@testset "decopositon" begin
    ef = [1, 0, 0, 0]
    e, f = decompose_kron(ef)
    @test e ≈ [1, 0]
    @test f ≈ [1, 0]

    u1 = rand_unitary(4)
    UA, UB, Ud, VA, VB = decomposition(u1)
    U = kron(UA,UB) * Ud * kron(VA,VB);
    println("det(input): $(det(u1))");
    println("det(U): $(det(U))");
    println("det(kron(UA,UB)): $(det(kron(UA,UB)))");
    println("det(kron(Ud)): $(det(Ud))");
    println("det(kron(VA,VB)): $(det(kron(VA,VB)))");
    @test U ≈ u1
end
