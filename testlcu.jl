using Test
include("non-unitary.jl")


@testset "non-unitary.jl" begin
    N = Diagonal(rand(4))
    circ = lcu(3,1,2,3,N)
    reg2 = rand_state(2)
    regt = copy(reg2)
    reg  = join(reg2,zero_state(1))
    @show statevec(reg)
    reg |> circ |> normalize!
    regt |> matblock(N) |> normalize!
    @test fidelity(reg,join(regt,zero_state(1)))[1] ≈ 1.0
end
