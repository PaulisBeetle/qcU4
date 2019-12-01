using LinearAlgebra
using Yao
using Test

# the first qubit is ancilla
nugate(Σ; ϵ::Float64) = matblock([cos(ϵ*Σ) -sin(ϵ*Σ); sin(ϵ*Σ) cos(ϵ*Σ)])
function nugate(Σ::Diagonal; ϵ::Float64)
    c = Diagonal(cos.(ϵ .* Σ.diag))
    s = Diagonal(sin.(ϵ .* Σ.diag))
    matblock([c -s; s c])
end

@testset "rot unitary" begin
    blk = nugate([1.0 0; 0 0])
    @test isunitary(blk)
    @test nqubits(blk) == 2
    @show mat(blk)
end

function reflect_circuit(gen::AbstractBlock{N}) where N
    reflect0 = control(N, -collect(1:N-1), N=>-Z)
    chain(gen', reflect0, gen)
end

function grover_step!(reg::AbstractRegister, U::AbstractBlock)
    nbit = nqubits(U)
    oracle = put(nbit, nbit=>Z)
    reg |> oracle |> reflect_circuit(U)
end

function grover_nonunitary(nbit::Int, pair::Pair; k, ϵ)
    locs = pair.first
    Σ = pair.second
    U = chain(repeat(nbit, H, 1:nbit-1), put(nbit, (locs..., nbit)=>nugate(Σ; ϵ=π/4)))
    reg = zero_state(nbit)
    reg |> U
    for i=1:k
        grover_step!(reg, U)
    end
    return reg
end

@testset "non-unitary" begin
    Σ = Diagonal(ComplexF64[1.0, 0])
    reg = grover_nonunitary(2, (1,)=>Σ; k=1, ϵ = π/4)
    @test all(measure(reg, 2, nshots=10) .== 1)
end

# k =1
# epsilon => π/4
# Σ = [1, 0]
function tstar(Σ::AbstractBlock{N}, ϵ) where N
    expect(cos(ϵ*Σ)^2, uniform_state(N))
end

function getϵ(k::Int; ϵ = 0.1)
    return 0.1
    #loss = (tstar(Σ, ϵ, U) - cos(π/(4k+2))^2)
end

@testset "tstar" begin
    k = 1
    @test ϵ ≈ π/4
end
