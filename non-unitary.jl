using Yao
using Yao.ConstGate
using BitBasis
using LinearAlgebra
using Plots

struct unidecomp
    coeff::ComplexF64
    unit::Array{ComplexF64,2}
    unidecomp(coeff,unit) = new(coeff,unit)
end

function run_()
    f = 1.0
    for i=1:10
        reg |> g |> p0
        f *= normalize_factor(reg)
        reg |> normalize!
    end
    return f
end

#=
non-unitary probabilitistic gate
N:non-unitary gate, ϵ:parameter
return fidelity and probability
=#

function nonunitfidprob(N,ϵ)
    m = exp([zeros(4,4) ϵ*N; -ϵ*N zeros(4,4)])
    reg2 = rand_state(2)
    reg = join(zero_state(1),reg2)
    #println(statevec(reg))
    p31 = put(3,3=>ConstGate.P0)
    regt = copy(reg)
    reg |> matblock(m)
    p0 = probs(reg)[1]
    reg |> matblock(m) |> p31 |> normalize!
    #println(statevec(reg))
    reg2 |> matblock(N) |> normalize!
    return fidelity(reg,join(zero_state(1),reg2)), p0
end

function unitarydecomposition(N) #N = n1.coeff*n1.unit+n2.coeff*n2.unit
    nfactor = opnorm(N)
    N = N / nfactor
    F = svd(N)
    D = Diagonal(F.S) + im*sqrt.(I-Diagonal((F.S).^2))
    return unidecomp(0.5*nfactor,F.U*D*F.Vt),unidecomp(0.5*nfactor,F.U*conj(D)*F.Vt)
end

function lcu(ntot,top,mid,bot,N)
    ua,ub = unitarydecomposition(N)
    chain(
    ntot,
    put(ntot, (top=>H)),
    put(ntot, (top=>X)),
    control(ntot,top,(mid,bot)=>matblock(ua.unit)),
    put(ntot, (top=>X)),
    control(ntot,top,(mid,bot)=>matblock(ub.unit)),
    put(ntot, (top=>H)),
    put(ntot,(top=>ConstGate.P0))
    )
end
