using Yao
using LinearAlgebra
using Roots
using LaTeXStrings
using Plots


function reflect_circuit(gen::AbstractBlock{N}) where N
    reflect0 = control(N, -collect(1:N-1),N=>-Z)
    chain(gen',reflect0,gen)
end

ΣW(i,W,τ) = W[i] > 0 ?  Diagonal(Complex.([1.,exp(-2*W[i]*τ),exp(-2*W[i]*τ),1.])) : Diagonal(Complex.([exp(2*W[i]*τ),1.,1.,exp(2*W[i]*τ)]))
UΣ(Σ,ϵ) = matblock([cos(ϵ*Σ) -sin(ϵ*Σ); sin(ϵ*Σ) cos(ϵ*Σ)])

function target_state(Σ1,Σ2,ϵ)
    nbit = 5
    zero_state(nbit) |> repeat(nbit,H,(1,2,4)) |> repeat(nbit,X,(3,5)) |> put(nbit,(1,2)=>matblock(sin(ϵ*Σ1))) |> put(nbit,(2,4)=>matblock(sin(ϵ*Σ2))) |> normalize!
end

function deps(Σ1,Σ2,reg,ϵ)
    c = chain(3,put((1,2)=>matblock(sin(ϵ*Σ1)^2)), put((2,3)=>matblock(sin(ϵ*Σ2)^2)))
    return abs(expect(c,reg))
end

function probsuc(Σ1,Σ2,reg,ϵ,k)
    d = Complex(deps(Σ1,Σ2,reg,ϵ))
    return abs((sqrt(d)+sqrt(d-1))^(2k+1)+(sqrt(d)-sqrt(d-1))^(2k+1))/4
end

function probsuc(d,k)
    d = Complex(d)
    return abs((sqrt(d)+sqrt(d-1))^(2k+1)+(sqrt(d)-sqrt(d-1))^(2k+1))/4
end

function doublegrover(Σ1,Σ2,d)
    ϵ = minimum(find_zeros(x->deps(Σ1,Σ2,uniform_state(3),x)-d,0,π))
    Ui = repeat(nbit,H,(1,2,4))
    oracle = control(5,(3,),5=>Z)
    gen = chain(Ui,put((1,2,3)=>UΣ(Σ1,ϵ)),put((2,4,5)=>UΣ(Σ2,ϵ)))
    reg = zero_state(5) |> gen
    prob = []
    t_state = target_state(Σ1,Σ2,ϵ)
    for i = 1:10
        reg |> oracle |> reflect_circuit(gen)
        p = abs(reg'*t_state)
        push!(prob,p)
    end
    prob
end

#1,2 and 4 are physical qubits, 3 and 5 are ancillary qubits.
#Σ1 acts on the (1,2). Σ2 acts on (2,4).

Σ1 = Diagonal(Complex.([1.,exp(-5),exp(-5),1.]))
Σ2 = Diagonal(Complex.([exp(-5),1.,1.,exp(-5)]))
#nbit = 5
#Ui = repeat(nbit,H,(1,2,4))
d = 1/4
#ϵ = minimum(find_zeros(x->deps(Σ1,Σ2,uniform_state(3),x)-d,0,π))
#gen = chain(Ui,put((1,2,3)=>UΣ(Σ1,ϵ)),put((2,4,5)=>UΣ(Σ2,ϵ)))
#reg = zero_state(nbit) |> gen |> control(5,(3,),5=>Z) |> reflect_circuit(gen)
#measure(reg,(3,5),nshots=10)
plot(doublegrover(Σ1,Σ2,d))
