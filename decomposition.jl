#=
Decomposition of 2-qubit unitary operator
=#
using LinearAlgebra;
using Yao
using QuantumInformation;

const MagicBasis = 1/sqrt(2)*[1 -im 0 0; 0 0 1 -im; 0 0 -1 -im; 1 im 0 0];

"""
    uaub(tpsi::AbstractMatrix) -> UA, UB, ξ

Decompose a unitary matrix `tpsi` to `kron(UA, UB)*?`.
"""
function uaub(tpsi::AbstractMatrix)
    γ = [get_angle(tpsi[:,i]) for i in 1:4];
    mtpsi = MagicBasis * tpsi;
    barpsi = tpsi .* exp.(-im.*γ)';

    ef_compbasis = 1/sqrt(2)*(barpsi[:,1] - 1*im*barpsi[:,2]);
    epfp_compbasis = 1/sqrt(2)*(barpsi[:,1] + 1*im*barpsi[:,2]);
    #convert ef and epfp into computational basis
    ef = MagicBasis*ef_compbasis;
    epfp = MagicBasis*epfp_compbasis;

    e, f = decompose_kron(ef)
    ep, fp = decompose_kron(epfp)

    barpsi_compbasis = MagicBasis * barpsi;
    efbasis = [kron(e,f) kron(e,fp) kron(ep,f) kron(ep,fp)];
    coeff = transpose(efbasis'*barpsi_compbasis)

    ξ, phase = modify_phase(γ, coeff)

    UA = kron([1;0], e') + kron([0;1], ep') * phase;
    UB = kron([1;0], f') + kron([0;1], fp') * inv(phase);
    return UA, UB, ξ
end

"""
    modify_phase(γ::Vector, coeff::AbstractMatrix) -> ξ, phase

...
"""
function modify_phase(γ::Vector{T}, coeff::AbstractMatrix) where T
    ξ = zeros(T, 4); ξ[1] = -γ[1]; ξ[2] = -γ[2] + π;
    index = 0;
    if 2*coeff[3,2]*coeff[3,3] ≈ -1
        phase = sqrt(2) * coeff[3,2];
        phase2 = im * sqrt(2) * coeff[4,2];
        if phase ≈ phase2
            ξ[3] = -γ[3];
            ξ[4] = -γ[4];
        else
            ξ[3] = -γ[3];
            ξ[4] = π - γ[4];
        end
    else
        phase = sqrt(2) * coeff[3,2] *im;
        phase2 = sqrt(2) * coeff[4,2];
        if phase ≈ phase2
            ξ[3] = -γ[3] + π/2;
            ξ[4] = -γ[4] - π/2;
        else
            ξ[3] = -γ[3] + π/2;
            ξ[4] = -γ[4] + π/2;
        end
        index = 1;
    end
    return ξ, phase
end

"""
"""
function get_angle(psi::Vector; tol=1e-4)
    norms = abs2.(psi)
    res = angle(psi[argmax(norms)])
    if !all(p->abs2(p) < tol || isapprox(mod(angle(p)-res+tol, π), 0; atol=2tol), psi)
        error("phase error, got $(angle.(psi))")
    end
    return res
end

"""
"""
function decomposition(u1::AbstractMatrix)
    u2 = MagicBasis' * u1 * MagicBasis;

    # ...
    uu = transpose(u2)*u2;
    eigensystem = eigen(uu);

    ϵ = 1/2*map(angle,eigensystem.values);
    psi = eigensystem.vectors;
    VA, VB, xi = uaub(psi);
    tildapsi = (u2 * psi) * Diagonal(exp.(-im.*ϵ))
    UAD,UBD, s = uaub(tildapsi);
    λ = s - ϵ - xi;
    Ud = MagicBasis * Diagonal(exp.(-im .* λ)) * MagicBasis'
    return UAD', UBD', Ud, VA, VB
end

"""
"""
function decompose_kron(ef::AbstractVector{T}) where T
    A = kron(ef, ef');
    crossff = ptrace(A,[2,2],[1]);
    f = eigen(crossff).values[1] ≈ 1 ? eigen(crossff).vectors[:,1] : eigen(crossff).vectors[:,2];
    e = isapprox(f[1], 0; atol=1e-2) ? [ef[2]/f[2], ef[4]/f[2]] : [ef[1]/f[1], ef[3]/f[1]]
    return e, f
end
