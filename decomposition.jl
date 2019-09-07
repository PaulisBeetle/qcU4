
#=
Decomposition of 2-qubit unitary operator
=#
using LinearAlgebra;
using QuantumInformation;

MagicBasis = 1/sqrt(2)*[1 -im 0 0; 0 0 1 -im; 0 0 -1 -im; 1 im 0 0];

function generation_u4()
    CNOT = [1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0];
    SWAP = [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1];
    iSWAP = [1 0 0 0; 0 0 im 0; 0 im 0 0; 0 0 0 1];
    halfSWAP = [1 0 0 0;0 1/2*(1+im) 1/2*(1-im) 0;0 1/2*(1-im) 1/2*(1+im) 0;0 0 0 1];
    X = [0 1; 1 0];
    Y = [0 -im; im  0];
    Z = [1 0; 0 -1];
    c1 = rand(3);
    c2 = rand(3);
    su2_1 = exp(-im * c1[1] * X - im * c1[2] * Y - im * c1[3] * Z);
    su2_2 = exp(-im * c2[2] * X - im * c2[2] * Y - im * c2[3] * Z);
    println("su2_1: $su2_1");
    println("su2_2: $su2_2");
    su2su2 = kron(su2_1,su2_2);
    println("su2su2: $su2su2");
    #return su2su2
    return halfSWAP;
end

function uaub(tpsi)
    normpsi = norm.(tpsi);
    γ = [0.0; 0.0; 0.0; 0.0];
    γ = [angle(tpsi[argmax(normpsi[:,i]),i]) for i in 1:4];
    #barpsi = barf(tpsi);
    mtpsi = MagicBasis * tpsi;
    barpsi = copy(tpsi);
    for i in 1:4
        barpsi[:,i]=tpsi[:,i]*exp(-im*γ[i]);
    end
    ef_compbasis = 1/sqrt(2)*(barpsi[:,1] - 1*im*barpsi[:,2]);
    epfp_compbasis = 1/sqrt(2)*(barpsi[:,1] + 1*im*barpsi[:,2]);
    #convert ef and epfp into computational basis
    ef = MagicBasis*ef_compbasis;
    epfp = MagicBasis*epfp_compbasis;
    A = kron(ef,adjoint(ef));
    crossff = ptrace(A,[2,2],[1]);
    f = eigen(crossff).values[1] ≈ 1 ? eigen(crossff).vectors[:,1] : eigen(crossff).vectors[:,2];
    e = [0.0 + 0.0*im ; 0.0*im];
    e[1] = (1+f[1])≈1 == true ? ef[2]/f[2] : ef[1]/f[1];
    e[2] = (1+f[1])≈1 == true ? ef[4]/f[2] : ef[3]/f[1];
    B = kron(epfp,adjoint(epfp));
    crossfpfp = ptrace(B,[2,2],[1]);
    fp = eigen(crossfpfp).values[1] ≈ 1 ? eigen(crossfpfp).vectors[:,1] : eigen(crossfpfp).vectors[:,2];
    ep = [0.0 + 0.0*im ; 0.0 + 0.0*im];
    ep[1] = (1+fp[1])≈1 == true ? epfp[2]/fp[2] : epfp[1]/fp[1];
    ep[2] = (1+fp[1])≈1 == true ? epfp[4]/fp[2] : epfp[3]/fp[1];
    #tpsi1 = 1/sqrt(2)*(kron(e,f) + kron(ep,fp));
    #tpsi2 = im/sqrt(2)*(kron(e,f) - kron(ep,fp));
    #barpsi1 = MagicBasis * barpsi[:,1];
    #barpsi2 = MagicBasis * barpsi[:,2];
    barpsi_compbasis = MagicBasis * barpsi;
    efbasis = [kron(e,f) kron(e,fp) kron(ep,f) kron(ep,fp)];
    coeff = barpsi;
    for i in 1:4
        for j in 1:4
            coeff[i,j] = adjoint(efbasis[:,j])*barpsi_compbasis[:,i];
        end
    end
    println("coeff:\n $coeff");
    ξ = zeros(4); ξ[1] = -γ[1]; ξ[2] = -γ[2] + π;
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
    println(phase);
    UA = kron([1;0],adjoint(e)) + kron([0;1],adjoint(ep)) * phase;
    UB = kron([1;0],adjoint(f)) + kron([0;1],adjoint(fp)) * inv(phase);
    Dξ = [exp(im*ξ[1]) 0 0 0; 0 exp(im*ξ[2]) 0 0;0 0 exp(im*ξ[3]) 0; 0 0 0 exp(im*ξ[4])];
    println(kron(UA,UB) * mtpsi * Dξ - MagicBasis)
    println("index:  $index")
    return UA, UB, ξ
end

function decomposition(u1)
    u2 = adjoint(MagicBasis) * u1 * MagicBasis;
    println("u2: $u2");
    uu = transpose(u2)*u2;
    println("uu: $uu");
    eigensystem = eigen(uu);
    ϵ = 1/2*map(angle,eigensystem.values);
    psi = eigensystem.vectors;
    println("psi: $psi");
    whypsichange = copy(psi);
    #return uaub(psi)
    VA, VB, xi = uaub(psi);
    #D = [exp(-im*ϵ[1]) 0 0 0 ; 0 exp(-im*ϵ[2]) 0 0 ; 0 0 exp(-im*ϵ[3]) 0 ; 0 0 0 exp(-im*ϵ[4])];
    tildapsi = copy(psi);
    for i in 1:4
        tildapsi[:,i] = exp(-im*ϵ[i]) * u2 * whypsichange[:,i];
    end
    println("tildapsi: $tildapsi");
    whytildachange = copy(tildapsi);
    UAD,UBD, s = uaub(tildapsi);
    UA = adjoint(UAD);
    UB = adjoint(UBD);
    λ = s - ϵ - xi;
    println("s: $s");
    println("ϵ: $ϵ");
    println("xi: $xi");
    println("λ: $λ");
    println("sum of λ: $(sum(λ))");
    Ud = zeros(Complex,4,4);
    for i in 1:4
        Ud = Ud + exp(-im * λ[i]) * kron(MagicBasis[:,i],adjoint(MagicBasis[:,i]));
    end
    U = kron(UA,UB) * Ud * kron(VA,VB);
    println("det(input): $(det(u1))");
    println("det(U): $(det(U))");
    println("det(kron(UA,UB)): $(det(kron(UA,UB)))");
    println("det(kron(Ud)): $(det(Ud))");
    println("det(kron(VA,VB)): $(det(kron(VA,VB)))");
    println(U-u1);
    return
end

decomposition(generation_u4())
