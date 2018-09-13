module TheoryOfCISS

using LinearAlgebra
using Molecules
using Molecules: Hydrogen, Carbon, Lead, Å, angularmomentum, countorbitals, nB, nF
using Unitful: Energy, eV, meV, ħ

import Molecules: bondenergy, onsiteenergy

export makemol, simulate

"""
    makemol(N, ϕ, l)

"""
function makemol(N, ϕ=π/4, l=1Å)
    δz = √3 * l * sin(ϕ) / 2
    θ = π - 2atan(√3 * cos(ϕ))
    𝐫0 = [√(l^2 - δz^2) / 2sin(θ/2), 0Å, 0Å]

    mol = Molecule{Float64}()

    axes = Molecules.makeaxes(Float64, 0.0, ϕ)
    A′ = Hydrogen(𝐫0, orbital"1s", axes=axes)
    push!(mol, A′)
    for i in 1:N
        rotation_matrix = LinearAlgebra.Givens(1, 2, cos(i*θ), sin(i*θ))
        shift_vector = [0Å, 0Å, i*δz]
        axes = Molecules.makeaxes(Float64, i*θ, ϕ)
        𝐫1 = shift_vector + (rotation_matrix * (𝐫0/1Å))*1Å
        𝐫2 = shift_vector + (rotation_matrix * ((𝐫0 + [l, 0Å, 0Å])/1Å))*1Å

        A = Carbon(𝐫1, orbital"2s", orbital"2p_x", orbital"2p_y", orbital"2p_z", axes=axes)
        B = Hydrogen(𝐫2, orbital"1s", axes=axes)
        push!(mol, A, B)
        push!(mol, Bond(A, A′), Bond(A, B))
        A′ = A
    end
    rotation_matrix = LinearAlgebra.Givens(1, 2, cos((N+1)*θ), sin((N+1)*θ))
    shift_vector = [0Å, 0Å, (N+1)*δz]
    A = Hydrogen(shift_vector + (rotation_matrix * (𝐫0/1Å))*1Å, orbital"1s")
    push!(mol, A)
    push!(mol, Bond(A, A′))
    return mol
end

# We model the lead using a dummy atom, with a four tunable levels of
# s, p_x, p_y and p_z type.

const DummyAtom = Atom{<:Any, :D, 1, 0}

onsiteenergy(::DummyAtom, ::Orbital"2s") = 0.0eV
onsiteenergy(::DummyAtom, ::Orbital"2p_{*}") = 0.0eV

bondenergy(::Val{:σ}, δr, ::DummyAtom, ::Carbon, ::Orbital"2s", ::Orbital"2s") = -1.0eV
bondenergy(::Val{:σ}, δr, ::Carbon, ::DummyAtom, ::Orbital"2s", ::Orbital"2s") = -1.0eV
bondenergy(::Val{:σ}, δr, ::DummyAtom, ::Carbon, ::Orbital"2s", ::Orbital"2p_{*}") = -1.0eV
bondenergy(::Val{:σ}, δr, ::Carbon, ::DummyAtom, ::Orbital"2p_{*}", ::Orbital"2s") = 1.0eV
bondenergy(::Val{:σ}, δr, ::Carbon, ::DummyAtom, ::Orbital"2s", ::Orbital"2p_{*}") = 1.0eV
bondenergy(::Val{:σ}, δr, ::DummyAtom, ::Carbon, ::Orbital"2p_{*}", ::Orbital"2s") = -1.0eV
bondenergy(::Val{:σ}, δr, ::DummyAtom, ::Carbon, ::Orbital"2p_{*}", ::Orbital"2p_{*}") = 1.0eV
bondenergy(::Val{:σ}, δr, ::Carbon, ::DummyAtom, ::Orbital"2p_{*}", ::Orbital"2p_{*}") = 1.0eV
bondenergy(::Val{:π}, δr, ::DummyAtom, ::Carbon, ::Orbital"2p_{*}", ::Orbital"2p_{*}") = -0.5eV
bondenergy(::Val{:π}, δr, ::Carbon, ::DummyAtom, ::Orbital"2p_{*}", ::Orbital"2p_{*}") = -0.5eV

function simulate(N; η = 1e-15eV, ω0 = 32meV, m = 1.0, λ = 4meV, axis = :z)
    # Construct molecule and leads
    mol = makemol(N)
    hydrogen_atoms = collect(filter(x -> x isa Hydrogen, copy(mol.atoms)))
    carbon_atoms = collect(filter(x -> x isa Carbon, copy(mol.atoms)))

    atom_left = DummyAtom(0.0Å, 0.0Å, carbon_atoms[1].position[3] - 1.0Å, orbital"2s", orbital"2p_x", orbital"2p_y", orbital"2p_z")
    atom_right = DummyAtom(0.0Å, 0.0Å, carbon_atoms[end].position[3] + 1.0Å, orbital"2s", orbital"2p_x", orbital"2p_y", orbital"2p_z")
    lead_left = Lead(Molecule([atom_left], Set(Bond[])), 1000.0eV)
    lead_right = Lead(Molecule([atom_right], Set(Bond[])), 1000.0eV)

    μ, Δ = chemicalpotential(mol)
    H0 = hamiltonian(mol)
    ΣL = selfenergy(μ, lead_left, mol, Set([Bond(atom_left, carbon_atoms[1])]))
    ΣR = selfenergy(μ, lead_right, mol, Set([Bond(atom_right, carbon_atoms[end])]))
    ΓL, ΓR = -imag(ΣL)/2, -imag(ΣR)/2

    k = rand(countorbitals(mol), countorbitals(mol))
    γL = ((1-m)*I + m*k'k).*ΓL

    V = sum(vibrationalcoupling(mol, atom, ω0 / ħ, atom.axes[:,3]*1e-6Å) for atom in hydrogen_atoms)
    L = angularmomentum(axis, mol)

    𝐈 = one(H0)
    G(E,n) = inv(@. (E + im*η)*𝐈 - H0 + im*(ΓL + ΓR) - λ*n*L/ħ)

    T0(E,n) = tr(γL*G(E,n)'ΓR*G(E,n))

    function T(E::Energy, n, α=1, β=1, γ=1)
        GE = G(E, n)
        GEω₊ = G(E + ω0, n)
        GEω₋ = G(E - ω0, n)

        X = γL*GE'
        Y = ΓR*GE
        Z₊ = GEω₊*V
        Z₋ = GEω₋*V

        # We use the identity `tr(X'Y) = dot(X, Y)`, where `dot(A, B)`
        # refers to the Frobenius inner product
        A0 = dot(X', Y)

        # Self energies
        WA = V'GE'Y*X

        A1 = dot(Z₊, WA)
        A2 = dot(Z₋, WA)

        WB = GE*X*Y*V

        B1 = dot(WB', Z₊)
        B2 = dot(WB', Z₋)

        # Vertex corrections
        WC = GE*X
        C1 = dot(Z₊, ΓR*Z₊*WC)
        C2 = dot(Z₋, ΓR*Z₋*WC)

        return A0 + nB(ω0) * (α*A1 + β*B1 + γ*C1) + (1 + nB(ω0)) * (α*A2 + β*B2 + γ*C2)
    end

    function f(E::Energy, α=1, β=1, γ=1)
        Tup = T(E, 1, α, β, γ)
        Tdown = T(E, -1, α, β, γ)
        return (Tup - Tdown) / (Tup + Tdown)
    end

    return f, T, V, μ, Δ
end

function vibration(mol, ω0, u)
    gs = [phononcoupling(mol, idx, ω0 / ħ, axes * u)
          for (idx, axes) in zip(find(a->isa(a, Hydrogen), mol.atoms)[2:end-1], localaxes(mol).axes)]
    ωs = [ω0 for _ in 1:length(gs)]
    return ωs, gs
end

end # module
