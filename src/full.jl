struct Simulation{T}
    N::Integer
    mol::Molecule
    μ::T
    Δ::T
    H0::AbstractMatrix
    ΣL::Function
    ΣR::Function
    ΣA::Function
    ΓL::Function
    ΓR::Function
    ΓA::Function
    L::AbstractMatrix
    Vs::AbstractVector
end

function Simulation(N::Integer; η = 1e-15, axis = :z, ϕ = π/4, l=2, handedness=1, α=1., Γ=1e1)
    # Construct molecule and leads
    mol = makemol(N, ϕ, l, handedness=handedness, vec=dvec)
    skt = loaddir(joinpath(PARAMDIR, "ob2-1-1", "base"))
    hydrogen_atoms = collect(filter(x -> x isa Hydrogen, mol.atoms))
    carbon_atoms = collect(filter(x -> x isa Carbon, mol.atoms))

    atom_left = Carbon(
        [carbon_atoms[1].position[1] - 0, carbon_atoms[1].position[2] + 0, carbon_atoms[1].position[3] - l],
        orbital"2s", orbital"2p_x", orbital"2p_y", orbital"2p_z")
    atom_right = Carbon(
        [carbon_atoms[end].position[1] - 0, carbon_atoms[end].position[2] + 0, carbon_atoms[end].position[3] + l],
        orbital"2s", orbital"2p_x", orbital"2p_y", orbital"2p_z")
    atom_extra = Carbon(
        [carbon_atoms[3].position[1] - 0, carbon_atoms[3].position[2] - 0, carbon_atoms[3].position[3] - l],
        orbital"2s", orbital"2p_x", orbital"2p_y", orbital"2p_z")
    lead_left = Lead(Molecule([atom_left], Set(Bond[])), Γ)
    lead_right = Lead(Molecule([atom_right], Set(Bond[])), Γ)
    lead_extra = Lead(Molecule([atom_extra], Set(Bond[])), Γ)

    # Calculate H and all derivatives w.r.t. positions of atoms in one
    # go, using automatic forward differentiation
    T = Dual{Nothing,Float64,3*(2N + 2)}
    HdH = hamiltonian(T, skt, mol)
    # Extract the unperturbed hamiltonian
    H0 = (x->x.value).(HdH)
    μ, Δ = real.(chemicalpotential(H0, mol))
    # Extract the perturbation components and rotate them into the
    # local axes
    V(i) = (x->x.partials.values[i]).(HdH)
    Vs = reduce(vcat, [(1/sqrt(2mass(skt, atom))*α*atom.axes*[V(3i - 2), V(3i - 1), V(3i)])[2:2] for (i,atom) in enumerate(mol.atoms) if atom isa Carbon])
    # Find the self energies
    ΣL = selfenergy(T, skt, lead_left, mol, Set([Bond(atom_left, carbon_atoms[1])]), f=(x->x.value))
    ΣR = selfenergy(T, skt, lead_right, mol, Set([Bond(atom_right, carbon_atoms[end])]), f=(x->x.value))
    ΣA = selfenergy(T, skt, lead_extra, mol, Set([Bond(atom_extra, carbon_atoms[1])]), f=(x->x.value))

    ΓL(E) = (x = ΣL(E); im*(x - x'))
    ΓR(E) = (x = ΣR(E); im*(x - x'))
    ΓA(E) = (x = ΣA(E); im*(x - x'))

    L = angularmomentum(axis, mol)

    return Simulation(N, mol, μ, Δ, H0, ΣL, ΣR, ΣA, ΓL, ΓR, ΓA, L, Vs)
end

function G(E; sim::Simulation, n::Real, η=1e-16)
    A = (E + im*η)*I - Matrix(sim.H0) - (sim.ΣL(E) + sim.ΣR(E)) - sim.λ*n*sim.L
    # A = (E + im*η)*I - Matrix(sim.H0) + im*(sim.ΓL(E) + sim.ΓR(E))/2 - (sim.λ)*n*sim.L
    return inv(A)
end

function G(E, Σ::AbstractArray; sim::Simulation, n::Real, η=1e-16)
    # A = (E + im*η)*I - Matrix(sim.H0) - (sim.ΣL(E) + sim.ΣR(E) + Σ) - sim.λ*n*sim.L
    A = (E + im*η)*I - Matrix(sim.H0) + im*(sim.ΓL(E) + sim.ΓR(E) + im*Σ*2)/2 - (sim.λ)*n*sim.L
    return inv(A)
end

##################################################
# 
##################################################

function Δ_as_function_of_length(N, ϕ; l=2, handedness=1)
    mol = makemol(N, ϕ, l, handedness=handedness)
    skt = loaddir(joinpath(PARAMDIR, "ob2-1-1", "base"))
    H0 = hamiltonian(Float64, skt, mol)
    μ, Δ = real.(chemicalpotential(H0, mol))
    return Δ
end

##################################################
# Normal NEGF transport calculation to verify model
##################################################

function theorem_T(E, n, m=0; sim::Simulation, ω0=0.05/27, T=9.5e-4, N=nB(ω0,T), λ=0.14e-3, N_iterations=5)
    ΓL = sim.ΓL(E)
    ΓR = sim.ΓR(E)
    ΓA = sim.ΓA(E)
    γL = (1 - m)*ΓL + m*ΓA
    H0 = Matrix(sim.H0)
    f = 1 / ω0

    G = inv(E*I - H0 + im*(ΓL + ΓR)/2 - λ*n*sim.L)

    ΓR′ = copy(ΓR)
    ΓR′′ = copy(ΓR)
    for _ in 1:N_iterations
        ΓR′′ .= ΓR
        for V in sim.Vs
            ΓR′′ .+= f*(2N + 1)*V*G'ΓR′*G*V
        end
        ΓR′ .= ΓR′′
        G .= inv(E*I - H0 + im*(ΓL + ΓR′)/2 - λ*n*sim.L)
    end

    return real(dot(G*γL*G', ΓR′))
end

function theorem_f(E, m=0; sim::Simulation, T=9.5e-4, ω0=0.05/27, λ=0.14e-3, N_iterations=5)
    Tup = theorem_T(E, 1, m, sim=sim, ω0=ω0, T=T, λ=λ, N_iterations=N_iterations)
    Tdown = theorem_T(E, -1, m, sim=sim, ω0=ω0, T=T, λ=λ, N_iterations=N_iterations)
    return (Tup - Tdown), (Tup + Tdown)
end

function theorem_scan(N, m; T=9.5e-4, ω0=0.05/27, ϕ=π/4, α=1., λ=0.14e-3, Γ=1e1, axis=:z, N_iterations=5)
    print("$N.")
    sim = Simulation(N, ϕ=ϕ, α=α, axis=axis, Γ=Γ)
    return theorem_f(sim.μ, m, sim=sim, T=T, ω0=ω0, λ=λ, N_iterations=N_iterations)
end

##################################################
# Normal NEGF transport calculation to verify model (old)
##################################################

function basic_T(E, n, m=1; sim::Simulation, ω0=0.05/27, Nrange=0:10, Ns=rand(Nrange, length(sim.Vs)))
    ΓL0 = sim.ΓL(E)
    ΓR0 = sim.ΓR(E)
    ΓR1₊ = ΓR1₋ = ΓR0
    # ΓR1₊ = sim.ΓR(E + ω0)
    # ΓR1₋ = sim.ΓR(E - ω0)
    γA = sim.ΓA(E)

    GE0 = G(E, n=n, sim=sim)
    GE1₊ = GE1₋ = GE2₊ = GE2₋ = GE0
    # GE1₊ = G(E + ω0, n=n, sim=sim)
    # GE1₋ = G(E - ω0, n=n, sim=sim)
    # GE2₊ = G(E + 2*0, n=n, sim=sim)
    # GE2₋ = G(E - 2*ω0, n=n, sim=sim)

    γL = ((1 - m)*ΓL0 + m*γA)

    f = 1 / ω0
    ΣV0 = zero(GE0)
    ΣV1₊ = zero(GE0)
    ΣV1₋ = zero(GE0)
    # Correct for the fact that for the ±1 cases, there will be one
    # phonon less/more in the system, distributed with equal
    # probability in each mode
    M = 1/length(sim.Vs)
    for (N,V) in zip(Ns,sim.Vs)
        ΣV0 += f*N*V'GE1₊*V + f*(1 + N)*V*GE1₋*V'
        ΣV1₊ += f*(N - M)*V'GE2₊*V + f*(1 + N - M)*V*GE0*V'
        ΣV1₋ += f*(N + M)*V'GE0*V + f*(1 + N + M)*V*GE2₋*V'
    end

    GE0full = G(E, ΣV0, n=n, sim=sim)
    GE1₊full = G(E, ΣV1₊, n=n, sim=sim)
    GE1₋full = G(E, ΣV1₋, n=n, sim=sim)
    # GE1₊full = G(E + ω0, ΣV1₊, n=n, sim=sim)
    # GE1₋full = G(E - ω0, ΣV1₋, n=n, sim=sim)

    ΣVv = zero(GE0)
    ΣVvfull = zero(GE0)
    for (N,V) in zip(Ns,sim.Vs)
        ΣVv += f*N*V'GE1₊'ΓR1₊*GE1₊*V + f*(1 + M)*V*GE1₋'ΓR1₋*GE1₋*V'
        ΣVvfull += f*N*V'GE1₊full'ΓR1₊*GE1₊full*V + f*(1 + N)*V*GE1₋full'ΓR1₋*GE1₋full*V'
    end

    x = GE0*γL*GE0'
    X = GE0full*γL*GE0full'
    # We use the identity `tr(X'Y) = dot(X, Y)`, where `dot(A, B)`
    # refers to the Frobenius inner product
    a = dot(X, ΓR0 + ΣVv) + dot(x, ΣVvfull) - dot(x, ΣVv)
    return real(a)
end

function basic_f(E, m=1; sim::Simulation, ω0=0.05/27, Nrange=0:10, Ns=rand(Nrange, length(sim.Vs)))
    # Tn = basic_T(E, 0, m, Ns=Ns, sim=sim, ω0=ω0)
    # return Tn
    Tup = basic_T(E, 1, m, Ns=Ns, sim=sim, ω0=ω0)
    Tdown = basic_T(E, -1, m, Ns=Ns, sim=sim, ω0=ω0)
    return (Tup - Tdown), (Tup + Tdown)
end

function basic_scan(N, ω0, m; ϕ=π/4, α=1., Nrange=1:1)
    print("$N.")
    sim = Simulation(N, ϕ=ϕ, α=α)
    return basic_f(sim.μ, m, sim=sim, ω0=ω0, Nrange=Nrange)
end

##################################################
# Proving the theorem
##################################################

function proof_T(E, n, m=1; sim::Simulation, ω0=0.05/27, N::Integer=0, Nmin=0, Nmax=Inf, λ=0.14e-3)
    ΓL0 = sim.ΓL(E - N*ω0)
    ΓR0 = sim.ΓR(E - N*ω0)
    ΓL1₊ = sim.ΓL(E - (N - 1)*ω0) * (N > Nmin)
    ΓL1₋ = sim.ΓL(E - (N + 1)*ω0) * (N < Nmax)
    ΓR1₊ = sim.ΓR(E - (N - 1)*ω0) * (N > Nmin)
    ΓR1₋ = sim.ΓR(E - (N + 1)*ω0) * (N < Nmax)
    ΓA0 = sim.ΓA(E - N*ω0)

    GE0 = inv((E - N*ω0)*I - Matrix(sim.H0) + im*(ΓL0 + ΓR0)/2 - λ*n*sim.L)
    GE1₊ = inv((E - (N - 1)*ω0)*I - Matrix(sim.H0) + im*(ΓL1₊ + ΓR1₊)/2 - λ*n*sim.L) * (N > Nmin)
    GE1₋ = inv((E - (N + 1)*ω0)*I - Matrix(sim.H0) + im*(ΓL1₋ + ΓR1₋)/2 - λ*n*sim.L) * (N < Nmax)

    γL = ((1 - m)*ΓL0 + m*ΓA0)

    ΣV = zero(GE0)
    ΣVv = zero(GE0)

    L = length(sim.Vs)
    f = 1 / ω0
    Na = N > 0 ? sum((N - i)*(L - 1)^i for i in 0:N-1) : 0
    Nb = sum((N + 1 - i)*(L - 1)^i for i in 0:N)
    for V in sim.Vs
        ΣV += f*Na*V'GE1₊*V + f*Nb*V*GE1₋*V'
        ΣVv += f*Na*V'GE1₊'ΓR1₊*GE1₊*V + f*Nb*V*GE1₋'ΓR1₋*GE1₋*V'
    end

    X = GE0*γL*GE0'
    Y = ΣV'GE0'ΓR0
    Nc = length(sim.Vs)^N

    a = dot(X, Nc*ΓR0 + Y + Y' + ΣVv)
    return real(a)
end

function proof_f(E, m=1; sim::Simulation, ω0=0.05/27, Nmin=0, Nmax=1)
    Tup = Tdown = 0.
    Tn = 0.
    Z = sum(length(sim.Vs)^i for i in Nmin:Nmax)
    for N in Nmin:Nmax
        Tup += proof_T(E, 1, m, sim=sim, ω0=ω0, N=N, Nmin=Nmin, Nmax=Nmax) / Z
        Tdown += proof_T(E, -1, m, sim=sim, ω0=ω0, N=N, Nmin=Nmin, Nmax=Nmax) / Z
    end
    return (Tup - Tdown)/2, (Tup + Tdown)/2
end

function proof_scan(N, ω0, m, Nmax; ϕ=π/4, α=1.)
    print("$N.")
    sim = Simulation(N, ϕ=ϕ)
    return proof_f(sim.μ, m, sim=sim, ω0=ω0, Nmax=Nmax, α=α)
end

##################################################
# Specific electron energy
##################################################

function specific_T(E, n; sim::Simulation, m=0, ω0=0.05/27, T=9.5e-4)
    ΓL0 = sim.ΓL(E)
    ΓR0 = sim.ΓR(E)
    ΓR1₊ = sim.ΓR(E + ω0)
    ΓR1₋ = sim.ΓR(E - ω0)
    γA = sim.ΓA(E)

    GE0 = G(E, n=n, sim=sim)
    GE1₊ = G(E + ω0, n=n, sim=sim)
    GE1₋ = G(E - ω0, n=n, sim=sim)
    GE2₊ = G(E + 2ω0, n=n, sim=sim)
    GE2₋ = G(E - 2ω0, n=n, sim=sim)

    γL = ((1 - m)*ΓL0 + m*γA)

    ΣV0 = zero(GE0)
    ΣV1₊ = zero(GE0)
    ΣV1₋ = zero(GE0)
    for V in sim.Vs(ω0)
        ΣV0 += nB(ω0, T)*V'GE1₊*V + (1 + nB(ω0, T))*V*GE1₋*V'
        ΣV1₊ += nB(ω0, T)*V'GE2₊*V + (1 + nB(ω0, T))*V*GE0*V'
        ΣV1₋ += nB(ω0, T)*V'GE0*V + (1 + nB(ω0, T))*V*GE2₋*V'
    end

    G0full = G(E, ΣV0, n=n, sim=sim)
    G1₊full = G(E + ω0, ΣV1₊, n=n, sim=sim)
    G1₋full = G(E - ω0, ΣV1₋, n=n, sim=sim)

    ΣVv = zero(GE0)
    for V in sim.Vs(ω0)
        ΣVv += nB(ω0, T)*V'G1₊full'ΓR1₊*G1₊full*V + (1 + nB(ω0, T))*V*G1₋full'ΓR1₋*G1₋full*V'
    end

    X = G0full*γL*G0full'
    # # We use the identity `tr(X'Y) = dot(X, Y)`, where `dot(A, B)`
    # # refers to the Frobenius inner product
    a = dot(X, (ΓR0 + ΣVv))
    return real(a)
end

function specific_f(E; sim::Simulation, m=0, ω0=0.05/27, T=9.5e-4)
    Tup = specific_T(E, 1, m=m, sim=sim, ω0=ω0, T=T)
    Tdown = specific_T(E, -1, m=m, sim=sim, ω0=ω0, T=T)
    return (Tup - Tdown)/2, (Tup + Tdown)/2
end

function specific_scan(N; ω0=0.05/27, m=0, T=9.5e-4)
    print("$N.")
    sim = Simulation(N, ϕ=π/4, handedness=-1)
    return specific_f(sim.μ, m=m, sim=sim, ω0=ω0, T=T)
end
