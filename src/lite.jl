function eigensort(H; by=identity)
    vals, vecs = eigen(H)
    s = sortperm(vals, by=by)
    return Eigen(vals[s], vecs[:, s])
end

function symmetrize!(A, B)
    size(A) != size(B) && error("A and B have different sizes")
    for (i, j) in zip(eachindex(A), eachindex(B))
        a, b = A[i], B[j]
        A[i] = (a + b) / 2
        B[j] = (a - b) / 2
    end
    return A, B
end

function symmetrize!(A::Array, B::Array)
    size(A) != size(B) && error("A and B have different sizes")
    for i in eachindex(A)
        a, b = A[i], B[i]
        A[i] = (a + b) / 2
        B[i] = (a - b) / 2
    end
    return A, B
end


struct LiteSimulation <: Simulation
    constructor::AbstractMoleculeConstructor
    mol::Molecule
    eig
    μ
    Δ
    H0
    S0
    H
    HSOC
    γL
    ΓL
    ΔΓL
    ΓR
    ΔΓR
    L
end
Base.show(io::IO, sim::LiteSimulation) = Base.print(io, "LiteSimulation(...)")

@export function LiteSimulation(
    m;
    V = 0.0,
    fixed_distance = nothing,
    correct_overlaps = true,
    align_molecule = false,
    xyz_file = nothing,

    ρ0 = Diagonal([0.0, 0.0, 1.0, 0.0]),
    ρL = 1.0I(4),
    ρR = 1.0I(4),
    seed = 47,
    rng = MersenneTwister(seed),
    coupling_type = :random,

    αL = 1.0,
    αR = 1.0,
    lead_orbitals = (orbital"2s", orbital"2p_x", orbital"2p_y", orbital"2p_z"),

    λ = DEFAULT_λ,

    kwargs...
)
    mol = m(; kwargs...)
    xyz_file != nothing && Molecules.read_xyz!(mol, xyz_file)
    align_molecule && Molecules.balance_position!(GLOBAL_SKT, mol)

    z1, z2 = extrema(map(x -> Molecules.position(x)[3], mol.atoms))
    z0 = (z1 + z2) / 2
    d = (z2 - z1) / 2
    electric_field = (x, y, z) -> iszero(d) ? 0.0 : -V * (z - z0) / d

    H0 = Hermitian(hamiltonian(Float64, SLATERKOSTER_HAMILTONIAN, mol, fixed_distance=fixed_distance, electric_field=electric_field))
    L0 = map(x -> angularmomentum(x, mol), [:x, :y, :z])

    left_atom, right_atom = mol.atoms[m.left_index], mol.atoms[end-m.right_index+1]

    f = Molecules.slaterkoster(
        CC_ssσ = -2.0,
        CC_spσ = -1.5,
        CC_psσ = 1.5,
        CC_ppσ = 1.5,
        CC_ppπ = -0.5,
    )

    γL = αL * makecoupling(Val(coupling_type), rng, ρ0, left_atom, mol, f=f, seed=seed+1, m=-1)
    ΓLup = αL * makecoupling(Val(coupling_type), rng, ρL, left_atom, mol, f=f, seed=seed+1, m=-1)
    ΓLdown = αL * makecoupling(Val(coupling_type), rng, ρL, left_atom, mol, f=f, seed=seed+1, m=-1)
    ΓRup = αR * makecoupling(Val(coupling_type), rng, ρR, right_atom, mol, f=f, seed=seed+2, m=+1)
    ΓRdown = αR * makecoupling(Val(coupling_type), rng, ρR, right_atom, mol, f=f, seed=seed+2, m=+1)

    ΓL, ΔΓL = symmetrize!(ΓLup, ΓLdown)
    ΓR, ΔΓR = symmetrize!(ΓRup, ΓRdown)

    if correct_overlaps
        S0 = Hermitian(overlap(Float64, SLATERKOSTER_OVERLAP, mol, fixed_distance=fixed_distance))
        Ssqrtinv = sqrt(inv(S0))
        H = Hermitian(Ssqrtinv*H0*Ssqrtinv)
        L = map(L -> Ssqrtinv*L*Ssqrtinv, L0)
        ΓL, ΓR, ΔΓL, ΔΓR = map(Γ -> Hermitian(Ssqrtinv*Γ*Ssqrtinv), (ΓL, ΓR, ΔΓL, ΔΓR))
    else
        S0 = I
        H = H0
        L = L0
    end
    eig = eigen(H, sortby=real)
    μ, Δ = chemicalpotential(mol, eig.values)

    HSOC = λ * sum(L ⊗ σ for (L, σ) in zip(L, σ))

    return LiteSimulation(m, mol, eig, μ, Δ, H0, S0, H, HSOC, γL, ΓL, ΔΓL, ΓR, ΔΓR, L)
end

function makecoupling(::Val{:identity}, J::UniformScaling, ρ, atom, mol; err=10, kwargs...)
    indices = Molecules.indices(mol, atom)
    N = Molecules.countorbitals(mol)
    U = spzeros(Float64, size(ρ, 2), N)
    U[:, indices] = Matrix(J, size(ρ, 2), length(indices))
    B = U'ρ*U
    rmul!(B, rank(B)*inv(tr(B)))
    if !isposdef!(B + 1e-16I)
        err == 0 && error("Invalid coupling!")
        return makecoupling(Val(:identity), rng, ρ, atom, mol, err=err-1)
    end
    return B
end

function makecoupling(::Val{:atom}, rng::AbstractRNG, args...; seed=nothing, kwargs...)
    !isnothing(seed) && Random.seed!(rng, seed+1)
    makecoupling(Val(:atom), randn(rng, 3), args...; kwargs...)
end

makecoupling(::Val{:atom}, J::Float64, args...; m=1, kwargs...) = makecoupling(Val(:atom), [0.0, 0.0, m*J], args...; kwargs...)

function makecoupling(::Val{:atom}, d, ρ, atom, mol; f=SLATERKOSTER_HAMILTONIAN, kwargs...)
    lead_atom = Carbon(atom.position + d)
    lead_bond = Bond(atom, lead_atom)
    lead = Molecule([lead_atom])
    coupling = Set([lead_bond])
    U = hamiltonian(Float64, f, lead, mol, coupling)
    return U'ρ*U
end

function makecoupling(::Val{:random}, rng, ρ, atom, mol; seed=nothing, err=10, kwargs...)
    indices = Molecules.indices(mol, atom)
    N = Molecules.countorbitals(mol)
    U = spzeros(Float64, size(ρ, 2), N)
    !isnothing(seed) && Random.seed!(rng, seed+1)
    U[:, indices] = randn(rng, size(ρ, 2), length(indices))
    B = U'ρ*U
    #rmul!(B, rank(B)*inv(tr(B)))
    if !isposdef!(B + 1e-16I)
        err == 0 && error("Invalid coupling!")
        return makecoupling(Val(:random), rng, ρ, atom, mol, err=err-1)
    end
    return B
end

@export function fullpropagator(
    sim::Simulation;
    E = sim.μ,
    λ = DEFAULT_λ,
    Γ = (sim.ΓL + sim.ΓR) ⊗ σ0,
    H = sim.H ⊗ σ0,
    Λσ = sim.HSOC,
)
    return inv(E*I - H - Λσ + im*Γ/2)
end

@export function propagator(
    sim::Simulation;
    E = sim.μ,
    Γ = sim.ΓL + sim.ΓR,
    H = sim.H,
)
    return inv(E*I - H + im*Γ/2)
end

@export function reduced_propagator(sim::Simulation; E=sim.μ, α=DEFAULT_α, λ=DEFAULT_λ,
                                    ΓL=sim.ΓL, ΓR=sim.ΓR, H=sim.H, L=sim.L, a=[0, 0, 0], η=0.0)
    return inv(E*I - H - sum(a[i]*L[i] for i in 1:3) + α*λ*im/2*(ΓL + ΓR) + im*η*I)
end

@export commutator(A, B) = A*B - B*A
@export timeevolution(H, dt) = exp(-im*H*dt)
@export timeevolution_1(H, dt) = I - im*H*dt# - H*H*dt^2/2
@export function update_ρ(ρ0, U, dt)
    #U = exp(-im*H*dt)
    ρ1 = U*ρ0*U'
    #ρ1 = ρ0 + im*dt*(ρ0*H - H*ρ0)# + 1e-16*rand(size(ρ0)...)
    return ρ1# * tr(ρ0) / tr(ρ1)
end

@export function correlation(ρ, H, A, B, t1, t2)
    U1 = timeevolution(H, t1)
    U2 = timeevolution(H, t2)
    return -im*tr(ρ*commutator(U1'A*U1, U2'B*U2))
end

@export function equilibrium_density(sim::Simulation, H)
    eig = eigen((H + H')/2)
    Es = eig.values
    μ = sim.μ
    vecs = eig.vectors[:, Es .< μ]
    ρ0 = vecs*vecs'

    # remvecs = eig.vectors[:, Es .> μ]
    # remvecs = remvecs * (Diagonal(randn(size(remvecs, 2) >> 1)) ⊗ σ[0])
    # ρ′ = remvecs*remvecs'
    ρ = ρ0# + 1e-1*ρ′
    return ρ * tr(ρ0) / tr(ρ)
end

@export function perturbed_density(sim::Simulation, H)
    eig = eigen((H + H')/2)
    Es = eig.values
    μ = sim.μ
    vecs = eig.vectors[:, Es .> μ][:, 1:2]
    #vec = fill(0.0, length(sim.eig.values))
    #vec[2] = 1.0
    return vecs*vecs'
end

@export function equilibrium_density(sim::Simulation)
    Es = sim.eig.values
    μ = sim.μ
    vecs = sim.eig.vectors[:, Es .< μ]
    return vecs*vecs'
end

@export function fullhamiltonian(
    sim::Simulation;
    λ = DEFAULT_λ,
    H = sim.H ⊗ σ0,
    Lσ = sum(L ⊗ σ for (L, σ) in zip(sim.L, σ)),
)
    return let
        H = H + λ*Lσ
        (H .+ H') ./ 2
    end
end

@export function adiabatic_polarization(sim::Simulation)
    ρ0 = equilibrium_density(sim)
    return imag(tr(sim.ΓL*ρ0*sim.Lz))
end

@export function adiabatic_phase(sim2::Simulation, sim1::Simulation, dt)
    a1 = sim1.eig.vectors
    Es = sim1.eig.values
    vecs = eachcol(a1)
    a2 = sim2.eig.vectors
    V = (sim2.H0 - sim1.H0) ./ dt
    #return commutator(sim1.H0 ⊗ σ0, sim1.Lz ⊗ σz)
    return commutator(sim1.H0 ⊗ σ0, commutator(sim1.H0 ⊗ σ0, 4e-3 * sim1.Lz ⊗ σz))
    ΔEs = (sim2.eig.values .- sim1.eig.values)
    ΔEs_ = diag(a1'V*a1) * dt
    return ΔEs ./ ΔEs_
    G(E) = inv((E + 1e-16*im)*I - sim1.H0)
    b = [G(En)*(V - Diagonal(Edots))*n for (En, n) in zip(Es, vecs)]
    return real.(reduce(hcat, b)*dt)

    p = Diagonal(exp.(-im .* sim1.eig.values .* dt) ./ diag(a2'a1) ./ dt)
    return (a1'a2*p - I / dt)
end

function calculate_s_and_t(G, γL, ΓR, Ls; λ=DEFAULT_λ)
    X = G*γL*G'
    t = 2real(dot(X, ΓR))
    s = map(L -> 4λ*dot(X, ΓR*G*L), Ls)
    return s, t
end

@export function nearest_eigenstates(eig::Eigen, E; n=2, N=length(eig.values), nmin=1)
    vals = collect(enumerate(abs.(inv.(eig.values .- E))))
    sort!(vals, by=x->x[2], rev=true)
    indices = map(val -> val[1], vals)
    eigvals = eig.values[indices[nmin:(nmin+n-1)]]
    U = reduce(hcat, (eig.vectors[:, indices[i]] for i in nmin:(nmin+n-1)))
    V = reduce(hcat, (eig.vectors[:, indices[i]] for i in (nmin+n):N))
    return eigvals, U, V, indices
end

@export function eigenstate(sim::Simulation; E = sim.μ, ϕ=0.0)
    index_below = findlast(y -> E - y > 0, sim.eig.values)
    index_above = findfirst(y -> E - y <= 0, sim.eig.values)
    vals = sim.eig.values[index_below], sim.eig.values[index_above]
    vecs = sim.eig.vectors[:, index_below], sim.eig.vectors[:, index_above]
    x = (E - vals[1]) / (vals[2] - vals[1])
    weights = sqrt(1 - x), exp(im*ϕ)*sqrt(x)
    return sum(vecs .* weights)
end

@export function extract_slater_koster_params(l)
    mol = Molecule()
    A = Carbon([0.0, 0.0, 0.0], orbital"2s", orbital"2p_x", orbital"2p_y", orbital"2p_z")
    B = Carbon([0.0, 0.0, 1.0*l], orbital"2s", orbital"2p_x", orbital"2p_y", orbital"2p_z")
    C = Hydrogen([0.0, 0.0, 0.0], orbital"1s")
    D = Hydrogen([0.0, 0.0, 1.0*l], orbital"1s")
    push!(mol, A, B, C, D)
    push!(mol, Bond(A, B))
    push!(mol, Bond(B, C))
    push!(mol, Bond(C, D))
    push!(mol, Bond(D, A))
    H = hamiltonian(Float64, SLATERKOSTER_HAMILTONIAN, mol)
    println("H-H")
    println("E_s = $(H[9,9])")
    println("V_ssσ = $(H[9,10])")

    println("C-C")
    println("E_s = $(H[1,1])")
    println("E_p = $(H[2,2])")
    println("V_ssσ = $(H[1,5])")
    println("V_spσ = $(H[1,8])")
    println("V_psσ = $(H[4,5]) ($(H[4,5] == -H[1,8] ? "" : "in")consistent)")
    println("V_ppσ = $(H[4,8])")
    println("V_ppπ = $(H[2,6])")

    println("H-C")
    println("V_ssσ = $(H[5,9])")
    println("V_spσ = $(H[9,8])")
    println("V_psσ = $(H[4,10]) ($(H[4,10] == -H[8,9] ? "" : "in")consistent)")
end
