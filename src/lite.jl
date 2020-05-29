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
    H
    H0
    S0
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
    ρL = 1.0I,
    ρR = 1.0I,
    seed = 47,
    rng = MersenneTwister(seed),

    αL = 1.0,
    αR = 1.0,

    kwargs...
)
    mol = m(; kwargs...)
    xyz_file != nothing && Molecules.read_xyz!(mol, xyz_file)
    align_molecule && Molecules.balance_position!(GLOBAL_SKT, mol)

    z1, z2 = extrema(map(x -> Molecules.position(x)[3], mol.atoms))
    z0 = (z1 + z2) / 2
    d = (z2 - z1) / 2
    electric_field = (x, y, z) -> iszero(d) ? 0.0 : -V * (z - z0) / d

    H0 = Hermitian(hamiltonian(Float64, GLOBAL_SKT, mol, fixed_distance=fixed_distance, electric_field=electric_field))
    if correct_overlaps
        S0 = Hermitian(overlap(Float64, GLOBAL_SKT, mol, fixed_distance=fixed_distance))
        Ssqrtinv = sqrt(inv(S0))
        H = Ssqrtinv*H0*Ssqrtinv
    else
        S0 = I
        H = H0
    end
    eig = eigen(H, sortby=real)
    μ, Δ = chemicalpotential(mol, eig.values)

    left_atom, right_atom = mol.atoms[m.left_index], mol.atoms[end-m.right_index+1]

    γL = αL .* makecoupling(rng, ρ0, H0, Molecules.indices(mol, left_atom))
    ΓLup = αL .* makecoupling(rng, ρL, H0, Molecules.indices(mol, left_atom))
    ΓLdown = αL .* makecoupling(rng, ρL, H0, Molecules.indices(mol, left_atom))
    ΓRup = αR .* makecoupling(rng, ρR, H0, Molecules.indices(mol, right_atom))
    ΓRdown = αR .* makecoupling(rng, ρR, H0, Molecules.indices(mol, right_atom))

    ΓL, ΔΓL = symmetrize!(ΓLup, ΓLdown)
    ΓR, ΔΓR = symmetrize!(ΓRup, ΓRdown)

    L = map(x -> angularmomentum(x, mol), [:x, :y, :z])

    return LiteSimulation(m, mol, eig, μ, Δ, H, H0, S0, γL, ΓL, ΔΓL, ΓR, ΔΓR, L)
end

makecoupling(rng, ρ::UniformScaling, A, indices; kwargs...) = makecoupling(rng, Matrix(1.0I, length(indices), length(indices)), A, indices; kwargs...)
makecoupling(rng::UniformScaling, ρ::UniformScaling, A, indices; kwargs...) = makecoupling(rng, Matrix(1.0I, length(indices), length(indices)), A, indices; kwargs...)

function makecoupling(rng::UniformScaling, ρ, A, indices; err=10)
    U = spzeros(Float64, size(ρ, 2), size(A, 1))
    U[:, indices] = Matrix(rng, size(ρ, 2), length(indices))
    B = U'ρ*U
    rmul!(B, rank(B)*inv(tr(B)))
    if !isposdef!(B + 1e-16I)
        err == 0 && error("Invalid coupling!")
        return makecoupling(rng, ρ, A, indices, err=err-1)
    end
    return B
end

function makecoupling(rng, ρ, A, indices; err=10)
    U = spzeros(Float64, size(ρ, 2), size(A, 1))
    U[:, indices] = randn(rng, size(ρ, 2), length(indices))
    B = U'ρ*U
    rmul!(B, rank(B)*inv(tr(B)))
    if !isposdef!(B + 1e-16I)
        err == 0 && error("Invalid coupling!")
        return makecoupling(rng, ρ, A, indices, err=err-1)
    end
    return B
end

@export function fullpropagator(
    sim::Simulation;
    E = sim.μ,
    α = DEFAULT_α,
    λ = DEFAULT_λ,
    Γ = (sim.ΓL + sim.ΓR) ⊗ σ0,
    H = sim.H ⊗ σ0,
    Λσ = λ * sum(L ⊗ σ for (L, σ) in zip(sim.L, σ)),
)
    return inv(E*I - H + α * λ * im * Γ / 2 - Λσ)
end

@export function propagator(sim::Simulation; E=sim.μ, α=DEFAULT_α, λ=DEFAULT_λ,
                            ΓL=sim.ΓL, ΓR=sim.ΓR, H=sim.H)
    return inv(E*I - H + α*λ*im/2*(ΓL + ΓR))
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

@export function fullhamiltonian(sim::Simulation; λ=DEFAULT_λ)
    H = complex(sim.H ⊗ σ0)
    #H += -im * sim.ΓR ⊗ σ0
    H += sum(λ * Li ⊗ σi for (Li, σi) in zip(sim.L, σ))
    return (H .+ H') ./ 2
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

function calculate_s_and_t(G, γL, ΓR, Ls)
    X = G*γL*G'
    t = 2real(dot(X, ΓR))
    s = map(L -> 4dot(X, ΓR*G*L), Ls)
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
    H = hamiltonian(Float64, GLOBAL_SKT, mol)
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
