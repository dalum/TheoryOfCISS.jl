function eigensort(H; by=identity)
    vals, vecs = eigen(H)
    s = sortperm(vals, by=by)
    return Eigen(vals[s], vecs[:, s])
end

#const GLOBAL_SKT = loaddir(joinpath(PARAMDIR, "ob2-1-1", "base"))
const GLOBAL_SKT = loaddir(joinpath(PARAMDIR, "mio-1-1"))

using Molecules: σ
@export const ⊗ = kron

struct LiteSimulation <: Simulation
    N::Integer
    mol::Molecule
    eig
    μ
    Δ
    H0
    γL
    ΓL
    ΓR
    Lx
    Ly
    Lz
end
Base.show(io::IO, sim::LiteSimulation) = Base.print(io, "LiteSimulation($(sim.N))")

@export function LiteSimulation(
    N::Integer;
    η = 1e-15,
    ϕ = π/4,
    l = 1.4,
    handedness = 1,
    α = 1.0,
    V = 0.0,
    θ = nothing,
    δz = nothing,
    fixed_distance = nothing,
    ρ0 = Diagonal([0.0, 1.0, 0.0]),
    ρL = Matrix(1.0I, 6, 6),
    ρR = Matrix(1.0I, 6, 6),
    rng = MersenneTwister(47)
)
    mol = makemol(N, ϕ, l, handedness=handedness, θ=θ, δz=δz)
    hydrogen_atoms = collect(filter(x -> x isa Hydrogen, mol.atoms))
    carbon_atoms = collect(filter(x -> x isa Carbon, mol.atoms))
    z1, z2 = extrema(map(x->Molecules.position(x)[3], mol.atoms))
    z0 = (z1 + z2)/2
    d = (z2 - z1)/2
    electric_field = (x, y, z) -> iszero(d) ? 0.0 : -V*(z - z0)/d

    H0 = hamiltonian(Float64, GLOBAL_SKT, mol, fixed_distance=fixed_distance, electric_field=electric_field)
    eig = eigensort(H0, by=real) # eigen(H0, sortby=real)
    μ, Δ = chemicalpotential(mol, eig.values)

    γL = sparse(α*makecoupling(rng, ρ0, H0, Molecules.indices(mol, carbon_atoms[1])))
    ΓL = sparse(α*makecoupling(rng, ρL, H0, Molecules.indices(mol, carbon_atoms[1])))
    ΓR = sparse(α*makecoupling(rng, ρR, H0, Molecules.indices(mol, carbon_atoms[end])))
    Lx = sparse(angularmomentum(:x, mol))
    Ly = sparse(angularmomentum(:y, mol))
    Lz = sparse(angularmomentum(:z, mol))

    return LiteSimulation(N, mol, eig, μ, Δ, H0, γL, ΓL, ΓR, Lx, Ly, Lz)
end

function makecoupling(rng, ρ, A, indices)
    U = fill(0.0, size(ρ, 2), size(A, 1))
    U[:, indices] = randn(rng, size(ρ, 2), length(indices))
    B = U'ρ*U
    return rmul!(B, inv(norm(B)))
end

@export propagator(sim::Simulation, E) = inv(E*I - sim.H0 + im*(sim.ΓL + sim.ΓR)/2)
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

@export function fullhamiltonian(sim::Simulation; λ=4e-3)
    H = complex(sim.H0 ⊗ σ[0])
    H += -im * sim.ΓR ⊗ σ[0]
    H += sum(λ * L ⊗ Molecules.σ[x] for (L, x) in zip((sim.Lx, sim.Ly, sim.Lz), (:x, :y, :z)))
    return H
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
    return (a1'a2*p - I/dt)
end

@export function polarization(G, γL, ΓR, Λ)
    X, Y = G*γL, G'ΓR
    t = real(dot(X', Y))
    s = 2*real(dot(X', Y*G*Λ))
    return s / t, t, s
end

@export function nearest_eigenstates(eig::Eigen, E)
    vals = collect(enumerate(abs.(inv.(eig.values .- E))))
    sort!(vals, by=x->x[2])
    eigvals = eig.values[vals[end][1]], eig.values[vals[end-1][1]]
    U = [eig.vectors[:,vals[end][1]] eig.vectors[:,vals[end-1][1]]]
    return eigvals, U
end

@export function extract_slater_koster_params(l)
    CCmol = Molecule()
    A = Carbon([0.0, 0.0, 0.0], orbital"2s", orbital"2p_x", orbital"2p_y", orbital"2p_z")
    B = Carbon([0.0, 0.0, 1.0*l], orbital"2s", orbital"2p_x", orbital"2p_y", orbital"2p_z")
    push!(CCmol, A, B)
    push!(CCmol, Bond(A, B))
    hamiltonian(Float64, GLOBAL_SKT, CCmol)
end
