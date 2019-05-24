projectionmatrix(us...) = sum(u*u' for u in map(normalize, us)) / length(us)

struct SimpleSimulation
    N::Integer
    mol::Molecule
    μ
    Δ
    H0
    ΓL
    ΓR
    ΔΓL
    ΔΓR
    Lx
    Ly
    Lz
    Lxs
    Lys
    Lzs
end

function SimpleSimulation(N::Integer; η = 1e-15, ϕ = π/4, l=1.4, handedness=1, α=1., Γ=1/27, θ=nothing, δz=nothing)
    # Construct molecule and leads
    mol = makemol(N, ϕ, l, handedness=handedness, θ=θ, δz=δz)
    skt = loaddir(joinpath(PARAMDIR, "ob2-1-1", "base"))
    hydrogen_atoms = collect(filter(x -> x isa Hydrogen, mol.atoms))
    carbon_atoms = collect(filter(x -> x isa Carbon, mol.atoms))

    H0 = hamiltonian(Float64, skt, mol)
    energies, eigvecs = molecularorbitals(H0)
    μ, Δ = real.(chemicalpotential(mol, energies))

    # ΓL = Γ * sum((u = Vector{Float64}(mol, carbon_atoms[1], i); u*u'/3) for i in 1:4)

    # ΓR = Γ * sum((u = Vector{Float64}(mol, carbon_atoms[end], i); u*u'/3) for i in 1:4)

    Random.seed!(47)
    ΓL = 1e-2*Γ * sum((u = Vector{Float64}(mol, carbon_atoms[1+j], i); u*u') for i in 1:4, j in 0:0)

    ΓR1 = 1e-2*Γ * sum((u = Vector{Float64}(mol, carbon_atoms[end-j], i); (rand() + 0.5)*u*u') for i in 1:4, j in 0:0)
    ΓR2 = 1e-2*Γ * sum((u = Vector{Float64}(mol, carbon_atoms[end-j], i); (rand() + 0.5)*u*u') for i in 1:4, j in 0:0)
    ΓR = (ΓR1 .+ ΓR2) ./ 2
    ΔΓR = (ΓR1 .- ΓR2) ./ 2

    Lxs = [angularmomentum(:x, mol, i) for (i, atom) in enumerate(mol) if atom isa Carbon]
    Lys = [angularmomentum(:y, mol, i) for (i, atom) in enumerate(mol) if atom isa Carbon]
    Lzs = [angularmomentum(:z, mol, i) for (i, atom) in enumerate(mol) if atom isa Carbon]
    Lx = angularmomentum(:x, mol)
    Ly = angularmomentum(:y, mol)
    Lz = angularmomentum(:z, mol)

    return SimpleSimulation(N, mol, μ, Δ, H0, ΓL, ΓR, 0I, ΔΓR, Lx, Ly, Lz, Lxs, Lys, Lzs)
end

function HeliceneSimulation(N::Integer; η = 1e-15, ϕ = 0.5, l=1.4, handedness=1, α=1., Γ=1/27, oxidation=0, electric_field=(x, y, z) -> 0.0)
    # Construct molecule and leads
    mol = makehelicene2(N, ϕ, l, handedness=handedness)
    skt = loaddir(joinpath(PARAMDIR, "ob2-1-1", "base"))
    hydrogen_atoms = collect(filter(x -> x isa Hydrogen, mol.atoms))
    carbon_atoms = collect(filter(x -> x isa Carbon, mol.atoms))
    oxygen_atoms = collect(filter(x -> x isa Oxygen, mol.atoms))
    nitrogen_atoms = collect(filter(x -> x isa Nitrogen, mol.atoms))

    H0 = hamiltonian(Float64, skt, mol, electric_field=electric_field)
    energies, eigvecs = molecularorbitals(H0)
    μ, Δ = real.(chemicalpotential(mol, energies, oxidation))

    Random.seed!(47)
    first_indices = Molecules.indices(mol, oxygen_atoms[1])
    last_indices = Molecules.indices(mol, nitrogen_atoms[end])
    HLM = begin
        U = fill(0.0, size(H0))
        A = rand(4, 4)
        U[first_indices, first_indices] = (A + A') / 2
        U
    end
    HRM = begin
        U = fill(0.0, size(H0))
        A = rand(4, 4)
        U[last_indices, last_indices] = (A .+ A') ./ 2
        U ./ tr(U)
    end

    ΓL = 1e-2*Γ * sum((u = Vector{Float64}(mol, oxygen_atoms[1+j], i); u*u') for i in 1:4, j in 0:0)
    #ΓL .= HLM'ΓL*HLM

    ΓR1 = 1e-2*Γ * sum((u = Vector{Float64}(mol, oxygen_atoms[end-j], i); u*u') for i in 1:4, j in 0:0)
    ΓR2 = 1e-2*Γ * sum((u = Vector{Float64}(mol, carbon_atoms[end-j], i); u*u') for i in 1:4, j in 0:0)
    #ΓR2 .= HRM'ΓR2*HRM
    ΓR = (ΓR1 .+ ΓR2) ./ 2
    ΔΓR = (ΓR1 .- ΓR2) ./ 2

    Lxs = [angularmomentum(:x, mol, i) for (i, atom) in enumerate(mol) if typeof(atom) <: Union{Carbon, Nitrogen, Oxygen}]
    Lys = [angularmomentum(:y, mol, i) for (i, atom) in enumerate(mol) if typeof(atom) <: Union{Carbon, Nitrogen, Oxygen}]
    Lzs = [angularmomentum(:z, mol, i) for (i, atom) in enumerate(mol) if typeof(atom) <: Union{Carbon, Nitrogen, Oxygen}]
    Lx = angularmomentum(:x, mol)
    Ly = angularmomentum(:y, mol)
    Lz = angularmomentum(:z, mol)

    return SimpleSimulation(N, mol, μ, Δ, H0, ΓL, ΓR, 0I, ΔΓR, Lx, Ly, Lz, Lxs, Lys, Lzs)
end

function new_coupling(rs, mol; γ = 1/27)
    carbon_atoms = collect(filter(x -> x isa Carbon, mol.atoms))
    uL = [Vector{Float64}(mol, carbon_atoms[1], i) for i in 2:4]
    k = sum(r*u for (u, r) in zip(uL, rs))
    return γ*k*k'
end

function theorem_f(E, m=0; sim, λ=1e-2 / 27, rs = normalize!([(0.5 - rand()) for _ in 1:3]), γ=1/27
)
    ΓL, ΓR, γL = sim.ΓL, sim.ΓR, new_coupling(rs, sim.mol, γ=γ)
    H0 = Matrix(sim.H0)
    ΓA = ΓL*(1 - m) + γL*m

    Λ = λ*sim.Lz
    G = inv(E*I - H0 + im*(ΓL + ΓR)/2)
    A = -im.*(G .- G') ./ 2
    t = real(dot(G*ΓA*G', ΓR))
    s = 2real(dot(ΓA, G'ΓR*G*Λ*G))
    # function T(b, Ln)
    #     G = inv(E*I - H0 + im*(ΓL + ΓR)/2 - λ*b*Ln)
    #     return real(dot(G*ΓA*G', ΓR))
    # end
    # p = derivative(x -> log(T(x, sim.Lz)), 0.0)

    return t, s, A
end

function theorem_f2(E, m=0; sim, λ=1e-2 / 27, rs = normalize!([(0.5 - rand()) for _ in 1:3]), γ=1/27
)
    ΓL, ΓR, γL = sim.ΓL, sim.ΓR, new_coupling(rs, sim.mol, γ=γ)
    H0 = Matrix(sim.H0)
    ΓA = ΓL*(1 - m) + γL*m

    Λ = λ*sim.Lz
    G = inv(E*I - H0 + im*(ΓL + ΓR)/2)
    t = real(dot(G*ΓA*G', ΓR))
    s = 2real(dot(ΓA, G'ΓR*G*Λ*G))
    # function T(b, Ln)
    #     G = inv(E*I - H0 + im*(ΓL + ΓR)/2 - λ*b*Ln)
    #     return real(dot(G*ΓA*G', ΓR))
    # end
    # p = derivative(x -> log(T(x, sim.Lz)), 0.0)

    return t, s / t
end

function theorem_f2_each(E, m=0; sim, λ=1e-2 / 27, rs = normalize!([(0.5 - rand()) for _ in 1:3])
)
    ΓL, ΓR, γL = sim.ΓL, sim.ΓR, new_coupling(rs, sim.mol)
    H0 = Matrix(sim.H0)
    ΓA = ΓL*(1 - m) + γL*m

    function T(b, Ln)
        G = inv(E*I - H0 + im*(ΓL + ΓR)/2 - λ*b*Ln)
        return real(dot(G*ΓA*G', ΓR))
    end

    t = T(0, sim.Lz)
    return t, [(
        derivative(x -> log(T(x, Lx)), 0.0),
        derivative(x -> log(T(x, Ly)), 0.0),
        derivative(x -> log(T(x, Lz)), 0.0))
            for (Lx, Ly, Lz) in zip(sim.Lxs, sim.Lys, sim.Lzs)]
end

function theorem_f2_approx(E, m=0; sim, λ=1e-2 / 27, rs = normalize!([(0.5 - rand()) for _ in 1:3])
)
    ΓL, ΓR, γL = sim.ΓL, sim.ΓR, new_coupling(rs, sim.mol)
    H0 = Matrix(sim.H0)
    G = inv(E*I - H0 + im*(ΓL + ΓR)/2)
    eig = eigen(E*I - H0)

    valmin =
        #(0, 36)
        findmin(inv.(eig.values))
    valmax =
        #(0, 37)
        findmax(inv.(eig.values))

    us = [
        #eig.vectors[:, valmin[2] - 1],
        eig.vectors[:, valmin[2]],
        eig.vectors[:, valmax[2]]
        #eig.vectors[:, valmax[2] + 1]
    ]
    #us = [eig.vectors[:, i] for i in 1:72]
    P = sum(u*u' for u in us)

    ΓA = ΓL*(1 - m) + γL*m

    local T, d
    let P = P, H0 = P*H0*P, ΓL = P*ΓL*P, ΓR = P*ΓR*P, ΓA = P*ΓA*P
        function T(b, Ln)
            Ln = P*Ln*P
            G = P*inv((E*I - H0 + im*(ΓL + ΓR)/2 - λ*b*Ln))*P
            return real(dot(G*ΓA*G', ΓR))
        end
        d = dot(ΓL, ΓA) / dot(ΓL, ΓL)#d = sqrt(tr((ΓL - ΓA).^2)) / sqrt(tr(ΓL.^2))
    end

    t = T(0, sim.Lz)
    p = derivative(x -> log(T(x, sim.Lz)), 0.0)
    return t, p*t, d, valmin, valmax
end

function theorem_f2_approx2(E, m=0; sim, λ=1e-2 / 27, rs = normalize!([(0.5 - rand()) for _ in 1:3])
)
    ΓL, ΓR, γL = sim.ΓL, sim.ΓR, new_coupling(rs, sim.mol)
    H0 = Matrix(sim.H0)
    G = inv(E*I - H0 + im*(ΓL + ΓR)/2)
    eig = eigen(E*I - H0)

    vals = collect(enumerate(abs.(inv.(eig.values))))
    sort!(vals, by=x->x[2])
    valmin =
        #(0, 36)
        #findmin(inv.(eig.values))
        vals[end]
    valmax =
        #(0, 37)
        #findmax(inv.(eig.values))
        vals[end-1]

    us = [
        #eig.vectors[:, valmin[2] - 1],
        eig.vectors[:, valmin[1]],
        eig.vectors[:, valmax[1]]
        #eig.vectors[:, valmax[2] + 1]
    ]
    P = sum(u*u' for u in us)

    ΓA = ΓL*(1 - m) + γL*m

    G = inv(E*I - H0 + im*(ΓL + ΓR)/2)
    Λ = λ*sim.Lz
    A = -im.*(G .- G') ./ 2

    #t = real(tr(ΓA*P*G'P*ΓR*P*G*P))


    Λ12 = us[1]'Λ*us[2]

    A1 = us[1]'A*us[1]
    A2 = us[2]'A*us[2]
    A12 = us[1]'A*us[2]
    A = [A1 A12; A12 A2]

    γ11 = us[1]'ΓA*us[1]
    γ12 = us[1]'ΓA*us[2]
    γ22 = us[2]'ΓA*us[2]
    γ = [γ11 γ12; γ12 γ22]

    Γ11 = us[1]'ΓR*us[1]
    Γ12 = us[1]'ΓR*us[2]
    Γ22 = us[2]'ΓR*us[2]
    Γ = [Γ11 Γ12; Γ12 Γ22]

    prefactor = real(sum(A[1,i]*A[j,k]*A[3-i,2]*(γ[j,i]*Γ[3-i,k] - γ[j,3-i]*Γ[i,k]) for i in 1:2, j in 1:2, k in 1:2))

    return 2imag(Λ12), prefactor
end

function theorem_magnet(E, m=0; sim, λ=1e-2 / 27, rs = normalize!([(0.5 - rand()) for _ in 1:3]), γ=1/27
)
    ΓL, ΓR, ΔΓR = sim.ΓL, sim.ΓR, sim.ΔΓR
    H0 = Matrix(sim.H0)
    ΓRup = (ΓR .* (2 - m) .+ m .* ΔΓR)
    ΓRdown = (ΓR .* (2 - m) .- m .* ΔΓR)

    Λ = λ*sim.Lz
    G = inv(E*I - H0 + im*(ΓL + ΓR)/2)
    t = real(dot(ΓL, G'ΓR*G))
    #s = real(dot(ΓL, G'ΔΓR*G*Λ*G))
    t += 2imag(dot(ΓL, G'ΔΓR*G*ΔΓR*G))
    s = real(dot(ΓRup, G'ΓRdown*G*Λ*G))

    return t, s
end

function β(γ, a, b)
    A = similar(γ)
    N, M = size(A)
    for j in 1:M, i in 1:N
        A[j, i] = (j == a && i == b) || (j == b && i == a) ? 1.0 : 0.0
    end
    return A
end

export theorem_derivative
function theorem_derivative(E, θ, ϕ; sim, λ=1e-2 / 27, rs = normalize!([(0.5 - rand()) for _ in 1:3]), γ=1/27
)
    ΓL, ΓR = sim.ΓL, sim.ΓR
    H0 = Matrix(sim.H0)
    #βs = [-sin(ϕ)*sin(θ), sin(ϕ)*cos(θ), cos(ϕ)]
    βs = [cos(ϕ)*cos(θ), cos(ϕ)*sin(θ), -sin(ϕ)]
    γ = sum(βs[i] * β(sim.ΓL, i+2, i+2) for i in 1:3)

    G = inv(E*I - H0 + im*(ΓL + ΓR)/2)
    Λ = λ*sim.Lz

    t = real(dot(γ, G'ΓR*G))
    s = 2*real(dot(γ, G'ΓR*G*Λ*G))

    return t, s, s / t
end

## Lattice

struct LatticeSimulation
    uc::Molecules.UnitCell
    H0
    Lx
    Ly
    Lz
end

function LatticeSimulation(; η = 1e-15, ϕ = π/4, l=1.4, handedness=1, α=1., Γ=1/27, θ=nothing, δz=nothing, ξs = range(0, π, length=200))
    skt = loaddir(joinpath(PARAMDIR, "ob2-1-1", "base"))
    uc = makelattice(ϕ, l)
    Lx = Molecules.angularmomentum(:x, uc.unit)
    Ly = Molecules.angularmomentum(:y, uc.unit)
    Lz = Molecules.angularmomentum(:z, uc.unit)

    H0(ξ) = hamiltonian(Complex{Float64}, skt, uc, ξ)

    return LatticeSimulation(uc, H0, Lx, Ly, Lz)
end
