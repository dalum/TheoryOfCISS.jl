struct SimData_f
    N
    ϕs
    Es
    Ts
    Ps
end

struct SimData_f_eigs
    N
    touched
    ϕs
    Es
    Ts
    Ss
    eigvals
    eigvecs
    Λs
    As
end

struct SimData_f_approx
    N
    ϕs
    Es
    Ts
    Ps
    A1s
    A2s
    Λs
    γs
    Γs
end

struct SimData_f_each
    N
    ϕs
    Es
    Ts
    Pxs
    Pys
    Pzs
end

struct SimData_f_contour
    N
    n
    ϕs
    Ts
    Pxs
    Pys
    Pzs
end

# ## dat14_hires
const data = Dict{Symbol,Any}()

@export function gen_datN_hires(
    N = 14;
    E_len = 1000,
    E_lower = -0.22,#0.26,
    E_upper = -0.17,#0.36,
    ϕ_len = 1000,
    ϕ_lower = 0.8,#1.0,
    ϕ_upper = 1.15,
    Es = range(E_lower, stop=E_upper, length=E_len),
    ϕs = range(ϕ_lower, stop=ϕ_upper, length=ϕ_len),
    clear = false
)
    if !(Symbol("dat$(N)_hires") in keys(data)) || clear
        Ts = fill(0.0, ϕ_len, E_len)
        Ps = fill(0.0, ϕ_len, E_len)
        data[Symbol("dat$(N)_hires")] = SimData_f(N, ϕs, Es, Ts, Ps)
    end
    dat = data[Symbol("dat$(N)_hires")]

    @showprogress 1 "" for (i,ϕ) in enumerate(ϕs)
        sim = SimpleSimulation(N, ϕ=ϕ)
        for (j,E) in enumerate(Es)
            t, p = theorem_f2(E, 1; sim=sim, rs=[0, 1, 0])
            dat.Ts[i,j] = t
            dat.Ps[i,j] = p
        end
    end
    JLD.save("dat$(N)_hires.jld", "dat", dat)
    return dat
end

@export function gen_helicene_datN_hires(
    N = 4;
    E_len = 50,
    E_lower = -0.25,#0.26,
    E_upper = -0.1,#0.36,
    ϕ_len = 100,
    ϕ_lower = 0.5,#1.0,
    ϕ_upper = 1.0,
    Es = range(E_lower, stop=E_upper, length=E_len),
    ϕs = range(ϕ_lower, stop=ϕ_upper, length=ϕ_len),
    clear = false
)
    if !(Symbol("helicene_dat$(N)_hires") in keys(data)) || clear
        Ts = fill(0.0, ϕ_len, E_len)
        Ps = fill(0.0, ϕ_len, E_len)
        data[Symbol("helicene_dat$(N)_hires")] = SimData_f(N, ϕs, Es, Ts, Ps)
    end
    dat = data[Symbol("helicene_dat$(N)_hires")]

    @showprogress 1 "" for (i,ϕ) in enumerate(ϕs)
        Ts, Ps = gen_helicene_datN_inner(Es, N, ϕ)
        dat.Ts[i,:] .= Ts
        dat.Ps[i,:] .= Ps
    end
    JLD.save("helicene_dat$(N)_hires.jld", "dat", dat)
    return dat
end

function gen_helicene_datN_inner(Es, N, ϕ)
    Ts = fill(0.0, length(Es))
    Ps = fill(0.0, length(Es))
    sim = HeliceneSimulation(N, ϕ=ϕ, Γ=1/27, oxidation=+1)
    for (j,E) in enumerate(Es)
        t, p = theorem_magnet(E, 1; sim=sim, rs=normalize!([0.0, 1.0, 0.0]))
        Ts[j] = t
        Ps[j] = p
    end
    return Ts, Ps
end


@export function gen_datN_hires_approx(
    N = 14;
    E_len = 200,
    E_lower = -0.22,
    E_upper = -0.17,
    ϕ_len = 200,
    ϕ_lower = 0.8,
    ϕ_upper = 1.15,
    Es = range(E_lower, stop=E_upper, length=E_len),
    ϕs = range(ϕ_lower, stop=ϕ_upper, length=ϕ_len),
    clear = false
)
    if !(Symbol("dat$(N)_hires_approx") in keys(data)) || clear
        Ts = fill(0.0, ϕ_len, E_len)
        Ps = fill(0.0, ϕ_len, E_len)
        data[Symbol("dat$(N)_hires_approx")] = SimData_f(N, ϕs, Es, Ts, Ps)
    end
    dat = data[Symbol("dat$(N)_hires_approx")]

    @showprogress 1 "" for (i,ϕ) in enumerate(ϕs)
        sim = SimpleSimulation(N, ϕ=ϕ)
        for (j,E) in enumerate(Es)
            t, p = theorem_f2_approx(E, 1; sim=sim, rs=[0, 1, 0])
            dat.Ts[i,j] = t
            dat.Ps[i,j] = p
        end
    end
    JLD.save("dat$(N)_hires_approx.jld", "dat", dat)
    return dat
end

@export function gen_datN_hires_approx2(
    N = 14;
    E_len = 200,
    E_lower = -0.22,
    E_upper = -0.17,
    ϕ_len = 200,
    ϕ_lower = 0.8,
    ϕ_upper = 1.15,
    Es = range(E_lower, stop=E_upper, length=E_len),
    ϕs = range(ϕ_lower, stop=ϕ_upper, length=ϕ_len),
    clear = false
)
    if !(Symbol("dat$(N)_hires_approx2") in keys(data)) || clear
        Ts = fill(0.0, ϕ_len, E_len)
        Ps = fill(0.0, ϕ_len, E_len)
        data[Symbol("dat$(N)_hires_approx2")] = SimData_f(N, ϕs, Es, Ts, Ps)
    end
    dat = data[Symbol("dat$(N)_hires_approx2")]

    @showprogress 1 "" for (i,ϕ) in enumerate(ϕs)
        sim = SimpleSimulation(N, ϕ=ϕ)
        for (j,E) in enumerate(Es)
            t, p = theorem_f2_approx2(E, 1; sim=sim, rs=[0, 1, 0])
            dat.Ts[i,j] = t
            dat.Ps[i,j] = p
        end
    end
    JLD.save("dat$(N)_hires_approx2.jld", "dat", dat)
    return dat
end

@export function gen_datN_hires_each(
    N = 14;
    E_len = 500,
    E_lower = -0.22,
    E_upper = -0.17,
    ϕ_len = 500,
    ϕ_lower = 0.82,
    ϕ_upper = 1.15,
    Es = range(E_lower, stop=E_upper, length=E_len),
    ϕs = range(ϕ_lower, stop=ϕ_upper, length=ϕ_len),
    clear = false
)
    if !(Symbol("dat$(N)_hires_each") in keys(data)) || clear
        Ts = fill(0.0, ϕ_len, E_len)
        Pxs = [fill(0.0, ϕ_len, E_len) for _ in 1:N]
        Pys = [fill(0.0, ϕ_len, E_len) for _ in 1:N]
        Pzs = [fill(0.0, ϕ_len, E_len) for _ in 1:N]
        data[Symbol("dat$(N)_hires_each")] = SimData_f_each(N, ϕs, Es, Ts, Pxs, Pys, Pzs)
    end
    dat = data[Symbol("dat$(N)_hires_each")]

    @showprogress 1 "" for (i,ϕ) in enumerate(ϕs)
        sim = SimpleSimulation(N, ϕ=ϕ)
        for (j,E) in enumerate(Es)
            t, ps = theorem_f2_each(E, 1; sim=sim, rs=[0, 1, 0])
            dat.Ts[i,j] = t
            for k in 1:N
                dat.Pxs[k][i,j] = ps[k][1]
                dat.Pys[k][i,j] = ps[k][2]
                dat.Pzs[k][i,j] = ps[k][3]
            end
        end
    end
    JLD.save("dat$(N)_hires_each.jld", "dat", dat)
    return dat
end

function gen_Ncontour(
    N = 14;
    n = floor(Int, N*2 + (N/2 + 1)),
    ϕ_len = 650,
    ϕ_lower = π/5,
    ϕ_upper = π/2.5,
    ϕs = range(ϕ_lower, stop=ϕ_upper, length=ϕ_len),
    clear = false
)
    if !(Symbol("dat$(N)_contour") in keys(data)) || clear
        Ts = fill(0.0, ϕ_len, E_len)
        Ps = fill(0.0, ϕ_len, E_len)
        data[Symbol("dat$(N)_hires")] = SimData_f(N, ϕs, Es, Ts, Ps)
    end
    dat = data[Symbol("dat$(N)_hires")]

    @showprogress 1 "" for (i,ϕ) in enumerate(ϕs)
        sim = SimpleSimulation(N, ϕ=ϕ)
        for (j,E) in enumerate(Es)
            p, t = theorem_f2(E + sim.μ, 1; sim=sim, rs=[0, 1, 0])
            dat.Ps[i,j] = p
            dat.Ts[i,j] = t
        end
    end
    JLD.save("dat$(N)_contour_$(ϕ_lower)_$(ϕ_upper).jld", "dat", dat)
    return dat
end

@export function gen_latticeΛs(ϕ, i; l=2.0, ξs = range(0, π, length=200), λ=1e-2 / 27)
    sim = LatticeSimulation(ϕ=ϕ, l=l)

    Λ = sim.Lz*λ

    eigs = Eigen[]
    for ξ in ξs
        eig = eigen(sim.H0(ξ))
        p = sortperm(eig.values, by=real)
        push!(eigs, Eigen(eig.values[p], eig.vectors[:, p]))
    end
    evs = hcat(map(x->x.vectors[:,i], eigs)...)
    vals = hcat(map(x->real(x.values[i]), eigs)...)'

    Λs = [2*real.(v1)'Λ*real.(v2) for v1 in eachcol(evs), v2 in eachcol(evs)]
    return Λs, vals, ξs
end

# function gen_near_bands(
#     N = 14;
#     l = 1.4,
#     E_len = 8000,
#     E_lower = 10/27.2,
#     E_upper = 22/27.2,
#     θ_len = 8000,
#     θ_lower = π/2,
#     θ_upper = π,
#     Es = range(E_lower, stop=E_upper, length=E_len),
#     θs = range(θ_lower, stop=θ_upper, length=θ_len),
#     clear = false,
#     threshold = 0.001,
#     Γ = 1/27.2
# )##:
#     if !(Symbol("dat$(N)_near_bands") in keys(data)) || clear
#         touched = fill(false, θ_len, E_len)
#         Ts = fill(0.0, θ_len, E_len)
#         Ss = fill(0.0, θ_len, E_len)
#         eigvals = Vector{Vector{Float64}}(undef, θ_len)
#         eigvecs = Vector{Matrix{Float64}}(undef, θ_len)
#         Λs = Vector{Matrix{Complex{Float64}}}(undef, θ_len)
#         As = Matrix{Matrix{Complex{Float64}}}(undef, θ_len, E_len)
#         data[Symbol("dat$(N)_near_bands")] = SimData_f_eigs(N, touched, θs, Es, Ts, Ss, eigvals, eigvecs, Λs, As)
#     end
#     dat = data[Symbol("dat$(N)_near_bands")]

#     @showprogress 1 "" for (i,θ) in enumerate(θs)
#         if all(dat.touched[i,:])
#             continue
#         end
#         sim = SimpleSimulationLite(N, θ=θ, δz=l, l=l, fixed_distance=l, Γ=Γ)
#         dat.eigvals[i] = sim.energies
#         dat.eigvecs[i] = sim.eigvecs
#         dat.Λs[i] = (1e-2/27.2)*sim.Lz
#         for (j,E) in enumerate(Es)
#             if any(abs.(E .- sim.energies) .< threshold) && !dat.touched[i,j]
#                 t, s, A = theorem_f(E, 1; sim=sim, rs=[0, 1, 0], γ=Γ)
#                 dat.Ts[i,j] = t
#                 dat.Ss[i,j] = s
#                 dat.As[i,j] = A
#                 dat.touched[i,j] = true
#             end
#         end
#     end
#     JLD.save("dat$(N)_near_bands.jld", "dat", dat)
#     return dat
# end

export gen_N_static_current
function gen_N_static_current(
    N;
    ϕs = range(0.8, 1.15, length=150),
    as = range(0.0, 1.0/27.2, length=5),
    Es_length = 10
)
    A = Matrix{Any}(undef, length(ϕs), length(as))
    skt = loaddir(joinpath(PARAMDIR, "ob2-1-1", "base"))

    @showprogress 1 "" for (i, ϕ) in enumerate(ϕs)
        sim = TheoryOfCISS.SimpleSimulation(N, ϕ=ϕ, l=2.0)
        z0 = (sum∘extrema)(map(x->Molecules.position(x)[3], sim.mol.atoms)) / 2
        d = - -(extrema(map(x->Molecules.position(x)[3], sim.mol.atoms))...)
        for (j, a) in enumerate(as)
            H0 = hamiltonian(Float64, skt, sim.mol, electric_field=(x, y, z) -> iszero(d) ? 0.0 : -a*(z - z0)/d)
            μ, Δ = real.(chemicalpotential(H0, sim.mol))
            sim_ = SimpleSimulation(sim.N, sim.mol, μ, Δ, H0, sim.ΓL, sim.ΓR, sim.ΔΓL, sim.ΔΓR, sim.Lx, sim.Ly, sim.Lz, sim.Lxs, sim.Lys, sim.Lzs)
            Es = sort!(eigvals(H0))

            #Es = range(E_lower, μ, length=Es_length)
            #ΔE = Es.step
            tmp = [sum_up_down(E -> theorem_magnet(E, 1, sim=sim_)[2], E0, 0.0002) for E0 in Es if E0 < μ]
            A[i, j] = (tmp, Es)
        end
    end
    return A
end

function sum_up_down(f, E0, step, threshold=1e-8)
    E = E0
    a = f(E)
    acc = 0.0
    while a > threshold
        acc += a
        E += step
        a = f(E)
    end
    Etop = E
    E = E0 - step
    a = f(E)
    while a > threshold
        acc += a
        E -= step
        a = f(E)
    end
    return acc*step
end

export gen_helicene_bandstructure
function gen_helicene_bandstructure(
    ;
    ϕs = range(0.0, 1.5, length=80),
    as = range(0.0, 5/27, length=40)
)
    A = Matrix{Any}(undef, length(ϕs), length(as))
    skt = loaddir(joinpath(PARAMDIR, "ob2-1-1", "base"))

    @showprogress 1 "" for (i, ϕ) in enumerate(ϕs)
        sim = TheoryOfCISS.HeliceneSimulation(4, ϕ=ϕ, oxidation=+1)
        z0 = (sum∘extrema)(map(x->Molecules.position(x)[3], sim.mol.atoms)) / 2
        for (j, a) in enumerate(as)
            H0 = hamiltonian(Float64, skt, sim.mol, electric_field=(x, y, z) -> -a*(z - z0))
            μ, Δ = real.(chemicalpotential(H0, sim.mol, +1))
            A[i, j] = 27.2 .* (sort!(eigvals(H0)), μ, Δ)
        end
    end
    return A
end

function integ(dat::TheoryOfCISS.SimData_f, Δ)
    di = round(Int, Δ / 2 / Float64(dat.Es.step))
    N = length(dat.Es)
    out = Matrix{Float64}(undef, length(dat.ϕs), N)
    for (i, E0) in enumerate(dat.Es)
        g = gaussian.(dat.Es, E0, Δ)'
        out[:,i] = sum((dat.Ps .* dat.Ts) .* g, dims=2) ./ sum(dat.Ts .* g, dims=2)
    end
    return out
end

function plot_integ(dat::TheoryOfCISS.SimData_f, args...)
    xs = dat.ϕs
    ys = sum(map(x -> x[1]*x[2], zip(dat.Ps, dat.Ts)), dims=2) ./ sum(dat6_hires.Ts, dims=2)
    plot(xs, ys, xlabel="ϕ", ylabel="integ. P", label="N = $(dat.N)", args...)
end

function plot_integ!(dat::TheoryOfCISS.SimData_f, args...)
    xs = dat.ϕs
    ys = sum(map(x -> x[1]*x[2], zip(dat.Ps, dat.Ts)), dims=2) ./ sum(dat6_hires.Ts, dims=2)
    plot!(xs, ys, xlabel="ϕ", ylabel="integ. P", label="N = $(dat.N)", args...)
end

using RecipesBase

@recipe function f(dat::SimData_f, g)
    xlims --> (dat.ϕs[1], dat.ϕs[end])
    ylims --> (dat.Es[1], dat.Es[end])
    data = g.(dat.Ts, dat.Ps)'
    lim = maximum(abs.(data))
    clims --> (-lim, lim)

    (dat.ϕs, dat.Es, data)
end

@recipe function f(dat::SimData_f_eigs, g)
    xlims --> (dat.ϕs[1], dat.ϕs[end])
    ylims --> (dat.Es[1], dat.Es[end])
    data = g.(dat.Ts, dat.Ps)
    lim = maximum(abs.(data))
    clims --> (-lim, lim)

    (dat.ϕs, dat.Es, data')
end

function nearest_Λ12(E, eigvals, eigvecs, Λ)
    vals = collect(enumerate(abs.(inv.(E .- eigvals))))
    sort!(vals, by=x->x[2])
    valmin = vals[end]
    valmax = vals[end-1]
    us = [eigvecs[:, valmin[1]], eigvecs[:, valmax[1]]]
    return -imag(us[1]'Λ*us[2])
end

function nearest_factor(E, eigvals, eigvecs, A, ΓR, γL)
    vals = collect(enumerate(abs.(inv.(E .- eigvals))))
    sort!(vals, by=x->x[2])
    valmin = vals[end]
    valmax = vals[end-1]
    us = [eigvecs[:, valmin[1]], eigvecs[:, valmax[1]]]

    A1 = us[1]'A*us[1]
    A2 = us[2]'A*us[2]
    A12 = us[1]'A*us[2]
    A = [A1 A12; A12 A2]

    γ11 = us[1]'γL*us[1]
    γ12 = us[1]'γL*us[2]
    γ22 = us[2]'γL*us[2]
    γ = [γ11 γ12; γ12 γ22]

    Γ11 = us[1]'ΓR*us[1]
    Γ12 = us[1]'ΓR*us[2]
    Γ22 = us[2]'ΓR*us[2]
    Γ = [Γ11 Γ12; Γ12 Γ22]

    prefactor = real(sum(A[1,i]*A[j,k]*A[3-i,2]*(γ[j,i]*Γ[3-i,k] - γ[j,3-i]*Γ[i,k]) for i in 1:2, j in 1:2, k in 1:2))

    return prefactor
end
