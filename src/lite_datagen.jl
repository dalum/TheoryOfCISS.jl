mutable struct DataSlice
    metadata
    y0
    data
end
Base.show(io::IO, col::DataSlice) = print(io, "DataSlice(...)")

mutable struct DataColumn
    metadata
    ys
    ysymbol
    data
    touched
end
DataColumn(metadata) = DataColumn(Float64[], metadata)
DataColumn(ys, metadata) = DataColumn(ys, nothing, metadata)
DataColumn(ys, ysymbol, metadata) = DataColumn(metadata, collect(ys), ysymbol, fill!(Vector(undef, length(ys)), NamedTuple()), fill(false, length(ys)))

function initialize!(col::DataColumn, ys, ysymbol)
    col.ys = ys
    col.ysymbol = ysymbol
    col.data = fill!(Vector(undef, length(ys)), NamedTuple())
    col.touched = fill(false, length(ys))
    return col
end

function Base.Vector{DataColumn}(xs::AbstractVector)
    self = Vector{DataColumn}(undef, length(xs))
    for (i, x) in enumerate(xs)
        self[i] = DataColumn((x=x,))
    end
    return self
end

function Base.Vector{DataColumn}(xsymbol, xs::AbstractVector)
    self = Vector{DataColumn}(undef, length(xs))
    for (i, x) in enumerate(xs)
        self[i] = DataColumn((xsymbol=xsymbol, x=x))
    end
    return self
end

Base.setindex!(col::DataColumn, value, ys...) = setindex!(col.data, map(indexsnap, ys))
Base.copy(col::DataColumn) = DataColumn(col.metadata, col.ys, col.ysymbol, copy(col.data), col.touched)
Base.show(io::IO, col::DataColumn) = print(io, "DataColumn(...)")

@export struct DataScatter end

function DataScatter(f, cols::Vector{DataColumn})
    xs = mapreduce(col -> [col.metadata.x for _ in col.ys], vcat, cols)
    ys = mapreduce(col -> col.ys, vcat, cols)
    zs = mapreduce(col -> map(f, col.data), vcat, cols)
    return xs, ys, zs
end

function indexsnap(xs, x0)
    x1, x2 = x0 .- inv.(extrema(inv.(x0 .- xs)))
    x0 = abs(x2 - x0) > abs(x1 - x0) ? x1 : x2
    return findfirst(x -> x == x0, xs)
end

@export function datamap(f, cols::Vector{DataColumn})
    ys = mapreduce(a -> collect(a.ys), hcat, cols)
    xs = mapreduce(a -> a.metadata.x, vcat, cols)
    xs = [x for _ in 1:size(ys, 1), x in xs]
    zs = mapreduce(a -> map(f, a.data), hcat, cols)
    return xs, ys, zs
end

function datamap(f, slices::Vector{DataSlice})
    xs = mapreduce(a -> a.metadata.x, vcat, slices)
    ys = mapreduce(a -> f(a.data), vcat, slices)
    return xs, ys
end

function energy_bands(cols::Vector{DataColumn}; step=1)
    eigvalues = mapreduce(col -> col.metadata.sim.eig.values, hcat, cols)'
    energy_range = range(floor(Int, minimum(eigvalues)), ceil(Int, maximum(eigvalues)), step=step)
    return eigvalues, energy_range
end

function energy_bounds(cols::Vector{DataColumn})
    gmin, gmax = 0.0, 0.0
    for col in cols
        gmin, gmax = energy_bounds(col)
    end
    return gmin, gmax
end

function energy_bounds(col::DataColumn)
    gmin = min(gmin, floor(Int, minimum(col.metadata.sim.eig.values)))
    gmax = max(gmax, ceil(Int, maximum(col.metadata.sim.eig.values)))
    return gmin, gmax
end

function ranges(cols::Vector{DataColumn})
    xs = mapreduce(a -> a.metadata.x, vcat, cols)
    ys = collect(cols[1].ys)
    return xs, ys
end

function Base.range(slices::Vector{DataSlice})
    return mapreduce(a -> a.metadata.x, vcat, slices)
end

function atpoint(
    f,
    x,
    y,
    args...;
    xsymbol,
    ysymbol,
    calc = calc_data1,
    v = Dict(),
    kwargs...
)
    sim = LiteSimulation(args...; (xsymbol => x,)..., kwargs...)
    return f(calc(sim; v..., (ysymbol => y,)...))
end

gaussian(x, x0, σ) = exp(-(x - x0)^2 / (2σ^2))
box(x, x0, σ) = abs(x - x0) < σ ? 1.0 : 0.0

function find_crossing(
    x0,
    y0,
    args...;
    xsymbol,
    threshold = 1e-12,
    max_iterations = 1e1,
    δx = 1e-6,
    attenuation = 1e-2,
    rounding = :average,
    kwargs...
)
    function optim(x, y)
        sim = LiteSimulation(args...; kwargs..., (xsymbol => x,)...)
        eigvals, _ = nearest_eigenstates(sim.eig, y, n=2)
        d = abs(eigvals[1] - eigvals[2])
        y = (eigvals[1] + eigvals[2])/2
        return d, y
    end

    x = x0
    _, y = optim(x, y0)
    for _ in 1:max_iterations
        (d, y), (d_, _) = optim(x, y), optim(x + δx, y)
        δd = d_ - d
        abs(d) < threshold && return x, y
        x = x - clamp(d*δx/δd, -attenuation, attenuation)
    end

    if rounding in (:upper, :lower, :both)
        sim = LiteSimulation(args...; kwargs..., (xsymbol => x,)...)
        eigvals, _ = nearest_eigenstates(sim.eig, y, n=2)
        rounding == :lower && return (x, minimum(eigvals))
        rounding == :upper && return (x, maximum(eigvals))
        rounding == :both && return (x, minimum(eigvals), maximum(eigvals))
    end

    return x, y
end

@export function track_crossing_length(
    points,
    args...;
    N0,
    pinned = true,
    N_lower = 2,
    N_upper = 52,
    N_step = 2,
    threshold = 0.001,
    f = prefactor,
    kwargs...
)
    Ns = vcat(N0:-N_step:N_lower, N0:N_step:N_upper)
    Ns_ordered = N_lower:N_step:N_upper
    indices = indexin(Ns, Ns_ordered)
    dat = fill!(Matrix{NamedTuple}(undef, length(Ns_ordered), length(points)), NamedTuple())
    points0 = copy(points)
    @showprogress 1 "" for (i, N) in zip(indices, Ns), (j, (x,y)) in enumerate(points)
        if N == N0 && !pinned
            x, y = points0[j]
        end
        x, y = find_crossing(x, y, N, threshold=threshold)
        if !pinned
            points[j] = (x, y)
        end

        a, b, c = atpoint(f, x, y; N = N, kwargs...)
        dat[i,j] = (x = x, y = y, p = a*b/c)
    end
    return dat
end

@export function load_columns(dirname)
    list = readdir(joinpath("data", dirname))
    cols = Vector{DataColumn}(undef, length(list))
    @showprogress 1 "" for (i, filename) in enumerate(list)
        @load joinpath("data", dirname, filename) col
        cols[i] = col
    end
    return cols
end

@export function gen_bands(
    args...;
    xsymbol,
    bounds,
    nsamples = 50,
    xs = range(bounds..., length=nsamples),
    kwargs...
)
    gen_bands!(Vector{DataColumn}(xsymbol, xs), args...; xs = xs, xsymbol = xsymbol, kwargs...)
end

function gen_bands!(
    cols::Vector{DataColumn},
    args...;
    xsymbol,
    xs,
    save = false,
    dirname = "$(now())_gen_bands()",
    kwargs...
)
    save && mkpath(joinpath("data", dirname))

    @showprogress 1 "" for (idx,col) in enumerate(cols)
        sim = LiteSimulation(args...; kwargs..., (xsymbol => col.metadata.x,)...)
        col.metadata = (col.metadata..., sim = sim, xlims=extrema(xs))
        save && @save joinpath("data", dirname, "col-$(string(idx, pad=6)).bson") col
    end

    return cols
end

@export function gen_on_bands!(
    cols;
    bounds = energy_bounds(cols),
    save = false,
    dirname = "$(now())_gen_on_bands()",
    kwargs...
)
    save && mkpath(joinpath("data", dirname))

    @showprogress 1 "" for (idx,col) in enumerate(cols)
        ys = filter(y -> bounds[1] <= y <= bounds[2], col.metadata.sim.eig.values)
        gen_col!(col; bounds=bounds, kwargs..., ys = ys, _progress_dt=Inf)
        save && @save joinpath("data", dirname, "col-$(string(idx, pad=6)).bson") col
    end

    return cols
end

@export function gen_near_bands!(
    cols;
    save = false,
    dirname = "$(now())_gen_near_bands()",
    kwargs...
)
    save && mkpath(joinpath("data", dirname))

    @showprogress 1 "" for (idx, col) in enumerate(cols)
        gen_col!(col; kwargs..., _progress_dt=Inf)
        save && @save joinpath("data", dirname, "col-$(string(idx, pad=6)).bson") col
    end

    return cols
end

@export function gen_col!(
    col::DataColumn;
    ysymbol,
    bounds = energy_bounds(col),
    nsamples = 50,
    ys = range(bounds..., length=nsamples),
    threshold = nothing,
    f = calc_data1,
    cond = (sim, y) -> isnothing(threshold) ? true : any(abs.(y .- sim.eig.values) .< threshold),
    clear = true,
    _progress_dt = 1,
    kwargs...

)
    col.metadata = (col.metadata..., ylims=extrema(ys))
    sim = col.metadata.sim
    if clear || col.ys != ys || col.ysymbol != ysymbol
        initialize!(col, ys, ysymbol)
    end
    @showprogress _progress_dt for (i, y) in enumerate(ys)
        if cond(sim, y) && !col.touched[i]
            col.touched[i] = true
            col.data[i] = f(sim; kwargs..., (ysymbol => y,)...)
        end
    end
    return col
end

function first_order_perturbation(sim::Simulation; E=sim.μ, α=DEFAULT_α, λ=DEFAULT_λ)
    γL = α*sim.γL
    ΓL = α*sim.ΓL
    ΓR = α*sim.ΓR

    G = propagator(sim, E=E, Γ=ΓL+ΓR)
    X = G*γL*G'
    t = 2real(dot(X, ΓR))
    s = map(L -> 4λ*dot(X, ΓR*G*L), sim.L)

    return (E=E, α=α, s=s, t=t)
end

function calc_data1(sim::Simulation; E=sim.μ, α=DEFAULT_α, λ=DEFAULT_λ, n=2)
    G = propagator(sim, E=E)
    eigvals, U = nearest_eigenstates(sim.eig, E, n=n)

    X = G*sim.γL*G'
    t = 2real(dot(X, sim.ΓR))
    s = map(L -> 4λ*dot(X, sim.ΓR*G*L), sim.L)

    return (
        E=E, α=α, n=n,

        s = s,
        t = t,
        eigvals = eigvals,
        G = U'G*U,
        γL = U'sim.γL*U,
        ΓL = U'sim.ΓL*U,
        ΓR = U'sim.ΓR*U,
        L = map(Li -> U'Li*U, sim.L)
    )
end

function calc_data2(sim::Simulation; E=sim.μ, α=DEFAULT_α, λ=DEFAULT_λ)
    G = propagator(sim, E=E)
    carbon_atoms = collect(filter(x -> x isa Carbon, sim.mol.atoms))

    Ls = collect(zip(
        (angularmomentum(:x, sim.mol, atom) for atom in carbon_atoms),
        (angularmomentum(:y, sim.mol, atom) for atom in carbon_atoms),
        (angularmomentum(:z, sim.mol, atom) for atom in carbon_atoms)
    ))

    eigvals, U = nearest_eigenstates(sim.eig, E)

    ss = NTuple{3,Complex{Float64}}[]
    ts = Float64[]

    for L in Ls
        s, t = calculate_s_and_t(G, sim.γL, sim.ΓR, L, λ=λ)
        push!(ss, v)
        push!(ts, t)
    end

    return (
        E=E, α=α, n=n,

        ss = ss,
        ts = ts,
        eigvals = eigvals,
        G = U'G*U,
        γL = U'sim.γL*U,
        ΓL = U'sim.ΓL*U,
        ΓR = U'sim.ΓR*U,
        Ls = [map(Li -> U'Li*U, L) for L in Ls]
    )
end

# Calculate the reduced

function calc_data3(sim::Simulation; E=sim.μ, α=DEFAULT_α, λ=DEFAULT_λ, n=2, nmin=1)
    G = propagator(sim, E=E)
    Gfull = fullpropagator(sim, E=E, λ=λ)
    eigvals, U = nearest_eigenstates(sim.eig, E, n=n, nmin=nmin)
    Ufull = U ⊗ σ0

    s, t = calculate_s_and_t(G, sim.γL, sim.ΓR, sim.L, λ=λ)

    return (
        E=E, α=α, λ=λ, n=n,

        s = s,
        t = t,
        eigvals = eigvals,

        G = U'G*U,
        Gfull = Ufull'Gfull*Ufull,
        γL = U'sim.γL*U,
        ΓL = U'sim.ΓL*U,
        ΓR = U'sim.ΓR*U,
        L = map(Li -> U'Li*U, sim.L),
    )
end

function calc_data4(sim::Simulation; E=sim.μ, α=DEFAULT_α, λ=DEFAULT_λ, n=2)
    G = propagator(sim, E=E)
    eigvals, U, V = nearest_eigenstates(sim.eig, E, n=n)

    s, t = calculate_s_and_t(G, sim.γL, sim.ΓR, sim.L, λ=λ)

    return (
        E=E, α=α, λ=λ, n=n,

        s = s,
        t = t,
        eigvals = eigvals,
        G = G,
        U = U,
        V = V,
        γL = sim.γL,
        ΓL = sim.ΓL,
        ΓR = sim.ΓR,
        L = sim.L
    )
end

function calc_data5(sim::Simulation; E=sim.μ, α=DEFAULT_α, λ=DEFAULT_λ, n=2, Γτ=1.0, γL=sim.γL)
    G = propagator(sim, E=E, ΓL=I / (Γτ*α*λ))
    eigvals, U = nearest_eigenstates(sim.eig, E, n=n)
    s, t = calculate_s_and_t(G, γL ./ (2π*Γτ), sim.ΓR, sim.L, λ=λ)

    return (
        E=E, α=α, n=n,

        s = s,
        t = t,
        eigvals = eigvals,
        G = U'G*U,
        γL = U'sim.γL*U,
        ΓL = U'sim.ΓL*U,
        ΓR = U'sim.ΓR*U,
        L = map(Li -> U'Li*U, sim.L)
    )
end

function calc_data6(sim::Simulation; E=sim.μ, α=DEFAULT_α, λ=DEFAULT_λ, c=4, cL=1.0)
    Γ = α/2 * sim.ΓR ⊗ σ0
    ΔΓ = α/2 * sim.ΓR * (c - 1) / (c + 1) ⊗ σ[:z]
    γL = sim.γL ⊗ σ0
    ΓL = α * sim.ΓL ⊗ σ0

    Gup = fullpropagator(
        sim,
        E = E,
        λ = λ,
        Γ = Γ + ΔΓ+ cL*ΓL
    )
    Gdown = fullpropagator(
        sim,
        E = E,
        λ = λ,
        Γ = Γ - ΔΓ + cL*ΓL
    )

    t_up = tr(γL*Gup'ΓL*Gup)
    t_down = tr(γL*Gdown'ΓL*Gdown)

    return (
        E=E, α=α, λ=λ,
        t_up = t_up,
        t_down = t_down,
        spintransmission = (t_up - t_down) / 2,
        transmission = (t_up + t_down) / 2,
    )
end

function simple_transport(sim::Simulation; E=sim.μ, α=DEFAULT_α, λ=DEFAULT_λ)
    ΓL = sim.ΓL ⊗ σ0
    ΓR = sim.ΓR ⊗ σ0

    G = fullpropagator(sim, E=E, λ=λ, Γ=ΓL+ΓR)
    transmission = tr(ΓL*G'*ΓR*G)

    return (E=E, α=α, λ=λ, transmission=transmission)
end

@export percentage(x) = 100x

@export function transmission(x::NamedTuple; default=0.0)
    return real(get(x, :transmission, get(x, :t, default)))
end

@export function spintransmission(x::NamedTuple)
    return real.(get(x, :spintransmission, get(x, :s, [0.0, 0.0, 0.0])))
end
spintransmission(i::Int; kwargs...) = (x::NamedTuple) -> spintransmission(x, i; kwargs...)
spintransmission(x::NamedTuple, i; kwargs...) = spintransmission(x; kwargs...)[i]

@export polarization(x::NamedTuple; kwargs...) = spintransmission(x; kwargs...) / transmission(x; default=1e-16, kwargs...)
polarization(x::NamedTuple, i; kwargs...) = spintransmission(x, i; kwargs...) / transmission(x; default=1e-16, kwargs...)
polarization(i::Int; kwargs...) = (x::NamedTuple) -> polarization(x, i; kwargs...)

@export function precession(x::NamedTuple)
    return imag.(get(x, :s, [0.0, 0.0, 0.0])) / transmission(x; default=1e-16)
end
precession(x::NamedTuple, i; kwargs...) = precession(x; kwargs...)[i]
precession(i::Int; kwargs...) = (x::NamedTuple) -> precession(x, i; kwargs...)

@export function transmission_reduced(x::NamedTuple; α=get(x, :α, DEFAULT_α), λ=get(x, :λ, DEFAULT_λ))
    G = get(x, :G, nothing)
    isnothing(G) && return 0.0

    γL, ΓR = (α*λ)*x[:γL], (α*λ)*x[:ΓR]
    return 2real(dot(γL, G'ΓR*G))
end

@export function spintransmission_reduced(x::NamedTuple, i=nothing; α=get(x, :α, DEFAULT_α), λ=get(x, :λ, DEFAULT_λ))
    G = get(x, :G, nothing)
    isnothing(G) && return (i == nothing ? [0.0, 0.0, 0.0] : 0.0)
    Ls, γL, ΓR = x[:L], (α*λ)*x[:γL], (α*λ)*x[:ΓR]

    X = G*γL
    Y = ΓR*G
    s = map(L -> 4real(dot(X, Y*L*G)), Ls)
    return (i == nothing ? s : s[i]) * λ
end
spintransmission_reduced(i::Int; kwargs...) = (x::NamedTuple) -> spintransmission_reduced(x, i; kwargs...)

@export function polarization_reduced(x::NamedTuple, i=nothing; α=get(x, :α, DEFAULT_α), λ=get(x, :λ, DEFAULT_λ))
    G = get(x, :G, nothing)
    isnothing(G) && return (i == nothing ? [0.0, 0.0, 0.0] : 0.0)
    Ls, γL, ΓR = x[:L], (α*λ)*x[:γL], (α*λ)*x[:ΓR]

    X = G*γL
    Y = ΓR*G
    s = map(L -> 4real(dot(X, Y*L*G)), Ls)
    t = 2real(dot(X, Y))
    return (i == nothing ? s : s[i]) * λ / t
end
polarization_reduced(i::Int; kwargs...) = (x::NamedTuple) -> polarization_reduced(x, i; kwargs...)

@export function precession_reduced(x::NamedTuple, i=nothing; α=get(x, :α, DEFAULT_α), λ=get(x, :λ, DEFAULT_λ))
    G = get(x, :G, nothing)
    isnothing(G) && return (i == nothing ? [0.0, 0.0, 0.0] : 0.0)
    Ls, γL, ΓR = x[:L], (α*λ)*x[:γL], (α*λ)*x[:ΓR]

    X = G*γL
    Y = ΓR*G
    s = map(L -> 4imag(dot(X, Y*L*G)), Ls)
    t = 2real(dot(X, Y))
    return (i == nothing ? s : s[i]) * λ / t
end
precession_reduced(i::Int; kwargs...) = (x::NamedTuple) -> precession_reduced(x, i; kwargs...)

@export function fulltransmission_reduced(x::NamedTuple; α=get(x, :α, DEFAULT_α), λ=get(x, :λ, DEFAULT_λ))
    G = get(x, :Gfull, fill(0.0, 4, 4))
    iszero(G) && return 0.0
    γL, ΓR = (α*λ)*x[:γL] ⊗ σ0, (α*λ)*x[:ΓR] ⊗ σ0

    return real(dot(γL, G'ΓR*G))
end

@export function fullspintransmission_reduced(x::NamedTuple, i=nothing; α=get(x, :α, DEFAULT_α), λ=get(x, :λ, DEFAULT_λ))
    G = get(x, :Gfull, nothing)
    isnothing(G) && return 0.0
    γL, ΓR = (α*λ)*x[:γL] ⊗ σ0, (α*λ)*x[:ΓR] ⊗ σ0
    M = Matrix(1.0I, size(G, 1) ÷ 2, size(G, 2) ÷ 2)

    s = map(σi -> real(dot(γL, G'ΓR*(M ⊗ σi)*G)), collect(σ))
    return (i == nothing ? s : s[i])
end
fullspintransmission_reduced(i::Int; kwargs...) = (x::NamedTuple) -> fullspintransmission_reduced(x, i; kwargs...)

@export function fullpolarization_reduced(x::NamedTuple, i=nothing; α=get(x, :α, DEFAULT_α), λ=get(x, :λ, DEFAULT_λ))
    G = get(x, :Gfull, nothing)
    isnothing(G) && return 0.0
    γL, ΓR = (α*λ)*x[:γL] ⊗ σ0, (α*λ)*x[:ΓR] ⊗ σ0
    M = Matrix(1.0I, size(G, 1) ÷ 2, size(G, 2) ÷ 2)

    s = map(σi -> real(dot(γL, G'ΓR*(M ⊗ σi)*G)), collect(σ))
    t = real(dot(γL, G'ΓR*G))
    return (i == nothing ? s : s[i]) / t
end
fullpolarization_reduced(i::Int; kwargs...) = (x::NamedTuple) -> fullpolarization_reduced(x, i; kwargs...)

@export function smoothen(f, cols::Vector{DataColumn}; kwargs...)
    out = Vector{DataColumn}(undef, length(cols))
    @showprogress 1 "" for (i, col) in enumerate(cols)
        out[i] = smoothen(f, col; kwargs...)
    end
    return out
end

function smoothen(f, cols::Vector{DataColumn}, y0; kwargs...)
    out = Vector{DataSlice}(undef, length(cols))
    @showprogress 1 "" for (i, col) in enumerate(cols)
        out[i] = smoothen(f, col, y0; kwargs...)
    end
    return out
end

function smoothen(f, col::DataColumn; kwargs...)
    out = copy(col)
    for (i, y0) in enumerate(col.ys)
        out.data[i] = smoothen(f, col, y0; kwargs...).data
    end
    return out
end

function smoothen(fs::NamedTuple, col::DataColumn, y0; σ = 0.1, w = gaussian)
    ws = w.(col.ys, y0, σ)
    slice = DataSlice(col.metadata, y0, map(f -> sum(x -> f(x[1]) * x[2], zip(col.data, ws)), fs))
    return slice
end

function smoothen(f, col::DataColumn, y0; σ = 0.1, w = gaussian)
    ws = w.(col.ys, y0, σ)
    slice = DataSlice(col.metadata, y0, (result = sum(f.(col.data) .* ws),))
    return slice
end

Λ12(x::NamedTuple; λ=8e-3) = 2*imag(λ*get(x, :L, fill(0.0, 2, 2))[1,2])
function prefactor(x::NamedTuple; i=nothing, α=get(x, :α, DEFAULT_α), λ=get(x, :λ, DEFAULT_λ))
    G = get(x, :G, nothing)
    isnothing(G) && return (i == nothing ? [0.0, 0.0, 0.0] : 0.0)
    Ls, γL, ΓR = x[:L], (α*λ)*x[:γL], (α*λ)*x[:ΓR]

    Λ12 = λ .* map(L -> L[1,2], Ls)
    s = 2 * (G'*(γL*G*ΓR - ΓR*G*γL)*G')[2,1]
    t = real(tr(γL*G'ΓR*G))
    return imag(i == nothing ? Λ12 : Λ12[i]), s / t

end

function polarization_prefactor(x::NamedTuple)
    a, b, c = prefactor(x)
    return a*b/c
end

function annotate_crossings!(f::Function, points; N=14, fontsize=6, xtransform=identity, ztransform=identity, kwargs...)
    for (x,y) in points
        x, y = find_crossing(x, y, N, threshold=1e-6)
        a, b, c = atpoint(f, x, y; N = N, kwargs...)
        a, b = sign(b)*a, ztransform(abs(b))
        data = "\\($(ScientificNotation.string(a / 1e-3, backend=:latex, digits=1, lower=-2, upper=2)), $(ScientificNotation.string(b/c * 1e-3, backend=:latex, digits=1, lower=-2, upper=2))\\)"
        Plots.annotate!([(xtransform(x), y, Plots.text(data, fontsize, :center))])
    end
    Plots.plot!()
end

@export function yslice(cols::Vector{DataColumn}, y0)
    idx = indexsnap(cols[1].ys, y0)
    return [DataSlice(col.metadata, y0, col.data[idx]) for col in cols]
end
