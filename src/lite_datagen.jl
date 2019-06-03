mutable struct DataSlice
    metadata
    y0
    data
end
Base.show(io::IO, col::DataSlice) = print(io, "DataSlice(...)")

mutable struct DataColumn
    metadata
    ys
    data
    touched
end
DataColumn(metadata) = DataColumn(Float64[], metadata)
DataColumn(ys, metadata) = DataColumn(metadata, collect(ys), fill!(Vector(undef, length(ys)), NamedTuple()), fill(false, length(ys)))

function initialize!(col::DataColumn, ys)
    col.ys = ys
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
Base.copy(col::DataColumn) = DataColumn(col.metadata, col.ys, copy(col.data), col.touched)
Base.show(io::IO, col::DataColumn) = print(io, "DataColumn(...)")

function indexsnap(xs, x0)
    x1, x2 = x0 .- inv.(extrema(inv.(x0 .- xs)))
    y0 = abs(y2 - y0) > abs(x1 - x0) ? x1 : x2
    return findfirst(x -> x == x0, xs)
end

function find_crossing(
    x0,
    y,
    N = 14;
    l = DEFAULT_LENGTH,
    α = 1.0,
    threshold = 1e-12,
    max_iterations = 1e1,
    δx = 1e-6,
    attenuation = 1e-2
)
    function optim(x, y)
        sim = LiteSimulation(N, ϕ=x, l=l, α=α)
        eigvals, _ = nearest_eigenstates(sim.eig, y)
        d = abs(eigvals[1] - eigvals[2])
        y = (eigvals[1] + eigvals[2])/2
        return d, y
    end

    x = x0
    _, y = optim(x, y)
    for _ in 1:max_iterations
        (d, y), (d_, _) = optim(x, y), optim(x + δx, y)
        δd = d_ - d
        abs(d) < threshold && return x, y
        x = x - clamp(d*δx/δd, -attenuation, attenuation)
    end
    return x, y
end

@export function track_crossing_length(
    points;
    N0,
    pinned = true,
    N_lower = 2,
    N_upper = 52,
    N_step = 2,
    threshold = 0.001,
    α = 1.0,
    l = DEFAULT_LENGTH,
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
    xsymbol = :ϕ,
    x_len = 50,
    x_lower = 0.0,
    x_upper = π/2,
    xs = range(x_lower, x_upper, length=x_len),
    kwargs...
)
    gen_bands!(Vector{DataColumn}(xsymbol, xs), args...; xs = xs, xsymbol = xsymbol, kwargs...)
end

@export function gen_bands!(
    cols::Vector{DataColumn},
    N = 14;
    l = DEFAULT_LENGTH,
    ϕ = 1.0,
    xsymbol = :ϕ,
    x_len = 1500,
    x_lower = 0.0,
    x_upper = π/2,
    xs = range(x_lower, x_upper, length=x_len),
    α = 1.0,
    dirname = "$(now())_gen_bands($(N), alpha=$α)"
)
    mkpath(joinpath("data", dirname))

    @showprogress 1 "" for (idx,col) in enumerate(cols)
        sim = LiteSimulation(N, ϕ=ϕ, l=l, α=α, ρ0 = Diagonal([0.0, 1.0, 0.0]); (xsymbol => col.metadata.x,)...)
        col.metadata = (col.metadata..., N=N, sim = sim)
        @save joinpath("data", dirname, "col-$(string(idx, pad=6)).bson") col
    end

    return cols
end

@export function gen_on_bands!(
    cols;
    clear = false,
    dirname = "$(now())_gen_on_bands()",
    f = calc_data1
)
    mkpath(joinpath("data", dirname))

    @showprogress 1 "" for (idx,col) in enumerate(cols)
        sim = col.metadata.sim
        initialize!(col, sim.eig.values)
        for (i,y) in enumerate(col.ys)
            if !col.touched[i]
                col.touched[i] = true
                col.data[i] = f(sim, y)
            end
        end

        @save joinpath("data", dirname, "col-$(string(idx, pad=6)).bson") col
    end

    return cols
end

@export function gen_near_bands!(
    cols,
    N = 14;
    y_len = 500,
    y_lower = 0.0,
    y_upper = 2.8,
    ys = range(y_lower, y_upper, length=y_len),
    threshold = 0.1,
    α = 1.0,
    dirname = "$(now())_gen_near_bands($(N), alpha=$α)",
    f = calc_data1
)
    mkpath(joinpath("data", dirname))

    @showprogress 1 "" for (idx,col) in enumerate(cols)
        sim = col.metadata.sim
        initialize!(col, ys)
        for (i,y) in enumerate(ys)
            if any(abs.(y .- sim.eig.values) .< threshold) && !col.touched[i]
                col.touched[i] = true
                col.data[i] = f(sim, y)
            end
        end

        @save joinpath("data", dirname, "col-$(string(idx, pad=6)).bson") col
    end

    return cols
end

function calc_data1(sim::Simulation, y)
    G = propagator(sim, y)
    eigvals, U = nearest_eigenstates(sim.eig, y)
    v, t, s = polarization(G, sim.γL, sim.ΓR, sim.Lz)

    return (
        v = v,
        t = t,
        s = s,
        eigvals = eigvals,
        G = U'G*U,
        Lx = U'sim.Lx*U,
        Ly = U'sim.Ly*U,
        Lz = U'sim.Lz*U,
        γL = U'sim.γL*U,
        ΓL = U'sim.ΓL*U,
        ΓR = U'sim.ΓR*U
    )
end

function calc_data2(sim::Simulation, y)
    G = propagator(sim, y)
    eigvals, U = nearest_eigenstates(sim.eig, y)
    vx, t, _ = polarization(G, sim.γL, sim.ΓR, sim.Lx)
    vy, _, _ = polarization(G, sim.γL, sim.ΓR, sim.Ly)
    vz, _, _ = polarization(G, sim.γL, sim.ΓR, sim.Lz)

    return (
        vx = vx,
        vy = vy,
        vz = vz,
        t = t,
        eigvals = eigvals,
    )
end

function datamap(f, cols::Vector{DataColumn})
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
    y;
    calc = calc_data1,
    N = 14,
    l = DEFAULT_LENGTH,
    α = 1.0
)
    sim = LiteSimulation(N, ϕ=x, l=l, α=α)
    return f(calc(sim, y))
end

gaussian(x, x0, σ) = exp(-(x - x0)^2 / (2σ^2))
box(x, x0, σ) = abs(x - x0) < σ ? 1.0 : 0.0

@export percentage(x) = 100x

@export transmission(x::NamedTuple; λ=8e-3) = get(x, :t, 0.0)
@export spintransmission(x::NamedTuple; λ=8e-3) = λ*get(x, :s, 0.0)
@export polarization(x::NamedTuple; λ=8e-3) = λ*get(x, :s, 0.0) / get(x, :t, 1e-16)

@export Cvec(x::NamedTuple; λ=8e-3, v = :vz) = λ*2imag(get(x, v, 0.0))# * get(x, :t, 1e-16)
@export Cveclen(x::NamedTuple; λ=8e-3) = sqrt(sum(Cvec(x, v = v)^2 for v in (:vx, :vy, :vz)))
@export Cvecproj(x::NamedTuple; λ=8e-3, v = :vz) = abs(Cvec(x, v = v)) / Cveclen(x)

@export Dvec(x::NamedTuple; λ=8e-3, v = :vz) = λ*2real(get(x, v, 0.0))# * get(x, :t, 1e-16)
@export Dveclen(x::NamedTuple; λ=8e-3) = sqrt(sum(Dvec(x, v = v)^2 for v in (:vx, :vy, :vz)))
@export Dvecproj(x::NamedTuple; λ=8e-3, v = :vz) = abs(Dvec(x, v = v)) / Dveclen(x)


function transmission_twoband(x::NamedTuple; λ=8e-3)
    G = get(x, :G, fill(0.0, 2, 2))
    iszero(G) && return 0.0
    γL = get(x, :γL, fill(0.0, 2, 2))
    ΓR = get(x, :ΓR, fill(0.0, 2, 2))
    return real(dot(γL, G'ΓR*G))
end

function spintransmission_twoband(x::NamedTuple; λ=8e-3)
    G = get(x, :G, fill(0.0, 2, 2))
    iszero(G) && return 0.0
    Λ = λ*get(x, :L, fill(0.0, 2, 2))
    γL = get(x, :γL, fill(0.0, 2, 2))
    ΓR = get(x, :ΓR, fill(0.0, 2, 2))
    return 2*real(dot(γL, G'ΓR*G*Λ*G))
end

function polarization_twoband(x::NamedTuple; λ=8e-3)
    G = get(x, :G, fill(0.0, 2, 2))
    iszero(G) && return 0.0
    Λ = λ*get(x, :L, fill(0.0, 2, 2))
    γL = get(x, :γL, fill(0.0, 2, 2))
    ΓR = get(x, :ΓR, fill(0.0, 2, 2))
    X = G*γL
    Y = ΓR*G
    return 2*real(dot(X, Y*Λ*G)) / real(dot(X, Y))
end

function smoothen(cols::Vector{DataColumn}; kwargs...)
    out = Vector{DataColumn}(undef, length(cols))
    @showprogress 1 "" for (i, col) in enumerate(cols)
        out[i] = smoothen(col; kwargs...)
    end
    return out
end

function smoothen(cols::Vector{DataColumn}, y0; kwargs...)
    out = Vector{DataSlice}(undef, length(cols))
    @showprogress 1 "" for (i, col) in enumerate(cols)
        out[i] = smoothen(col, y0; kwargs...)
    end
    return out
end

function smoothen(
    col::DataColumn;
    kwargs...
)
    out = copy(col)
    for (i, y0) in enumerate(ys)
        out.data[i] = smoothen(col, y0, kwargs...).data
    end
    return out
end

function smoothen(
    col::DataColumn,
    y0;
    σ = 0.1,
    s = spintransmission,
    t = transmission,
    w = gaussian
)
    ws = w.(col.ys, y0, σ)
    slice = DataSlice(
        col.metadata,
        y0,
        (p = sum(s.(col.data) .* ws) / sum(t.(col.data) .* ws),)
    )
    return slice
end

Λ12(x::NamedTuple; λ=8e-3) = 2*imag(λ*get(x, :L, fill(0.0, 2, 2))[1,2])
function prefactor(x::NamedTuple; λ=8e-3, L = :Lz)
    G = get(x, :G, fill(0.0, 2, 2))
    iszero(G) && return 0.0, 0.0, 1e-16
    A = -im .* (G .- G') ./ 2

    Λ = λ*get(x, L, fill(0.0, 2, 2))
    γL = get(x, :γL, fill(0.0, 2, 2))
    ΓL = get(x, :ΓL, fill(0.0, 2, 2))
    ΓR = get(x, :ΓR, fill(0.0, 2, 2))

    return 2*imag(Λ[1,2]), real((A*(γL*A'ΓR - ΓR*A'γL)*A)[1,2]), real(dot(γL, A'ΓR*A))
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

function yslice(cols::Vector{DataColumn}, y0)
    idx = indexsnap(cols[1].ys, y0)
    return [DataSlice(col.metadata, y0, col.data[idx]) for col in cols]
end
