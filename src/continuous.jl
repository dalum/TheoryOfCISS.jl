"""
    PeriodicPolyacetylene(N; ϕ=π/4, l=DEFAULT_LENGTH, handedness=1, δz=nothing, θ=nothing, r0=nothing)

Molecule constructor for (twisted) Polyacetylene with periodic boundary conditions and twist angle, `ϕ`.

"""
@export struct PeriodicPolyacetylene <: AbstractMoleculeConstructor
    θ
    l
    handedness
    δz
    ϕ
    r0

    function PeriodicPolyacetylene(
        ;
        θ = π/2,
        l = DEFAULT_LENGTH,
        handedness = 1,
        δz = nothing,
        ϕ = nothing,
        r0 = nothing
    )
        return new(θ, l, handedness, δz, ϕ, r0)
    end
end

function (m::PeriodicPolyacetylene)(; kwargs...)
    θ = get(kwargs, :θ, m.θ)
    l = get(kwargs, :l, m.l)
    handedness = get(kwargs, :handedness, m.handedness)
    δz = get(kwargs, :δz, m.δz)
    ϕ = get(kwargs, :ϕ, m.ϕ)
    r0 = get(kwargs, :r0, m.r0)

    ρ = atan(2tan(θ/2))

    if δz !== nothing && θ !== nothing
        r0 = [l, 0, 0]
    end

    if δz === nothing
        δz = √3*l*sin(ρ)/2
    end
    if ϕ === nothing
        ϕ = π - 2atan(√3*cos(ρ))*handedness
    end

    if r0 === nothing
        r0 = [√(l^2 - δz^2) / 2sin(ϕ/2), 0, 0]
    end
    r0::Vector

    mols = [Molecule(), Molecule()]

    for i in 1:2
        rotation_matrix = LinearAlgebra.Givens(1, 2, cos((i-1)*ϕ), -sin((i-1)*ϕ))
        shift_vector = [0, 0, (i-1)*δz]
        axes = makeaxes(π/2 + (i-1)*ϕ, ρ*handedness)
        r1 = shift_vector + rotation_matrix*r0
        r2 = shift_vector + rotation_matrix*(r0 + [abs(l), 0, 0])
        A = Carbon(vec(r1...), orbital"2s", orbital"2p_x", orbital"2p_y", orbital"2p_z", axes=axes)
        B = Hydrogen(vec(r2...), orbital"1s", axes=axes)
        push!(mols[i], A, B)
        push!(mols[i], Bond(A, B))
    end

    return Molecules.UnitCell(mols..., Set([Bond(mols[1].atoms[1], mols[2].atoms[1])])), nothing, nothing
end
