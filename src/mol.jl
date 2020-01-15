@export ϕtoθ(ϕ) = π - 2atan(√3 * cos(ϕ))
@export ϕtoδz(ϕ, l=DEFAULT_LENGTH) = √3/2 * l * sin(ϕ)

abstract type AbstractMoleculeConstructor end

"""
    Cumulene(N; θ=π/2*(N % 2), l=DEFAULT_LENGTH, orbitals=(orbital"2s", orbital"2p_x", orbital"2p_y", orbital"2p_z"))

Molecule constructor for Cumulene with relative twist angle, `θ`, of the
hydrogen positions at either end.  The orbitals of the Carbon atoms
can be restricted by providing an iterable of `orbitals` as a keyword
argument.

"""
@export struct Cumulene <: AbstractMoleculeConstructor
    N
    θ
    l_ch
    l_cc
    orbitals
    valence
    left_index::Int
    right_index::Int

    function Cumulene(
        ;
        N = 5,
        θ = π/2*(N % 2),
        l_ch = CH_BOND_LENGTH,
        l_cc = CC_TRIPLE_BOND_LENGTH,
        orbitals = (orbital"2s", orbital"2p_x", orbital"2p_y", orbital"2p_z"),
        valence = missing,
        left_index = 3,
        right_index = 3
    )
        return new(N, θ, l_ch, l_cc, orbitals, valence, left_index, right_index)
    end
end

function (m::Cumulene)(; kwargs...)
    N = get(kwargs, :N, m.N)
    θ = get(kwargs, :θ, m.θ)
    l_ch = get(kwargs, :l_ch, m.l_ch)
    l_cc = get(kwargs, :l_cc, m.l_cc)
    orbitals = get(kwargs, :orbitals, m.orbitals)
    valence = get(kwargs, :valence, m.valence)

    mol = Molecule()

    r0 = [0.0, 0.0, 0.0]
    r1 = [-√3/2, 0, -0.5] * l_ch
    r2 = [√3/2, 0, -0.5] * l_ch
    axes = makeaxes(0.0, 0.0)

    A = Hydrogen(r1, axes=axes)
    B = Hydrogen(r2, axes=axes)

    C = Carbon(r0, orbitals..., axes=axes, valence=valence)

    push!(mol, A, B, C)
    push!(mol, Bond(A, C), Bond(B, C))

    C′ = C
    for i in 1:(N-1)
        axes = makeaxes(N > 1 ? θ*i/(N-1) : 0.0, 0.0)
        r = [0, 0, i*l_cc]
        C = Carbon(r, orbitals..., axes=axes, valence=valence)
        push!(mol, C)
        i > 0 && push!(mol, Bond(C, C′))
        C′ = C
    end

    r1 = [-√3/2 * cos(θ), -√3/2 * sin(θ), (N - 0.5)] * l_ch
    r2 = [√3/2 * cos(θ), √3/2 * sin(θ), (N - 0.5)] * l_ch
    axes = makeaxes(θ, 0.0)

    A = Hydrogen(r1, axes=axes)
    B = Hydrogen(r2, axes=axes)

    push!(mol, A, B)
    push!(mol, Bond(A, C′), Bond(B, C′))
    return mol
end

"""
    Polyacetylene(N; ϕ=π/4, l=DEFAULT_LENGTH, handedness=1, δz=nothing, θ=nothing, r0=nothing)

Molecule constructor for (twisted) Polyacetylene with twist angle, `ϕ`.

"""
@export struct Polyacetylene <: AbstractMoleculeConstructor
    N
    θ
    l_ch
    l_cc
    handedness
    δz
    ϕ
    r0
    left_index::Int
    right_index::Int

    function Polyacetylene(
        ;
        N = 12,
        θ = π/2,
        l_ch = CH_BOND_LENGTH,
        l_cc = (CC_SINGLE_BOND_LENGTH + CC_DOUBLE_BOND_LENGTH)/2,
        handedness = 1,
        δz = nothing,
        ϕ = nothing,
        r0 = nothing,
        left_index = 2,
        right_index = 3
    )
        return new(N, θ, l_ch, l_cc, handedness, δz, ϕ, r0, left_index, right_index)
    end
end

function (m::Polyacetylene)(; kwargs...)
    N = get(kwargs, :N, m.N)
    θ = get(kwargs, :θ, m.θ)
    l_ch = get(kwargs, :l_ch, m.l_ch)
    l_cc = get(kwargs, :l_cc, m.l_cc)
    handedness = get(kwargs, :handedness, m.handedness)
    δz = get(kwargs, :δz, m.δz)
    ϕ = get(kwargs, :ϕ, m.ϕ)
    r0 = get(kwargs, :r0, m.r0)

    ρ = atan(2tan(θ/2)) + π * (θ > π)

    if δz !== nothing && θ !== nothing
        r0 = [l_cc, 0, 0]
    end

    if δz === nothing
        δz = √3*l_cc*sin(ρ)/2
    end
    if ϕ === nothing
        ϕ = π - 2atan(√3*cos(ρ))*handedness
    end

    if r0 === nothing
        r0 = [√(l_cc^2 - δz^2) / 2sin(ϕ/2), 0, 0]
    end
    r0::Vector

    mol = Molecule()

    rotation_matrix = LinearAlgebra.Givens(1, 2, cos(ϕ), -sin(ϕ))
    shift_vector = [0, 0, δz]
    r1 = shift_vector + rotation_matrix*r0

    axes = makeaxes(π/2, ρ*handedness)
    A′ = Hydrogen(r1 + (r0 - r1)*l_ch/l_cc, axes=axes)
    push!(mol, A′)

    for i in 1:N
        rotation_matrix = LinearAlgebra.Givens(1, 2, cos(i*ϕ), -sin(i*ϕ))
        shift_vector = [0, 0, i*δz]
        axes = makeaxes(π/2 + i*ϕ, ρ*handedness)
        r1 = shift_vector + rotation_matrix*r0
        r2 = shift_vector + rotation_matrix*(r0 + [l_ch, 0, 0])
        A = Carbon(r1, axes=axes)
        B = Hydrogen(r2, axes=axes)
        push!(mol, A, B)
        push!(mol, Bond(A, A′), Bond(A, B))
        A′ = A
    end
    rotation_matrix = LinearAlgebra.Givens(1, 2, cos((N+1)*ϕ), -sin((N+1)*ϕ))
    shift_vector = [0, 0, (N+1)*δz]
    axes = makeaxes(π/2 + (N+1)*ϕ, ρ*handedness)
    r3 = shift_vector + rotation_matrix*r0
    A = Hydrogen(r1 + (r3 - r1)*l_ch/l_cc, axes=axes)
    push!(mol, A)
    push!(mol, Bond(A, A′))
    return mol
end

###########################################
###########################################
###########################################

"""
    Helicene(N; ϕ=π/4, l=DEFAULT_LENGTH, handedness=1, δz=nothing, θ=nothing, r0=nothing)

Molecule constructor for Helicene with twist angle, `ϕ`.

"""
@export struct Helicene <: AbstractMoleculeConstructor
    N
    l_ch
    l_cc
    δz
    left_index::Int
    right_index::Int

    function Helicene(
        ;
        N = 4,
        l_ch = CH_BOND_LENGTH,
        l_cc = (CC_SINGLE_BOND_LENGTH + CC_DOUBLE_BOND_LENGTH)/2,
        δz = 0.5,
        left_index = 4,
        right_index = 3
    )
        return new(N, l_ch, l_cc, δz, left_index, right_index)
    end
end

function (m::Helicene)(; kwargs...)
    N = get(kwargs, :N, m.N)
    δz = get(kwargs, :δz, m.δz)
    l_ch = get(kwargs, :l_ch, m.l_ch)
    l_cc = get(kwargs, :l_cc, m.l_cc)

    shift_vector = [0, 0, δz]
    θ = 2π/6
    ϕ = 0

    r0 = [l_cc, 0, 0]
    r_12 = [l_cc, 0, 0]
    r_23 = [l_cc/2, √3/2*l_cc, 0]
    r_34 = [-l_cc/2, √3/2*l_cc, 0]
    r_35 = [l_ch, 0, 0]
    r_46 = [l_ch/2, √3/2*l_ch, 0]

    r1 = r0
    r2 = r1 + r_12
    r3 = r2 + r_23
    r4 = r3 + r_34
    r5 = r3 + r_35
    r6 = r4 + r_46

    ϕ1 = atan(r1[2], r1[1])
    ϕ2 = atan(r2[2], r2[1])
    ϕ3 = atan(r3[2], r3[1])
    ϕ4 = atan(r4[2], r4[1])
    ϕ5 = atan(r5[2], r5[1])
    ϕ6 = atan(r6[2], r6[1])

    mol = Molecule()

    i = -1
    axes = makeaxes(π/2 + i*θ, 0)
    rotation_matrix = LinearAlgebra.Givens(1, 2, cos(i*θ + ϕ), -sin(i*θ + ϕ))

    r1_ = (i + 3ϕ1/π)*shift_vector + rotation_matrix*r1
    r4_ = (i + 3ϕ4/π)*shift_vector + rotation_matrix*r4
    A′ = Hydrogen(r1 + (r1_ - r1)*l_ch/l_cc, axes=axes)
    D′ = Hydrogen(r2 + (r4_ - r2)*l_ch/l_cc, axes=axes)
    push!(mol, A′)
    push!(mol, D′)
    for i in 0:(N-1)
        rotation_matrix = LinearAlgebra.Givens(1, 2, cos(i*θ + ϕ), -sin(i*θ + ϕ))
        axes = makeaxes(π/2 + i*θ, 0)

        r1_ = (i + 3ϕ1/π)*shift_vector + rotation_matrix*r1
        r2_ = (i + 3ϕ2/π)*shift_vector + rotation_matrix*r2
        r3_ = (i + 3ϕ3/π)*shift_vector + rotation_matrix*r3
        r4_ = (i + 3ϕ4/π)*shift_vector + rotation_matrix*r4
        r5_ = (i + 3ϕ5/π)*shift_vector + rotation_matrix*r5
        r6_ = (i + 3ϕ6/π)*shift_vector + rotation_matrix*r6

        A = Carbon(r1_, axes=axes)
        B = Carbon(r2_, axes=axes)
        C = Carbon(r3_, axes=axes)
        D = Carbon(r4_, axes=axes)
        E = Hydrogen(r5_, axes=axes)
        F = Hydrogen(r6_, axes=axes)
        push!(mol, A, B, C, D, E, F)
        push!(mol, Bond(A, A′), Bond(A, B), Bond(B, C), Bond(C, D), Bond(B, D′), Bond(C, E), Bond(D, F))
        A′ = A
        D′ = D
    end
    i = N
    rotation_matrix = LinearAlgebra.Givens(1, 2, cos(i*θ + ϕ), -sin(i*θ + ϕ))
    axes = makeaxes(π/2 + i*θ, 0)

    r1_ = (i + 3ϕ1/π)*shift_vector + rotation_matrix*r1
    r2_ = (i + 3ϕ2/π)*shift_vector + rotation_matrix*r2
    r3_ = (i + 3ϕ3/π)*shift_vector + rotation_matrix*r3

    A = Carbon(r1_, axes=axes)
    B = Carbon(r2_, axes=axes)
    C = Hydrogen(r2_ + (r3_ - r2_)*l_ch/l_cc, axes=axes)

    i = N + 1
    axes = makeaxes(π/2 + i*θ, 0)

    rotation_matrix = LinearAlgebra.Givens(1, 2, cos(i*θ + ϕ), -sin(i*θ + ϕ))
    r7_ = (i + 3ϕ1/π)*shift_vector + rotation_matrix*r1

    D = Hydrogen(r1_ + (r7_ - r1_)*l_ch/l_cc, axes=axes)

    push!(mol, A, B, C, D)
    push!(mol, Bond(A, A′), Bond(A, B), Bond(B, C), Bond(B, D′), Bond(A, D))

    return mol
end

function makehelicene2(N, a=0.5, l=1.4; handedness=1)
    shift_vector = [0, 0, a]
    θ = 2π/6#π - 2atan(√3 * cos(ϕ))*handedness
    ϕ1 = 0
    ϕ2 = 0
    ϕ3 = atan(√3, 5)
    ϕ4 = atan(√3, 2)
    ϕ5 = atan(√3, 7)
    ϕ6 = atan(3*√3, 5)
    r0 = [l, 0, 0]

    mol = Molecule()

    axes = makeaxes(π/2, 0.0)
    i = 0
    rotation_matrix1 = LinearAlgebra.Givens(1, 2, cos(i*θ + ϕ1), -sin(i*θ + ϕ1))
    rotation_matrix4 = LinearAlgebra.Givens(1, 2, cos(i*θ + ϕ4), -sin(i*θ + ϕ4))
    rotation_matrix5 = LinearAlgebra.Givens(1, 2, cos((i-1)*θ + ϕ1), -sin((i-1)*θ + ϕ1))

    r1 = (i + 3ϕ1/π)*shift_vector + rotation_matrix1*r0
    r4 = (i + 3ϕ4/π)*shift_vector + rotation_matrix4*(√(5 + √3)*r0)
    r5 = ((i-1) + 3ϕ1/π)*shift_vector + rotation_matrix5*r0
    A′ = Oxygen(vec(r1...), orbital"2s", orbital"2p_x", orbital"2p_y", orbital"2p_z", axes=axes)
    D′ = Hydrogen(vec(r4...), orbital"1s", axes=axes)
    F′ = Hydrogen(vec(r5...), orbital"1s", axes=axes)
    push!(mol, A′, D′, F′)
    push!(mol, Bond(A′, F′))

    for i in 1:N
        rotation_matrix1 = LinearAlgebra.Givens(1, 2, cos(i*θ + ϕ1), -sin(i*θ + ϕ1))
        rotation_matrix2 = LinearAlgebra.Givens(1, 2, cos(i*θ + ϕ2), -sin(i*θ + ϕ2))
        rotation_matrix3 = LinearAlgebra.Givens(1, 2, cos(i*θ + ϕ3), -sin(i*θ + ϕ3))
        rotation_matrix4 = LinearAlgebra.Givens(1, 2, cos(i*θ + ϕ4), -sin(i*θ + ϕ4))
        rotation_matrix5 = LinearAlgebra.Givens(1, 2, cos(i*θ + ϕ5), -sin(i*θ + ϕ5))
        rotation_matrix6 = LinearAlgebra.Givens(1, 2, cos(i*θ + ϕ6), -sin(i*θ + ϕ6))
        axes = makeaxes(π/2 + i*θ, 0.0)

        r1 = (i + 3ϕ1/π)*shift_vector + rotation_matrix1*r0
        r2 = (i + 3ϕ2/π)*shift_vector + rotation_matrix2*(2r0)
        r3 = (i + 3ϕ3/π)*shift_vector + rotation_matrix3*(√7*r0)
        r4 = (i + 3ϕ4/π)*shift_vector + rotation_matrix4*(√7*r0)
        r5 = (i + 3ϕ5/π)*shift_vector + rotation_matrix5*(√13*r0)
        r6 = (i + 3ϕ6/π)*shift_vector + rotation_matrix6*(√13*r0)

        A = Carbon(vec(r1...), orbital"2s", orbital"2p_x", orbital"2p_y", orbital"2p_z", axes=axes)
        B = Carbon(vec(r2...), orbital"2s", orbital"2p_x", orbital"2p_y", orbital"2p_z", axes=axes)
        if i == 2
            C = Nitrogen(vec(r3...), orbital"2s", orbital"2p_x", orbital"2p_y", orbital"2p_z", axes=axes)
            F = Carbon(vec(r6...), orbital"2s", orbital"2p_x", orbital"2p_y", orbital"2p_z", axes=axes)
        else
            C = Carbon(vec(r3...), orbital"2s", orbital"2p_x", orbital"2p_y", orbital"2p_z", axes=axes)
            F = Hydrogen(vec(r6...), orbital"1s", axes=axes)
        end

        if i == 3
            D = Nitrogen(vec(r4...), orbital"2s", orbital"2p_x", orbital"2p_y", orbital"2p_z", axes=axes)
            E = Carbon(vec(r5...), orbital"2s", orbital"2p_x", orbital"2p_y", orbital"2p_z", axes=axes)
        else
            D = Carbon(vec(r4...), orbital"2s", orbital"2p_x", orbital"2p_y", orbital"2p_z", axes=axes)
            E = Hydrogen(vec(r5...), orbital"1s", axes=axes)
        end

        push!(mol, A, B, C, D, E, F)
        push!(mol, Bond(A, A′), Bond(A, B), Bond(B, C), Bond(C, D), Bond(B, D′), Bond(C, E), Bond(D, F))

        if i == 3
            ϕ9 = atan(√3, 4)
            ϕ10 = -ϕ9
            rotation_matrix9 = LinearAlgebra.Givens(1, 2, cos(i*θ + ϕ9), -sin(i*θ + ϕ9))
            rotation_matrix10 = LinearAlgebra.Givens(1, 2, cos(i*θ + ϕ10), -sin(i*θ + ϕ10))
            r7 = (i + 3ϕ1/π)*shift_vector + rotation_matrix1*(4r0)
            r8 = (i + 3ϕ1/π)*shift_vector + rotation_matrix1*(5r0)
            r9 = (i + 3ϕ9/π)*shift_vector + rotation_matrix9*(√19*r0)
            r10 = (i + 3ϕ10/π)*shift_vector + rotation_matrix10*(√19*r0)

            G = Carbon(vec(r7...), orbital"2s", orbital"2p_x", orbital"2p_y", orbital"2p_z", axes=axes)
            H = Hydrogen(vec(r8...), orbital"1s", axes=axes)

            I = Hydrogen(vec(r9...), orbital"1s", axes=axes)
            J = Hydrogen(vec(r10...), orbital"1s", axes=axes)
            push!(mol, G, H, I, J)
            push!(mol, Bond(E, G), Bond(F′, G), Bond(G, H), Bond(E, I), Bond(F′, J))
        end

        A′ = A
        D′ = D
        F′ = F
    end

    i = N + 1
    rotation_matrix1 = LinearAlgebra.Givens(1, 2, cos(i*θ + ϕ1), -sin(i*θ + ϕ1))
    rotation_matrix2 = LinearAlgebra.Givens(1, 2, cos(i*θ + ϕ2), -sin(i*θ + ϕ2))
    rotation_matrix3 = LinearAlgebra.Givens(1, 2, cos(i*θ + ϕ3), -sin(i*θ + ϕ3))
    rotation_matrix4 = LinearAlgebra.Givens(1, 2, cos((i+1)*θ + ϕ1), -sin((i+1)*θ + ϕ1))
    rotation_matrix5 = LinearAlgebra.Givens(1, 2, cos((i+2)*θ + ϕ1), -sin((i+2)*θ + ϕ1))
    axes = makeaxes(π/2 + i*θ, 0.0)

    r1 = (i + 3ϕ1/π)*shift_vector + rotation_matrix1*r0
    r2 = (i + 3ϕ2/π)*shift_vector + rotation_matrix2*(r0 + [abs(l), 0, 0])
    r3 = (i + 3ϕ3/π)*shift_vector + rotation_matrix3*([abs(l)*√7, 0, 0])
    r4 = ((i+1) + 3ϕ1/π)*shift_vector + rotation_matrix4*r0
    r5 = ((i+2) + 3ϕ1/π)*shift_vector + rotation_matrix5*r0

    A = Carbon(vec(r1...), orbital"2s", orbital"2p_x", orbital"2p_y", orbital"2p_z", axes=axes)
    B = Carbon(vec(r2...), orbital"2s", orbital"2p_x", orbital"2p_y", orbital"2p_z", axes=axes)
    C = Hydrogen(vec(r3...), orbital"1s", axes=axes)
    D = Oxygen(vec(r4...), orbital"2s", orbital"2p_x", orbital"2p_y", orbital"2p_z", axes=axes)
    E = Hydrogen(vec(r5...), orbital"1s", axes=axes)
    push!(mol, A, B, C, D, E)
    push!(mol, Bond(A, A′), Bond(A, B), Bond(B, C), Bond(B, D′), Bond(A, D), Bond(D, E))
    return mol
end
