@export ϕtoθ(ϕ) = π - 2atan(√3 * cos(ϕ))
@export ϕtoδz(ϕ, l=1.4) = √3 * l * sin(ϕ) / 2

"""
    makemol(N, ϕ, l)
"""
@export function makemol(N, ϕ=π/4, l=1.4; handedness=1, δz=nothing, θ=nothing, r0=nothing)
    if δz !== nothing && θ !== nothing
        r0 = [l, 0, 0]
    end

    if δz === nothing
        δz = √3 * l * sin(ϕ) / 2
    end
    if θ === nothing
        θ = π - 2atan(√3 * cos(ϕ))*handedness
    end

    if r0 === nothing
        r0 = [√(l^2 - δz^2) / 2sin(θ/2), 0, 0]
    end
    r0::Vector

    mol = Molecule()

    axes = makeaxes(π/2, ϕ*handedness)
    A′ = Hydrogen(vec(r0...), orbital"1s", axes=axes)
    push!(mol, A′)
    for i in 1:N
        rotation_matrix = LinearAlgebra.Givens(1, 2, cos(i*θ), -sin(i*θ))
        shift_vector = [0, 0, i*δz]
        axes = makeaxes(π/2 + i*θ, ϕ*handedness)
        r1 = shift_vector + rotation_matrix*r0
        r2 = shift_vector + rotation_matrix*(r0 + [abs(l), 0, 0])
        A = Carbon(vec(r1...), orbital"2s", orbital"2p_x", orbital"2p_y", orbital"2p_z", axes=axes)
        B = Hydrogen(vec(r2...), orbital"1s", axes=axes)
        push!(mol, A, B)
        push!(mol, Bond(A, A′), Bond(A, B))
        A′ = A
    end
    rotation_matrix = LinearAlgebra.Givens(1, 2, cos((N+1)*θ), -sin((N+1)*θ))
    shift_vector = [0, 0, (N+1)*δz]
    axes = makeaxes(π/2 + (N+1)*θ, ϕ*handedness)
    r1 = shift_vector + rotation_matrix*r0
    A = Hydrogen(vec(r1...), orbital"1s", axes=axes)
    push!(mol, A)
    push!(mol, Bond(A, A′))
    return mol
end

function makehelicene(N, ϕ=π/4, l=1.4; handedness=1)
    δz = √3 * l * sin(ϕ) / 2
    shift_vector = [0, 0, δz]
    θ = 2π/6#π - 2atan(√3 * cos(ϕ))*handedness
    ϕ1 = 0
    ϕ2 = 0
    ϕ3 = π/9
    ϕ4 = 2π/9
    ϕ5 = atan(7/2, √3/2) - 2π/6
    ϕ6 = atan(5/2, 3*√3/2)
    r0 = [l, 0, 0]#[√(l^2 - δz^2) / 2sin(θ/2), 0, 0]

    mol = Molecule()

    axes = makeaxes(π/2, ϕ)
    rotation_matrix1 = LinearAlgebra.Givens(1, 2, cos(ϕ1), -sin(ϕ1))
    rotation_matrix4 = LinearAlgebra.Givens(1, 2, cos(ϕ4), -sin(ϕ4))

    i = 0
    r1 = (i + 3ϕ1/π)*shift_vector + rotation_matrix1*r0
    r4 = (i + 3ϕ4/π)*shift_vector + rotation_matrix4*([abs(l)*√(5 + √3), 0, 0])
    A′ = Hydrogen(vec(r1...), orbital"1s", axes=axes)
    D′ = Hydrogen(vec(r4...), orbital"1s", axes=axes)
    push!(mol, A′)
    push!(mol, D′)
    for i in 1:N
        rotation_matrix1 = LinearAlgebra.Givens(1, 2, cos(i*θ + ϕ1), -sin(i*θ + ϕ1))
        rotation_matrix2 = LinearAlgebra.Givens(1, 2, cos(i*θ + ϕ2), -sin(i*θ + ϕ2))
        rotation_matrix3 = LinearAlgebra.Givens(1, 2, cos(i*θ + ϕ3), -sin(i*θ + ϕ3))
        rotation_matrix4 = LinearAlgebra.Givens(1, 2, cos(i*θ + ϕ4), -sin(i*θ + ϕ4))
        rotation_matrix5 = LinearAlgebra.Givens(1, 2, cos(i*θ + ϕ5), -sin(i*θ + ϕ5))
        rotation_matrix6 = LinearAlgebra.Givens(1, 2, cos(i*θ + ϕ6), -sin(i*θ + ϕ6))
        axes = makeaxes(π/2 + i*θ, ϕ)

        r1 = (i + 3ϕ1/π)*shift_vector + rotation_matrix1*r0
        r2 = (i + 3ϕ2/π)*shift_vector + rotation_matrix2*(r0 + [abs(l), 0, 0])
        r3 = (i + 3ϕ3/π)*shift_vector + rotation_matrix3*([abs(l)*√7, 0, 0])
        r4 = (i + 3ϕ4/π)*shift_vector + rotation_matrix4*([abs(l)*√7, 0, 0])
        r5 = (i + 3ϕ5/π)*shift_vector + rotation_matrix5*([abs(l)*√13, 0, 0])
        r6 = (i + 3ϕ6/π)*shift_vector + rotation_matrix6*([abs(l)*√13, 0, 0])

        A = Carbon(vec(r1...), orbital"2s", orbital"2p_x", orbital"2p_y", orbital"2p_z", axes=axes)
        B = Carbon(vec(r2...), orbital"2s", orbital"2p_x", orbital"2p_y", orbital"2p_z", axes=axes)
        C = Carbon(vec(r3...), orbital"2s", orbital"2p_x", orbital"2p_y", orbital"2p_z", axes=axes)
        D = Carbon(vec(r4...), orbital"2s", orbital"2p_x", orbital"2p_y", orbital"2p_z", axes=axes)
        E = Hydrogen(vec(r5...), orbital"1s", axes=axes)
        F = Hydrogen(vec(r6...), orbital"1s", axes=axes)
        push!(mol, A, B, C, D, E, F)
        push!(mol, Bond(A, A′), Bond(A, B), Bond(B, C), Bond(C, D), Bond(B, D′), Bond(C, E), Bond(D, F))
        A′ = A
        D′ = D
    end
    i = N + 1
    rotation_matrix1 = LinearAlgebra.Givens(1, 2, cos(i*θ + ϕ1), -sin(i*θ + ϕ1))
    rotation_matrix2 = LinearAlgebra.Givens(1, 2, cos(i*θ + ϕ2), -sin(i*θ + ϕ2))
    rotation_matrix3 = LinearAlgebra.Givens(1, 2, cos(i*θ + ϕ3), -sin(i*θ + ϕ3))
    rotation_matrix4 = LinearAlgebra.Givens(1, 2, cos((i+1)*θ + ϕ1), -sin((i+1)*θ + ϕ1))
    axes = makeaxes(π/2 + i*θ, ϕ)

    r1 = (i + 3ϕ1/π)*shift_vector + rotation_matrix1*r0
    r2 = (i + 3ϕ2/π)*shift_vector + rotation_matrix2*(r0 + [abs(l), 0, 0])
    r3 = (i + 3ϕ3/π)*shift_vector + rotation_matrix3*([abs(l)*√7, 0, 0])
    r4 = ((i+1) + 3ϕ1/π)*shift_vector + rotation_matrix4*r0

    A = Carbon(vec(r1...), orbital"2s", orbital"2p_x", orbital"2p_y", orbital"2p_z", axes=axes)
    B = Carbon(vec(r2...), orbital"2s", orbital"2p_x", orbital"2p_y", orbital"2p_z", axes=axes)
    C = Hydrogen(vec(r3...), orbital"1s", axes=axes)
    D = Hydrogen(vec(r4...), orbital"1s", axes=axes)
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
