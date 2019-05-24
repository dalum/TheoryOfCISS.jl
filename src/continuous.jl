function makelattice(ϕ=π/4, l=1.4; handedness=1, θ=nothing, δz=nothing)
    mol = Molecule()
    mol′ = Molecule()
    if δz === nothing
        δz = √3 * l * sin(ϕ) / 2
    end
    shift_vector = [0, 0, δz]
    if θ === nothing
        θ = π - 2atan(√3 * cos(ϕ))*handedness
    end

    r1 = [√(l^2 - δz^2) / 2sin(θ/2), 0, 0]
    r2 = r1 + [abs(l), 0, 0]
    axes = makeaxes(π/2, ϕ)
    A = Carbon(r1, orbital"2s", orbital"2p_x", orbital"2p_y", orbital"2p_z", axes=axes)
    B = Hydrogen(r2, orbital"1s", axes=axes)
    push!(mol, A, B)
    push!(mol, Bond(A, B))

    rotation_matrix = LinearAlgebra.Givens(1, 2, cos(θ), -sin(θ))
    r1′ = shift_vector + rotation_matrix*r1
    r2′ = shift_vector + rotation_matrix*r2
    axes = makeaxes(π/2 + θ, ϕ)

    A′ = Carbon(r1′, orbital"2s", orbital"2p_x", orbital"2p_y", orbital"2p_z", axes=axes)
    B′ = Hydrogen(r2′, orbital"1s", axes=axes)
    push!(mol′, A′, B′)
    push!(mol′, Bond(A′, B′))

    return Molecules.UnitCell(mol, mol′, Set([Bond(A, A′)]))
end
