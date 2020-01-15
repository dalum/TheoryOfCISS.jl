struct CoupledMolecule <: Simulation
    constructor::AbstractMoleculeConstructor
    mol::Molecule
    eig
    μ
    Δ
    H
    H0
    S0
    γ
    Γ
    L
end
Base.show(io::IO, sim::LiteSimulation) = Base.print(io, "LiteSimulation(...)")
