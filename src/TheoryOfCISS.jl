module TheoryOfCISS

using InlineExports

# using AtomicSimulationEnvironment
# @export const ase = AtomicSimulationEnvironment

using BSON: @save, @load
using Dates: day, month, now, today, year
using LinearAlgebra
using Molecules
using Molecules: Molecule, Bond,
    Hydrogen, Carbon, Nitrogen, Oxygen,
    Lead, angularmomentum, hamiltonian, makeaxes, mass, nB, overlap
using ProgressMeter
using Random
using SlaterKoster: loaddir
using SparseArrays
using StaticArrays
using Statistics

@export const σ0 = [1.0 0.0; 0.0 1.0]
@export const σ = (
    x = [0.0 1.0; 1.0 0.0],
    y = [0.0 -1.0im; 1.0im 0.0],
    z = [1.0 0.0; 0.0 -1.0]
)
@export const ⊗ = kron

vec(x, y, z) = SVector(x, y, z)

abstract type Simulation end

include("globals.jl")

include("mol.jl")
include("continuous.jl")

include("lite.jl")
include("lite_datagen.jl")
include("plots_recipes.jl")

# include("../extra/scratch.jl")

# include("basic.jl")

end # module
