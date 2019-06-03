module TheoryOfCISS

using InlineExports

using BSON: @save, @load
using Dates: day, month, now, today, year
using LinearAlgebra
using Molecules
using Molecules: Molecule, Bond,
    Hydrogen, Carbon, Nitrogen, Oxygen,
    Lead, angularmomentum, hamiltonian, makeaxes, mass, nB
using ProgressMeter
using Random
using SlaterKoster: loaddir
using SparseArrays
using StaticArrays
using Statistics

import Plots

@export const ⊗ = kron

const PARAMDIR = abspath(@__DIR__, "..", "params")
DEFAULT_LENGTH = 1.4

vec(x, y, z) = SVector(x, y, z)

abstract type Simulation end

include("mol.jl")

include("lite.jl")
include("lite_datagen.jl")
include("plots_recipes.jl")

# include("basic.jl")
# include("continuous.jl")

end # module
