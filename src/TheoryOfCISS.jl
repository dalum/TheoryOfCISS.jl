module TheoryOfCISS

using InlineExports

using BSON: @save, @load
using Dates: day, month, now, today, year
using ForwardDiff: Dual, derivative
using LinearAlgebra
using Molecules
using Molecules: Molecule, Bond,
    Hydrogen, Carbon, Nitrogen, Oxygen,
    Lead, angularmomentum, hamiltonian, makeaxes, mass, nB
using ProgressMeter
using Random
using Roots: Roots, find_zero
using ScientificNotation
using SlaterKoster: loaddir
using SparseArrays
using StaticArrays
using Statistics

import Plots

const PARAMDIR = abspath(@__DIR__, "..", "params")
DEFAULT_LENGTH = 1.4

vec(N, n, x, y, z) = SVector(x, y, z)

function dvec(N, n, x, y, z)
    return SVector(
        Dual(x, Tuple(1.0*(3n - 2 == i) for i in 1:3N)),
        Dual(y, Tuple(1.0*(3n - 1 == i) for i in 1:3N)),
        Dual(z, Tuple(1.0*(3n == i) for i in 1:3N)))
end

abstract type Simulation end

include("mol.jl")

include("lite.jl")
include("lite-datagen.jl")

include("basic.jl")
include("continuous.jl")
include("datagen.jl")

end # module
