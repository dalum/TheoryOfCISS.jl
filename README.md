# TheoryOfCISS

Code repository accompanying the paper: *Theory of Chiral Induced Spin
Selectivity* by Sakse Dalum and Per Hedegård.

# Installation

From the Julia REPL, run:
```julia
]dev https://github.com/dalum/SlaterKoster.jl
]dev https://github.com/dalum/Molecules.jl
]dev https://github.com/dalum/TheoryOfCISS.jl
```
This will install the package and its dependencies in your `$HOME/.julia/dev/` directory.
The `Manifest.toml` file contains the relative paths `../Molecules/` and `../SlaterKoster/` for the location of the dependencies.
If, for some reason, the dependencies are in another location, please modify the `Manifest.toml` to point to the correct paths.

# Usage

Generate data for a molecule of length N = 24, for `x` (geometry parameter `ϕ`) on the interval `[0.0, π/2]`, and `y` (energy `E`) on the interval `[0.0, 3.0]`, sampling 100 points along each axis:
```julia
# Generate data
cols = gen_bands(24, xsymbol = :ϕ, x_len = 100, x_lower = 0.7, x_upper = 1.2)
gen_near_bands!(cols, y_len = 100, y_lower = 0.0, y_upper = 3.0, f = TheoryOfCISS.calc_data1)

# Plot
heatmap(percentage∘polarization, cols, xtransform=ϕtoθ)
plot!(cols, xtransform=ϕtoθ, xlabel="φ [rad.]", ylabel="E [eV]")
```
