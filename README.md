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

Load the packages:
```julia
using TheoryOfCISS
```
Generate data for a helicene molecule with `N = 7` rings ([7] helicene), for `x` (geometry parameter `δz` in units of Ångström) on the interval `[0.5, 1.0]`, and `y` (energy `E` in units of electronvolts) on the interval `[-5.0, 0.0]`, sampling 100 points along each axis:
```julia
cols = gen_bands(Helicene(N=7), xsymbol=:δz, bounds=(0.5, 1.0), nsamples=100);
gen_near_bands!(cols, ysymbol=:E, bounds=(-5, 0), nsamples=100);
```
Plot the polarization along the z-axis as a heatmap, and overlay the energies of the molecular states on top:
```julia
using Plots

heatmap(percentage∘polarization(3), cols);
plot!(cols, xguide="δz [Å]", yguide="E [eV]")
```
