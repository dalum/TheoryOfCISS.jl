# TheoryOfCISS

Code repository accompanying the paper: "Theory of Chiral Induced Spin
Selectivity".

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
