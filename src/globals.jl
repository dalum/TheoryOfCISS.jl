# Units
const ENERGY_UNIT = "eV"
const LENGTH_UNIT = "Å"
const TIME_UNIT = "ps"

const PARAMDIR = abspath(@__DIR__, "..", "params")
#const GLOBAL_SKT = loaddir(joinpath(PARAMDIR, "ob2-1-1", "base"))
const GLOBAL_SKT = loaddir(joinpath(PARAMDIR, "mio-1-1"))
const SLATERKOSTER_HAMILTONIAN = SlaterKoster.hamiltonian(GLOBAL_SKT)
const SLATERKOSTER_OVERLAP = SlaterKoster.overlap(GLOBAL_SKT)

const CC_SINGLE_BOND_LENGTH = 1.54
const CC_DOUBLE_BOND_LENGTH = 1.34
const CC_TRIPLE_BOND_LENGTH = 1.20
const CH_BOND_LENGTH = 1.09

const DEFAULT_α = 1.0
# In units of `electronvolt`
const DEFAULT_λ = 6e-3

# In units of `electronvolt * picosecond`
const HBAR = 6.582119569e-4
