const PARAMDIR = abspath(@__DIR__, "..", "params")
#const GLOBAL_SKT = loaddir(joinpath(PARAMDIR, "ob2-1-1", "base"))
const GLOBAL_SKT = loaddir(joinpath(PARAMDIR, "mio-1-1"))

const CC_SINGLE_BOND_LENGTH = 1.54
const CC_DOUBLE_BOND_LENGTH = 1.34
const CC_TRIPLE_BOND_LENGTH = 1.20
const CH_BOND_LENGTH = 1.09

const DEFAULT_α = 100.0
const DEFAULT_λ = 6e-3
