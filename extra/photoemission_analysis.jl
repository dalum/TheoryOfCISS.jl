using Measurements

d(a_in, a_out) = (1 - a_in*a_in) / (1 - a_in*a_out) * a_out - a_in
d̃(a_in, a_out) = a_out - a_in

# Au(111)
# a_in, a_out...

P_in = [
    # Au(111)
     22    0  -22
    # Au(poly.)
     -4    2    6
     -4    2    6
     -4    2    6
    # Au(poly.) (not reported)
     -4    2    6
    # Al(poly.)
      0    0    0
    # Au(poly.) (not reported)
     -4    2    6
     -4    2    6
     -4    2    6
    # Cu(332)
      0    0    0
      0    0    0
    # Ag(110)
      3    0   -3
      3    0   -3
    # Au(110)
     27    0  -24
     27    0  -24
    # Au(110)
    missing  0  missing
    missing  0  missing
] ./ 100

σ_in = [
    5  5  5
    5  5  5
    5  5  5
    5  5  5
    5  5  5
    5  5  5
    5  5  5
    5  5  5
    5  5  5
    5  5  5
    5  5  5
    5  5  5
    5  5  5
    5  5  5
    5  5  5
    missing  5  missing
    missing  5  missing
] ./ 100


P_out = [
    # Au(111)
    -29 -31 -35
    # Au(poly.)
    -40 -38 -35
    -39 -36 -32
    -61 -57 -55
    # Au(poly.)
      8  14  17
    # Al(poly.)
     14.9 15 15.1
    # Au(poly.)
    -17 -11 -6
    -20 -14 -10
    -22 -18 -14
    # Cu(332)
    -6.3 -6.7 -6
     11  12  12
    # Ag(110)
     -7  -9 -12
    6.2 7.1 6.3
    # Au(111)
     22  -8 -35
     34   8 -16
    # Au(111)
    missing  -12  missing
    missing    4  missing
] ./ 100

σ_out = [
    5  5  5
    5  5  5
    5  5  5
    5  5  5
    5  5  5
    5  5  5
    5  5  5
    5  5  5
    5  5  5
    5  5  5
    5  5  5
    5  5  5
    5  5  5
    5  5  5
    5  5  5
    missing  5  missing
    missing  5  missing
] ./ 100

# Polycrystalline
# a_in, a_out...

D(a_in, a_out, idx=2) = (a_out[idx] - a_in[idx]) / (1 - a_in[idx]*a_out[idx])
A(idx, a_in, a_out) = (1 - a_in[idx]*a_out[idx]) * D(a_in, a_out) - (a_out[idx] - a_in[idx])
ΔT(a_in, a_out) = (A(3, a_in, a_out) + A(1, a_in, a_out)) / (a_out[3] - a_out[1])
ΔD(a_in, a_out) = ΔT(a_in, a_out) * (a_out[3] + a_out[1])/2 - (A(3, a_in, a_out) - A(1, a_in, a_out)) / 2

D2(idx, a_in, a_out) = (a_out[idx] - a_in[idx] + (idx-2)*(ΔT(a_in, a_out)*a_out[idx] - ΔD(a_in, a_out))) / (1 - a_in[idx]*a_out[idx])

Dcons(a_in, a_out) = (D(a_in, a_out, 1) + D(a_in, a_out, 3)) / 2 - D(a_in, a_out, 2)
