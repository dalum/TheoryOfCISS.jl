using Measurements

# Figure 2A-D should be ordered:
#   ccw, lin, cw
# to be consistent with Figure 4.

# N-matrix, manually counted from bins:
# rows: cw, lin, ccw
# cols: in, out-50-bp, out-78-bp
N = [  63  54  48
      125 109  96
       61  55  48 ]

a = [ 0.22  -0.29  -0.545
      0.0   -0.31  -0.572
     -0.22  -0.35  -0.608 ]

σ = [ 0.05  0.03  0.07
      0.05  0.04  0.058
      0.05  0.03  0.059 ]

σ_mean = σ ./ sqrt.(N)

ā = a .± σ_mean

d(a_in, a_out) = (1 - a_in*a_in) / (1 - a_in*a_out) * a_out - a_in
d̃(a_in, a_out) = a_out - a_in

D = d.(ā[:, 1], ā[:, 2:3])
D̃ = d̃.(ā[:, 1], ā[:, 2:3])

Δ = D[2, :] .- (D[1, :] .+ D[3, :]) ./ 2
Δ̃ = D̃[2, :] .- (D̃[1, :] .+ D̃[3, :]) ./ 2

