let ntraces = 10, max_N = 20, dat = fill(0.0, max_N, ntraces)
    for (idx, seed) in enumerate((abs(rand(Int)) for _ in 1:ntraces))
        cols = gen_bands(Cumulene(θ = π/2 + π/100), xsymbol = :N, bounds=(1, max_N), nsamples = 20, seed=seed)
        gen_near_bands!(cols, ysymbol = :E, bounds=(-10, 0), nsamples = 200, f = TheoryOfCISS.calc_data3, threshold = 1.0)
        xs, ys, zs = TheoryOfCISS.datamap(x -> x[:p], TheoryOfCISS.smoothen(cols, σ=10, w=TheoryOfCISS.box, p = x -> fullpolarization_reduced(x, 3)))
        dat[:, idx] = zs[100, :]
    end
    dat
end

#plot(100 .* zs[4, :], xlabel="N", ylabel="polarization (%)")

let ntraces = 10, Δ = π/8, max_N = 5, dat = fill(0.0, max_N, ntraces)
    for (idx, θ) in enumerate(range(π/2 - Δ, π/2 + Δ, length=ntraces))
        cols = gen_bands(Cumulene(θ = θ), xsymbol = :N, bounds=(1, max_N), nsamples = max_N, rng=1.0I)
        gen_near_bands!(cols, ysymbol = :E, bounds=(-10, 0), nsamples = 200, f = TheoryOfCISS.calc_data3, threshold = 1.0)
        xs, ys, zs = TheoryOfCISS.datamap(x -> x[:p], TheoryOfCISS.smoothen(cols, σ=10, w=TheoryOfCISS.box, p = x -> fullpolarization_reduced(x, 3)))
        dat[:, idx] = zs[100, :]
    end
    dat
end
