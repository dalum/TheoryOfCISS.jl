@export function scr_test(sim; λ=DEFAULT_λ, indices=[37], weights=[1])
    ρ0 = sim.γL
    ρ0 = sum(zip(indices, weights)) do (i, w)
        v0 = eigen(sim.H).vectors[:,i]
        w*v0*v0'
    end
    Γ = sim.ΓR

    J = sim.H - im*Γ/2
    j, U = eigen(J)
    N = size(J, 1)

    _ρ_ = inv(U)*ρ0*inv(U')
    _Γ_ = U'Γ*U
    _L_ = [λ * U'L*inv(U') for L in sim.L]

    # function β(t)
    #     As = Matrix{Complex{Float64}}[]
    #     for L in _L_
    #         A = Matrix{Complex{Float64}}(undef, N, N)

    #         for idx in CartesianIndices(A)
    #             l, m = idx[1], idx[2]
    #             if j[l] == j[m]
    #                 A[l,m] = -im * t * exp(-im*j[l]*t) * L[l,m]
    #             else
    #                 A[l,m] = inv(j[l]- j[m]) * (exp(-im*j[l]*t)*L[l,m] - L[l,m]*exp(-im*j[m]*t))
    #             end
    #         end

    #         push!(As, U*A*inv(U))
    #     end
    #     return As
    # end

    function f(t)
        integrated_T = sum(CartesianIndices((N, N))) do idx
            k, l = idx[1], idx[2]
            x = j[k] - j[l]'
            if iszero(x)
                return _ρ_[k,l] * _Γ_[l,k] * t
            else
                return _ρ_[k,l] * _Γ_[l,k] * im*(exp(-im*x*t) - 1) / x
            end
        end

        integrated_S = sum(CartesianIndices((N, N, N))) do idx
            l, m, n = idx[1], idx[2], idx[3]
            x = j[l]' - j[m]'
            if iszero(x)
                a = j[n] - j[m]'
                return _ρ_[n,l] * _L_[3][l,m] * _Γ_[m,n] * im*(
                    im*t*exp(-im*a*t) / a +
                    (exp(-im*a*t) - 1) / a^2
                )
            else
                a = j[n] - j[l]'
                b = j[n] - j[m]'
                return _ρ_[n,l] * _L_[3][l,m] * _Γ_[m,n] * inv(b - a) * im*(
                    (exp(-im*a*t) - 1) / a - (exp(-im*b*t) - 1) / b
                )
                # return _ρ_[n,l] * _L_[3][l,m] * _Γ_[m,n] * inv(j[l]' - j[m]') * im*(
                #     (exp(-im*(j[n] - j[l]')*t) - 1) / (j[n] - j[l]') -
                #     (exp(-im*(j[n] - j[m]')*t) - 1) / (j[n] - j[m]')
                # )
            end
        end

        return integrated_T, integrated_S
        u = LinearAlgebra.exp!(-im*J*t)
        A = ρ0*u'Γ
        return tr(A*u), integrated_T, [tr(A*b) for b in β(t)], integrated_S
    end

    return f
end


@export function scr_test_full(sim; λ=DEFAULT_λ, indices=[Int(length(sim.eig.values) / 2 + 1)], weights=[1], α=1.0)
    ρ0 = sum(zip(indices, weights)) do (i, w)
        v0 = eigen(sim.H).vectors[:,i]
        w*v0*v0'
    end
    Γ = sim.ΓL
    Σ = - im * (α*sim.ΓR + sim.ΓL) / 2

    J = (sim.H + Σ) ⊗ σ0 + sum(λ * map(⊗, sim.L, σ))
    j, U = eigen(J)

    _ρ_ = inv(U) * (ρ0 ⊗ σ0) * inv(U')
    _Γ_ = U' * (Γ ⊗ σ0) * U
    _Γi_ = [U' * (Γ ⊗ σi) * U for σi in σ]

    N = size(J, 1)

    function f(t)
        integrated_T = sum(CartesianIndices((N, N))) do idx
            k, l = idx[1], idx[2]
            x_kl = im*(j[k] - j[l]')
            if iszero(x_kl)
                return _ρ_[k,l] * _Γ_[l,k] * t
            else
                return _ρ_[k,l] * _Γ_[l,k] * (1 - exp(-x_kl * t)) / x_kl
            end
        end

        integrated_S = sum(CartesianIndices((N, N))) do idx
            k, l = idx[1], idx[2]
            x_kl = im*(j[k] - j[l]')
            if iszero(x_kl)
                return _ρ_[k,l] * _Γi_[3][l,k] * t
            else
                return _ρ_[k,l] * _Γi_[3][l,k] * (1 - exp(-x_kl * t)) / x_kl
            end
        end

        # u = LinearAlgebra.exp!(-im*J*t)
        # A = (ρ0 ⊗ σ0) * u' * (Γ ⊗ σ0)
        return integrated_T, integrated_S#tr(A*u)
    end

    return f
end

# Magnet
@export function scr_test_magnet(sim; λ=DEFAULT_λ, indices=[Int(length(sim.eig.values) / 2 + 1)], weights=[1], α=1.0, m=0.5)
    ρ0 = sum(zip(indices, weights)) do (i, w)
        v0 = eigen(sim.H).vectors[:,i]
        w*v0*v0'
    end
    Γ = α*sim.ΓR
    Σ = - im * (α*sim.ΓR + sim.ΓL) / 2
    ΔΣ = - im * (m*sim.γL) / 2

    J = ((sim.H + Σ) ⊗ σ0 + ΔΣ ⊗ σ[:z] + sum(λ * map(⊗, sim.L, σ)),
         (sim.H + Σ) ⊗ σ0 - ΔΣ ⊗ σ[:z] + sum(λ * map(⊗, sim.L, σ)))
    j, U = begin
        tmp = eigen.(J)
        map(x -> x.values, tmp), map(x -> x.vectors, tmp)
    end

    ρ = map(U -> inv(U) * (ρ0 ⊗ σ0) * inv(adjoint(U)), U)
    Γ = map(U -> adjoint(U) * (Γ ⊗ σ0) * U, U)

    N = size(J[1], 1)

    function f(t)
        integrated_T = map(j, ρ, Γ) do j, ρ, Γ
            T = sum(CartesianIndices((N, N))) do idx
                k, l = idx[1], idx[2]
                x_kl = im*(j[k] - j[l]')
                real(x_kl) < 0 && error("negative `real(x_kl)` found: $(real(x_kl))")
                return ρ[k,l] * Γ[l,k] * (1 - exp(-x_kl * t)) / x_kl
            end
            real(T) < 0 && error("negative `real(T)` found: $(real(T))")
            real(T)
        end

        return integrated_T
    end

    return f
end

# Calculate the photon emission
@export function scr_photon_emission_full(
    sim;
    λ=DEFAULT_λ,
    indices=[37],
    weights=[1],
    seed = 824,
    W = 1e-1,
    rng = MersenneTwister(seed)
)
    ρ0 = sum(zip(indices, weights)) do (i, w)
        v0 = sim.eig.vectors[:,i]
        w*v0*v0'
    end
    Γ = sim.ΓR

    J = (sim.H - im*Γ/2) ⊗ σ0 + sum(λ * Li ⊗ σi for (Li, σi) in zip(sim.L, σ))
    j, U = eigen(J)
    N = size(J, 1)

    V = let
        local i, j = rand(rng, 1:size(sim.H, 1), 2)
        while i == j i, j = rand(rng, 1:N, 2) end
        W*sim.eig.vectors[:,i]*sim.eig.vectors[:,j]'
    end

    _ρ_ = inv(U) * (ρ0 ⊗ σ0) * inv(U')
    _Γ_ = U' * (Γ ⊗ σ0) * U
    _V_ = inv(U) * (V ⊗ σ0) * U

    function f()
        K(a, b) = inv(j[a]' - j[b])
        K_(a, b) = inv(j[a]' - j[b]')

        output = sum(CartesianIndices((N, N, N, N))) do idx
            k, l, m, n = idx[1], idx[2], idx[3], idx[4]

            X = 2π*im * _ρ_[k,l] * _V_'[l,m] * _Γ_[m,n] * _V_[m,n]

            if iszero(j[l]' - j[m]')
                return -X * K(l,n) * K(m,k) * (K(l,n) + K(m,k))
            else
                return X * K_(l,m) * (K(l,n)*K(l,k) - K(m,n)*K(m,k))
            end
        end

        return output
    end

    return f
end

##################################################

@export function scr_magnet(sim; E=sim.μ, λ=DEFAULT_λ, m=0.5, n=1, α=0.0, ϕ=0)
    ρ0 = let
        v0 = eigenstate(sim, E = E, ϕ=ϕ)
        0.5 * v0*v0'
    end

    Γ = sim.ΓL
    Σ = - im * (α*sim.ΓL + sim.ΓR) / 2
    ΔΣ = - im * (m*sim.ΔΓR) / 2

    J = ((sim.H + Σ) ⊗ σ0 + ΔΣ ⊗ σ[:z] + sum(λ * map(⊗, sim.L, σ)),
         (sim.H + Σ) ⊗ σ0 - ΔΣ ⊗ σ[:z] + sum(λ * map(⊗, sim.L, σ)))
    j, U = begin
        tmp = eigen.(J)
        map(x -> x.values, tmp), map(x -> x.vectors, tmp)
    end

    ρ = map(U -> inv(U) * (ρ0 ⊗ σ0) * inv(adjoint(U)), U)
    Γ = map(U -> adjoint(U) * (Γ ⊗ σ0) * U, U)

    N = size(J[1], 1)

    Ts = map(j, ρ, Γ) do j, ρ, Γ
        T = sum(CartesianIndices((N, N))) do idx
            k, l = idx[1], idx[2]
            x_kl = im*(j[k] - j[l]')
            real(x_kl) < 0 && error("negative `real(x_kl)` found: $(real(x_kl))")
            return ρ[k,l] * Γ[l,k] / x_kl
        end
        real(T) < 0 && error("negative `real(T)` found: $(real(T))")
        real(T)
    end

    return (Ts = Ts,)
end

@export function time_dependent_escape(
    sim;
    E = sim.μ,
    λ = DEFAULT_λ,
    m = 1.0,
    α = DEFAULT_α,
    ϕ = 0,
    hbar = HBAR
)
    ρ0 = let
        v0 = eigenstate(sim, E = E, ϕ=ϕ)
        v0*v0'
    end

    Γ = α*sim.ΓR
    Σ = - im * α*sim.ΓR / 2

    J = (sim.H + Σ) ⊗ σ0 + sum(λ * map(⊗, sim.L, σ))
    j, U = eigen(J)
    N = length(j)

    ρ = inv(U) * (ρ0 ⊗ ((σ0 + m * σ[:z]) / 2)) * inv(adjoint(U))
    Γ = adjoint(U) * (Γ ⊗ σ0) * U

    ξs = map(CartesianIndices((N, N))) do idx
        k, l = idx[1], idx[2]
        x_kl = im*(j[k] - j[l]')
        real(x_kl) < 0 && error("negative `real(x_kl)` found: $(real(x_kl))")
        return x_kl
    end
    Ts = map(CartesianIndices((N, N))) do idx
        k, l = idx[1], idx[2]
        x_kl = im*(j[k] - j[l]')
        real(x_kl) < 0 && error("negative `real(x_kl)` found: $(real(x_kl))")
        return ρ[k,l] * Γ[l,k] / x_kl
    end
    fs = map(CartesianIndices((N, N))) do idx
        k, l = idx[1], idx[2]
        x_kl = im*(j[k] - j[l]')
        real(x_kl) < 0 && error("negative `real(x_kl)` found: $(real(x_kl))")
        return f(t) = 1 - exp(-x_kl * t / hbar)
    end

    return (ξs = ξs, Ts = Ts, fs = fs)
end

@export function scr_f(
    δx = 0.0,
    δy = 0.1;
    t0 = 0,
    t1 = 1,
    nsamples = 500,
    ts = range(t0, t1, length=nsamples),
    seed = 42,
    α = 1,
)
    x0, y0, y1 = let
        xsymbol = :δz
        x0 = 0.8
        y0 = -4.2
        mol = Helicene(N=7)

        x0, y0, y1 = TheoryOfCISS.find_crossing(
            x0, y0, mol,
            xsymbol = xsymbol,
            threshold=1e-12,
            attenuation=1e-2,
            rounding = :both
        )
    end

    x = x0
    y = (y0 + y1) / 2

    cols_up = gen_bands(Helicene(N=7), xsymbol = :δz, bounds=(x + δx, x + δx), nsamples=1, seed=seed);
    cols_down = gen_bands(Helicene(N=7), xsymbol = :δz, bounds=(x + δx, x + δx), nsamples=1, seed=seed);

    gen_on_bands!(cols_up, ysymbol=:E, f = time_dependent_escape, m=1.0, bounds=(y - δy, y + δy), α=α);
    gen_on_bands!(cols_down, ysymbol=:E, f = time_dependent_escape, m=-1.0, bounds=(y - δy, y + δy), α=α);

    _f(t, x) = sum(map((T, f) -> T * f(t), x.Ts, x.fs))
    f_x(t, i) = real(_f(t, cols_up[1].data[i])), real(_f(t, cols_down[1].data[i]))
    f_xs = map(i -> mapreduce(t -> [f_x(t, i)...]', vcat, ts), 1:2)
    pol(a, b) = 100 * !iszero(a - b) * ((a - b) / (a + b))
    dat = map(x -> map(i -> pol(i...), eachrow(x)), f_xs)
    return ts, f_xs, dat
end
