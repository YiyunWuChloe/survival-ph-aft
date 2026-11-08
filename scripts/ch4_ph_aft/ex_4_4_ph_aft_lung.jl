# Section 4.4 (example after Right-Censored Data and Maximum Likelihood Estimation) — phase-type accelerated failure time model for the lung data
#
# Run from the code_package folder:
#     julia --project=. scripts/ch4_ph_aft/ex_4_4_ph_aft_lung.jl
#
# Model: T_i = exp(η_i) τ_i,  η_i = β_1 (age_i - 62) + β_2 female_i  (no intercept),  τ_i with baseline survival S_0 and density f_0.
#        S_i(t | x_i) = S_0(t e^{-η_i}),   f_i(t | x_i) = e^{-η_i} f_0(t e^{-η_i}).
# Right-censored log-likelihood (Chapter 4):  ℓ = Σ_i [δ_i log f_i(y_i | x_i) + (1 - δ_i) log S_i(y_i | x_i)].
# Baselines: exponential (= Erlang with m = 1), Erlang(2), Erlang(3), Coxian(2) (all phase-type), and Weibull (for comparison).
# β > 0 lengthens the survival time (T = e^η τ); this is the opposite sign of the simulate_data(effect = :aft) of the Julia package.
# All parameters (β and the log of the rates) are maximised together by Nelder-Mead from several random starting values.
# Standard errors: inverse of the numerical Hessian of the negative log-likelihood at the maximum.
# Checked against R: survreg (exponential and Weibull AFT) and a multi-start optim/optimHess of the same log-likelihood
# for the phase-type baselines (reference_checks_R/ex_4_4_ph_aft_lung_check.R).
#
# Output: figures/ex_4_4_ph_aft_lung.pdf (vector);  output/ex_4_4_ph_aft_lung.csv

using Plots, Printf, Random, LinearAlgebra, RDatasets, DataFrames, CSV

const BLUE   = RGB(0.18, 0.37, 0.52)
const ORANGE = RGB(0.88, 0.47, 0.22)
const GREEN  = RGB(0.18, 0.42, 0.18)
const FIGDIR = normpath(joinpath(@__DIR__, "..", "..", "figures"))
const OUTDIR = normpath(joinpath(@__DIR__, "..", "..", "output"))
mkpath(FIGDIR); mkpath(OUTDIR)

raw = dataset("survival", "lung")
y = Vector{Float64}(raw.Time)
δ = raw.Status .== 2
age62 = Vector{Float64}(raw.Age) .- 62
female = Float64.(raw.Sex .== 2)
n = length(y)
@assert n == 228 && count(δ) == 165
X = hcat(age62, female)

# ---- baseline survival and density -------------------------------------------------------------------
fact(k) = prod(1.0:k; init = 1.0)
function coxian2(u, l1, mu1, mu2)
    a = l1 + mu1; b = mu2
    if abs(a - b) < 1e-9 * max(a, b)
        S = exp(-a * u) * (1 + l1 * u); f = exp(-a * u) * (mu1 + l1 * mu2 * u)
    else
        S = exp(-a * u) + l1 * (exp(-b * u) - exp(-a * u)) / (a - b)
        f = mu1 * exp(-a * u) + l1 * mu2 * (exp(-b * u) - exp(-a * u)) / (a - b)
    end
    return S, f
end
erlang(m) = (u, r) -> (exp(-r[1] * u) * sum((r[1] * u)^j / fact(j) for j in 0:m-1), r[1]^m * u^(m - 1) * exp(-r[1] * u) / fact(m - 1))
baselines = Dict(
    :exp     => (name = "Exponential", nrates = 1, fun = erlang(1)),
    :erlang2 => (name = "Erlang(2)",    nrates = 1, fun = erlang(2)),
    :erlang3 => (name = "Erlang(3)",    nrates = 1, fun = erlang(3)),
    :coxian2 => (name = "Coxian(2)",    nrates = 3, fun = (u, r) -> coxian2(u, r[1], r[2], r[3])),
    :weibull => (name = "Weibull",      nrates = 2, fun = (u, r) -> (exp(-(u / r[2])^r[1]), (r[1] / r[2]) * (u / r[2])^(r[1] - 1) * exp(-(u / r[2])^r[1]))))
# Weibull: r = (shape, scale)

function negloglik(θ, key)
    b = baselines[key]
    β = θ[1:2]; r = exp.(θ[3:end])
    η = X * β
    ℓ = 0.0
    for i in 1:n
        u = y[i] * exp(-η[i])
        S, f = b.fun(u, r)
        ℓ += δ[i] ? log(max(exp(-η[i]) * f, 1e-300)) : log(max(S, 1e-300))
    end
    return -ℓ
end

function nelder_mead(f, x0; step = 0.3, maxit = 40_000, tol = 1e-14)
    d = length(x0)
    Xs = [copy(x0) for _ in 1:d+1]
    for i in 1:d; Xs[i+1][i] += step; end
    F = [f(x) for x in Xs]
    for it in 1:maxit
        idx = sortperm(F); Xs = Xs[idx]; F = F[idx]
        abs(F[end] - F[1]) < tol * (1 + abs(F[1])) && break
        c = sum(Xs[1:d]) ./ d
        xr = c .+ (c .- Xs[end]); fr = f(xr)
        if fr < F[1]
            xe = c .+ 2 .* (c .- Xs[end]); fe = f(xe)
            fe < fr ? (Xs[end] = xe; F[end] = fe) : (Xs[end] = xr; F[end] = fr)
        elseif fr < F[d]
            Xs[end] = xr; F[end] = fr
        else
            xc = fr < F[end] ? c .+ 0.5 .* (xr .- c) : c .+ 0.5 .* (Xs[end] .- c); fc = f(xc)
            if fc < min(fr, F[end])
                Xs[end] = xc; F[end] = fc
            else
                for i in 2:d+1; Xs[i] = Xs[1] .+ 0.5 .* (Xs[i] .- Xs[1]); F[i] = f(Xs[i]); end
            end
        end
    end
    i = argmin(F)
    return Xs[i], F[i]
end

function numerical_hessian(f, x; h = 1e-3)
    d = length(x); H = zeros(d, d)
    for i in 1:d, j in i:d
        ei = zeros(d); ei[i] = h; ej = zeros(d); ej[j] = h
        H[i, j] = (f(x .+ ei .+ ej) - f(x .+ ei .- ej) - f(x .- ei .+ ej) + f(x .- ei .- ej)) / (4h^2)
        H[j, i] = H[i, j]
    end
    return H
end

rng = MersenneTwister(5)
fits = Dict{Symbol, NamedTuple}()
for key in (:exp, :erlang2, :erlang3, :coxian2, :weibull)
    nr = baselines[key].nrates
    best = (Inf, zeros(2 + nr))
    for s in 1:(key == :coxian2 ? 25 : 12)
        θ0 = key == :weibull ? [randn(rng, 2) .* 0.2; log(0.8 + rand(rng)); log(400 + 200 * rand(rng))] :
                               [randn(rng, 2) .* 0.2; log.(10 .^ (rand(rng, nr) .* 1.5 .- 3.3))]
        θ, v = nelder_mead(θ -> negloglik(θ, key), θ0)
        θ, v = nelder_mead(θ -> negloglik(θ, key), θ; step = 0.05)
        θ, v = nelder_mead(θ -> negloglik(θ, key), θ; step = 0.01)
        v < best[1] && (best = (v, θ))
    end
    v, θ = best
    H = numerical_hessian(θ -> negloglik(θ, key), θ)
    @assert all(eigvals(Symmetric(H)) .> 0)
    se = sqrt.(diag(inv(H)))
    fits[key] = (ℓ = -v, θ = θ, se = se, k = 2 + nr)
    @printf("%-12s loglik %.3f  age %.5f (%.5f)  female %.4f (%.4f)  rates %s\n", baselines[key].name, -v, θ[1], se[1], θ[2], se[2],
            string(round.(exp.(θ[3:end]); sigdigits = 5)))
end

# ---- checks against R ----------------------------------------------------------------------------------
ref = Dict(:exp => (-1156.099, -0.01561871, 0.00910568, 0.4809349, 0.1670943), :weibull => (-1147.054, -0.01225703, 0.006957472, 0.3820851, 0.1274768),
           :erlang2 => (-1151.821, -0.01232825, 0.006226864, 0.3806746, 0.1130582), :erlang3 => (-1178.189, -0.01104937, 0.005021693, 0.3431839, 0.0908179),
           :coxian2 => (-1146.937, -0.01135813, 0.0069728, 0.3910823, 0.1298659))
for (key, (ℓr, a, sa, f, sf)) in ref
    F = fits[key]
    @assert isapprox(F.ℓ, ℓr; atol = 5e-3) "loglik $key"
    @assert isapprox(F.θ[1], a; atol = 5e-4) && isapprox(F.θ[2], f; atol = 5e-3) "beta $key"
    @assert isapprox(F.se[1], sa; rtol = 5e-2) && isapprox(F.se[2], sf; rtol = 5e-2) "se $key"
end

CSV.write(joinpath(OUTDIR, "ex_4_4_ph_aft_lung.csv"),
          DataFrame(baseline = [baselines[k].name for k in (:exp, :weibull, :erlang2, :erlang3, :coxian2)],
                    k = [fits[k].k for k in (:exp, :weibull, :erlang2, :erlang3, :coxian2)],
                    beta_age = [fits[k].θ[1] for k in (:exp, :weibull, :erlang2, :erlang3, :coxian2)],
                    se_age = [fits[k].se[1] for k in (:exp, :weibull, :erlang2, :erlang3, :coxian2)],
                    beta_female = [fits[k].θ[2] for k in (:exp, :weibull, :erlang2, :erlang3, :coxian2)],
                    se_female = [fits[k].se[2] for k in (:exp, :weibull, :erlang2, :erlang3, :coxian2)],
                    loglik = [fits[k].ℓ for k in (:exp, :weibull, :erlang2, :erlang3, :coxian2)]))

# ---- predicted survival curves for a man and a woman aged 62, with Kaplan-Meier by sex ---------------------
function S_pred(key, t, x)
    θ = fits[key].θ
    η = dot(x, θ[1:2])
    return baselines[key].fun(t * exp(-η), exp.(θ[3:end]))[1]
end
function km(yy, dd)
    t = sort(unique(yy[dd])); S = 1.0; out = Float64[]
    for u in t
        S *= 1 - count(dd .& (yy .== u)) / count(yy .>= u); push!(out, S)
    end
    return t, out
end
styles = [(:exp, RGB(0.45, 0.45, 0.45), :dashdot, "Exponential"), (:weibull, GREEN, :dot, "Weibull"),
          (:erlang2, ORANGE, :dash, "Erlang(2)"), (:erlang3, RGB(0.50, 0.30, 0.62), :dashdotdot, "Erlang(3)"),
          (:coxian2, RGB(0.72, 0.15, 0.15), :solid, "Coxian(2)")]
tg = range(0, maximum(y); length = 600)
panels = []
for (sx, x, ttl) in ((1, [0.0, 0.0], "(a) Men aged 62"), (2, [0.0, 1.0], "(b) Women aged 62"))
    idx = raw.Sex .== sx
    tt, SS = km(y[idx], δ[idx])
    p = plot([0; tt; maximum(y[idx])], [1.0; SS; SS[end]]; seriestype = :steppost, color = BLUE, lw = 1.3, alpha = 0.7,
             label = "Kaplan-Meier", xlabel = "Time (days)", ylabel = "S(t)", title = ttl, xlims = (0, 1050), ylims = (0, 1.02),
             legend = :topright)
    for (key, col, ls, lab) in styles
        plot!(p, tg, [S_pred(key, t, x) for t in tg]; color = col, lw = 1.8, ls = ls, label = lab)
    end
    push!(panels, p)
end
fig = plot(panels...; layout = (1, 2), size = (1100, 440), margin = 4Plots.mm)
savefig(fig, joinpath(FIGDIR, "ex_4_4_ph_aft_lung.pdf"))
println("Saved ", joinpath(FIGDIR, "ex_4_4_ph_aft_lung.pdf"))
