# Section 3.2.4 (Coxian distribution) — Coxian distribution with m = 2 phases fitted to the lung data (right-censored).
#
# Run from the code_package folder:
#     julia --project=. scripts/ch3_phase_type/ex_3_2_4_coxian_lung.jl
#
# Coxian, m = 2:  start in phase 1; from phase 1 move to phase 2 at rate λ1 or get absorbed at rate μ1; from phase 2 get absorbed at rate μ2.
#   α = (1, 0),  T = [-(λ1+μ1)  λ1; 0  -μ2],  d = (μ1, μ2)ᵀ,   S(y) = α e^{Ty} 1,   f(y) = α e^{Ty} d.
# Log-likelihood with right-censoring:  ℓ = Σ_i [δ_i log f(y_i) + (1 - δ_i) log S(y_i)]  (three parameters).
# No closed-form maximiser: ℓ is maximised numerically (own Nelder-Mead on the log-rates, 30 random starting values).
# Checked against R (closed-form S and f of the 2-phase Coxian, 300 random starts of optim/Nelder-Mead; reference_checks_R/ex_3_2_4_coxian_lung_check.R)
# and against the matrix-exponential formulas below.
#
# Output: figures/ex_3_2_4_coxian_lung.pdf   (vector);  output/ex_3_2_4_coxian_lung.csv

using Plots, Printf, Random, LinearAlgebra, RDatasets, DataFrames, CSV

const BLUE   = RGB(0.18, 0.37, 0.52)
const ORANGE = RGB(0.88, 0.47, 0.22)
const FIGDIR = normpath(joinpath(@__DIR__, "..", "..", "figures"))
const OUTDIR = normpath(joinpath(@__DIR__, "..", "..", "output"))
mkpath(FIGDIR); mkpath(OUTDIR)

raw = dataset("survival", "lung")
y = Vector{Float64}(raw.Time)
δ = raw.Status .== 2
n = length(y)
@assert n == 228 && count(δ) == 165

Tmat(λ1, μ1, μ2) = [-(λ1 + μ1) λ1; 0.0 -μ2]
S_f(t, p) = (E = exp(Tmat(p...) * t); E[1, 1] + E[1, 2])                       # α e^{Tt} 1
f_f(t, p) = (E = exp(Tmat(p...) * t); E[1, 1] * p[2] + E[1, 2] * p[3])           # α e^{Tt} d
loglik(p) = sum(δ[i] ? log(max(f_f(y[i], p), 1e-300)) : log(max(S_f(y[i], p), 1e-300)) for i in 1:n)
negll(θ) = -loglik(exp.(θ))

# closed form for a ≠ b (a = λ1 + μ1, b = μ2) as a check of the matrix-exponential formulas
function S_closed(t, (λ1, μ1, μ2))
    a = λ1 + μ1; b = μ2
    return exp(-a * t) + λ1 * (exp(-b * t) - exp(-a * t)) / (a - b)
end
let p = (0.004, 0.001, 0.006)
    for t in (50.0, 200.0, 500.0)
        @assert isapprox(S_f(t, p), S_closed(t, p); atol = 1e-12)
        h = 1e-5
        @assert isapprox(f_f(t, p), (S_f(t - h, p) - S_f(t + h, p)) / (2h); rtol = 1e-6)    # f = -S'
    end
end

# Nelder-Mead (minimisation)
function nelder_mead(f, x0; step = 0.5, maxit = 20_000, tol = 1e-13)
    d = length(x0)
    X = [copy(x0) for _ in 1:d+1]
    for i in 1:d; X[i+1][i] += step; end
    F = [f(x) for x in X]
    for it in 1:maxit
        idx = sortperm(F); X = X[idx]; F = F[idx]
        abs(F[end] - F[1]) < tol * (1 + abs(F[1])) && break
        c = sum(X[1:d]) ./ d
        xr = c .+ (c .- X[end]); fr = f(xr)
        if fr < F[1]
            xe = c .+ 2 .* (c .- X[end]); fe = f(xe)
            fe < fr ? (X[end] = xe; F[end] = fe) : (X[end] = xr; F[end] = fr)
        elseif fr < F[d]
            X[end] = xr; F[end] = fr
        else
            xc = fr < F[end] ? c .+ 0.5 .* (xr .- c) : c .+ 0.5 .* (X[end] .- c); fc = f(xc)
            if fc < min(fr, F[end])
                X[end] = xc; F[end] = fc
            else
                for i in 2:d+1; X[i] = X[1] .+ 0.5 .* (X[i] .- X[1]); F[i] = f(X[i]); end
            end
        end
    end
    i = argmin(F)
    return X[i], F[i]
end

rng = MersenneTwister(2)
results = Tuple{Float64, Vector{Float64}}[]
for s in 1:30
    θ0 = log.(10 .^ (rand(rng, 3) .* 2.5 .- 4))
    θ, v = nelder_mead(negll, θ0)
    θ, v = nelder_mead(negll, θ; step = 0.05)                  # restart once to refine
    push!(results, (-v, exp.(θ)))
end
sort!(results; by = r -> -r[1])
ℓ, p = results[1]
nbest = count(r -> r[1] > ℓ - 1e-3, results)
@printf("best log-likelihood %.4f; %d of 30 starts within 1e-3 of it\n", ℓ, nbest)
λ1, μ1, μ2 = p
@printf("lambda1 = %.6f, mu1 = %.6f, mu2 = %.6f;  λ1+μ1 = %.6f\n", λ1, μ1, μ2, λ1 + μ1)
# checks against R
@assert isapprox(ℓ, -1153.593; atol = 2e-3)
@assert isapprox(λ1, 0.003681629; rtol = 2e-3) && isapprox(μ1, 0.0009547780; rtol = 2e-3) && isapprox(μ2, 0.004636407; rtol = 2e-3)

a = λ1 + μ1
mean_fit = 1 / a + (λ1 / a) / μ2
function median_of(p)
    lo, hi = 0.0, 5000.0
    for _ in 1:200
        mid = (lo + hi) / 2
        S_f(mid, p) > 0.5 ? (lo = mid) : (hi = mid)
    end
    return (lo + hi) / 2
end
median_fit = median_of(p)
@printf("mean = %.1f days, median = %.1f days, P(absorbed from phase 1) = %.3f\n", mean_fit, median_fit, μ1 / a)
CSV.write(joinpath(OUTDIR, "ex_3_2_4_coxian_lung.csv"),
          DataFrame(lambda1 = λ1, mu1 = μ1, mu2 = μ2, mean = mean_fit, median = median_fit, loglik = ℓ))

# figure
tk = sort(unique(y[δ]))
Skm = cumprod([1 - count(δ .& (y .== t)) / count(y .>= t) for t in tk])
KM(t) = (i = searchsortedlast(tk, t); i == 0 ? 1.0 : Skm[i])
cens = sort(y[.!δ])
tg = range(0, maximum(y); length = 800)
fig = plot([0; tk; maximum(y)], [1.0; Skm; Skm[end]]; seriestype = :steppost, color = BLUE, lw = 2.2, label = "Kaplan-Meier",
           xlabel = "Time (days)", ylabel = "S(t)", title = "Lung cancer data (n = $n)",
           xlims = (0, 1050), ylims = (0, 1.02), size = (600, 420), margin = 4Plots.mm, legend = :topright)
scatter!(fig, cens, KM.(cens); markershape = :cross, color = BLUE, markersize = 4, markerstrokewidth = 1,
         label = "censored (on Kaplan-Meier curve)")
plot!(fig, tg, [S_f(t, p) for t in tg]; color = ORANGE, lw = 2, ls = :dash, label = "Coxian fit, m = 2")
savefig(fig, joinpath(FIGDIR, "ex_3_2_4_coxian_lung.pdf"))
println("Saved ", joinpath(FIGDIR, "ex_3_2_4_coxian_lung.pdf"))
