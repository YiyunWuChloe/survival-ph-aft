# Section 3.2.2 (Hyperexponential distribution) — fitted to the lung data (right-censored), by the EM algorithm.
#
# Run from the code_package folder:
#     julia --project=. scripts/ch3_phase_type/ex_3_2_2_hyperexponential_lung.jl
#
# Model: with probability α_j the individual comes from component j (j = 1, 2) with exponential rate λ_j.
#   f(y) = Σ α_j λ_j e^{-λ_j y},  S(y) = Σ α_j e^{-λ_j y};  log-likelihood  Σ_i [δ_i log f(y_i) + (1 - δ_i) log S(y_i)].
# EM: E-step (responsibilities)  r_ij = α_j f_j(y_i)^δ_i S_j(y_i)^(1-δ_i) / Σ_k α_k f_k(y_i)^δ_i S_k(y_i)^(1-δ_i);
#     M-step                      α_j = mean_i r_ij,   λ_j = Σ_i δ_i r_ij / Σ_i y_i r_ij.
# The log-likelihood never decreases in an EM step. 
# Output: figures/ex_3_2_2_hyperexponential_lung.pdf   (vector: Kaplan-Meier estimate and the fitted survival function)
# Checked against a multi-start Nelder-Mead maximisation of the same log-likelihood in R (reference_checks_R/ex_3_2_2_hyperexponential_lung_check.R):
#   the maximum is -1162.338 with λ_1 = λ_2, that is, the exponential fit.

using Plots, Printf, Random, RDatasets

raw = dataset("survival", "lung")
y = Vector{Float64}(raw.Time)
δ = Float64.(raw.Status .== 2)
n = length(y)
@assert n == 228 && sum(δ) == 165

loglik(α, λ) = sum(δ[i] * log(sum(α[j] * λ[j] * exp(-λ[j] * y[i]) for j in 1:2)) +
                   (1 - δ[i]) * log(sum(α[j] * exp(-λ[j] * y[i]) for j in 1:2)) for i in 1:n)

function em(α, λ; iters = 200_000, tol = 1e-12)
    ℓold = loglik(α, λ)
    for it in 1:iters
        r = zeros(n, 2)
        for i in 1:n, j in 1:2
            f = λ[j] * exp(-λ[j] * y[i]); S = exp(-λ[j] * y[i])
            r[i, j] = α[j] * (δ[i] == 1 ? f : S)
        end
        r ./= sum(r; dims = 2)
        α = vec(sum(r; dims = 1)) ./ n
        λ = [sum(δ .* r[:, j]) / sum(y .* r[:, j]) for j in 1:2]
        ℓnew = loglik(α, λ)
        @assert ℓnew >= ℓold - 1e-9                 # EM never decreases the log-likelihood
        abs(ℓnew - ℓold) < tol && return α, λ, ℓnew, it
        ℓold = ℓnew
    end
    return α, λ, ℓold, iters
end

ℓexp = 165 * log(165 / sum(y)) - 165
rng = MersenneTwister(1)
best = (zeros(2), zeros(2), -Inf, 0)
for s in 1:30
    α0 = [rand(rng, 0.05:0.05:0.95)]; α0 = [α0[1], 1 - α0[1]]
    λ0 = sort(exp.(rand(rng, 2) .* (log(1e-2) - log(1e-4)) .+ log(1e-4)))
    α, λ, ℓ, it = em(α0, λ0)
    if ℓ > best[3]
        global best = (α, λ, ℓ, it)
    end
end
α, λ, ℓ, it = best
@printf("best log-likelihood: %.4f  (exponential: %.4f)\n", ℓ, ℓexp)
@printf("alpha = (%.4f, %.4f), lambda = (%.6f, %.6f)\n", α..., λ...)
@assert isapprox(ℓ, -1162.338; atol = 1e-3) && isapprox(ℓexp, -1162.338; atol = 1e-3)
@assert abs(λ[1] - λ[2]) < 1e-5                    # the two rates coincide: the fit is the exponential distribution

# ---- figure: Kaplan-Meier estimate and fitted survival function ----------------------------------------
const BLUE   = RGB(0.18, 0.37, 0.52)
const ORANGE = RGB(0.88, 0.47, 0.22)
const FIGDIR = normpath(joinpath(@__DIR__, "..", "..", "figures"))
mkpath(FIGDIR)
δb = δ .== 1
tk = sort(unique(y[δb]))
Skm = cumprod([1 - count(δb .& (y .== t)) / count(y .>= t) for t in tk])
KM(t) = (i = searchsortedlast(tk, t); i == 0 ? 1.0 : Skm[i])
Sfit(t) = sum(α[j] * exp(-λ[j] * t) for j in 1:2)
cens = sort(y[.!δb])
tg = range(0, maximum(y); length = 800)
p = plot([0; tk; maximum(y)], [1.0; Skm; Skm[end]]; seriestype = :steppost, color = BLUE, lw = 2.2, label = "Kaplan-Meier",
         xlabel = "Time (days)", ylabel = "S(t)", title = "Lung cancer data (n = $n)",
         xlims = (0, 1050), ylims = (0, 1.02), size = (600, 420), margin = 4Plots.mm, legend = :topright)
scatter!(p, cens, KM.(cens); markershape = :cross, color = BLUE, markersize = 4, markerstrokewidth = 1,
         label = "censored (on Kaplan-Meier curve)")
plot!(p, tg, Sfit.(tg); color = ORANGE, lw = 2, ls = :dash, label = "Hyperexponential fit")
savefig(p, joinpath(FIGDIR, "ex_3_2_2_hyperexponential_lung.pdf"))
println("Saved ", joinpath(FIGDIR, "ex_3_2_2_hyperexponential_lung.pdf"))
