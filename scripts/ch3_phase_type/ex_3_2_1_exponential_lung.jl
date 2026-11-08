# Section 3.2.1 (Exponential distribution) — the exponential distribution fitted to the lung data, which have right-censoring.
#
# Run from the code_package folder:
#     julia --project=. scripts/ch3_phase_type/ex_3_2_1_exponential_lung.jl
#
# Likelihood with right-censoring: an observed death contributes the density f(y) = λ exp(-λ y), a censored observation
# contributes the survival probability S(y) = exp(-λ y), so   L(λ) = λ^d exp(-λ Σ y_i),   d = number of deaths.
# log L = d log λ - λ Σ y_i,   λ_hat = d / Σ y_i  (closed form, no optimisation),   se(λ_hat) = λ_hat / sqrt(d).
# Checked against R: survreg(Surv(time, status == 2) ~ 1, dist = "exponential").
#
# Output: figures/ex_3_2_1_exponential_lung.pdf   (vector)

using Plots, Printf, RDatasets

const BLUE   = RGB(0.18, 0.37, 0.52)
const ORANGE = RGB(0.88, 0.47, 0.22)
const FIGDIR = normpath(joinpath(@__DIR__, "..", "..", "figures"))
mkpath(FIGDIR)

raw = dataset("survival", "lung")
y = Vector{Float64}(raw.Time)
δ = raw.Status .== 2
n = length(y); d = count(δ); T = sum(y)
@assert n == 228 && d == 165 && T == 69593

λ = d / T
se = λ / sqrt(d)
ℓ = d * log(λ) - λ * T
@printf("lambda_hat = %.6f per day, se = %.6f, mean = %.2f days, median = %.2f days, loglik = %.3f\n", λ, se, 1 / λ, log(2) / λ, ℓ)
@assert isapprox(λ, 0.002370928; atol = 1e-9) && isapprox(1 / λ, 421.7758; atol = 1e-3)
@assert isapprox(log(2) / λ, 292.3527; atol = 1e-3) && isapprox(ℓ, -1162.338; atol = 1e-3) && isapprox(se, 0.0001845765; atol = 1e-9)
# check that λ_hat maximises the log-likelihood (the derivative is zero and the log-likelihood is lower on both sides)
@assert d / λ - T ≈ 0 atol = 1e-8
@assert all(d * log(l) - l * T < ℓ for l in (0.9λ, 0.99λ, 1.01λ, 1.1λ))
ci = (λ - 1.96se, λ + 1.96se)
@printf("95%% CI for lambda (normal approximation): (%.6f, %.6f)\n", ci...)

# Kaplan-Meier estimate of the survival function of the whole data set
tk = sort(unique(y[δ]))
S = cumprod([1 - count(δ .& (y .== t)) / count(y .>= t) for t in tk])
KM(t) = (i = searchsortedlast(tk, t); i == 0 ? 1.0 : S[i])
Sfit(t) = exp(-λ * t)
@assert isapprox(KM(300), 0.5306081; atol = 1e-6)
# where does the fitted curve cross the Kaplan-Meier curve?
cross = tk[findfirst(t -> t > 100 && Sfit(t) >= KM(t), tk)]
@printf("the fitted curve is below the Kaplan-Meier curve up to %g days and above it afterwards (checked at event times)\n", cross)
@printf("t = 100: KM %.3f, fit %.3f; t = 300: KM %.3f, fit %.3f; t = 600: KM %.3f, fit %.3f\n", KM(100), Sfit(100), KM(300), Sfit(300), KM(600), Sfit(600))

tg = range(0, maximum(y); length = 800)
cens = sort(y[.!δ])
p = plot([0; tk; maximum(y)], [1.0; S; S[end]]; seriestype = :steppost, color = BLUE, lw = 2.2, label = "Kaplan-Meier",
         xlabel = "Time (days)", ylabel = "S(t)", title = "Lung cancer data (n = $n)",
         xlims = (0, 1050), ylims = (0, 1.02), size = (600, 420), margin = 4Plots.mm, legend = :topright)
scatter!(p, cens, KM.(cens); markershape = :cross, color = BLUE, markersize = 4, markerstrokewidth = 1,
         label = "censored (on Kaplan-Meier curve)")
plot!(p, tg, Sfit.(tg); color = ORANGE, lw = 2, ls = :dash, label = "Exponential fit")
savefig(p, joinpath(FIGDIR, "ex_3_2_1_exponential_lung.pdf"))
println("Saved ", joinpath(FIGDIR, "ex_3_2_1_exponential_lung.pdf"))
