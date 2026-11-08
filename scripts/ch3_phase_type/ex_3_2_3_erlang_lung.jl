# Section 3.2.3 (Erlang distribution) — Erlang(m, λ) with m = 2 and m = 3 fitted to the lung data (right-censored).
#
# Run from the code_package folder:
#     julia --project=. scripts/ch3_phase_type/ex_3_2_3_erlang_lung.jl
#
# Erlang(m, λ):  f(y) = λ^m y^(m-1) e^(-λy) / (m-1)!,   S(y) = e^(-λy) Σ_{j=0}^{m-1} (λy)^j / j!.
# Log-likelihood with right-censoring:  ℓ(λ) = Σ_i [δ_i log f(y_i) + (1 - δ_i) log S(y_i)];  m is fixed, λ is the only parameter,
# so ℓ is maximised by a one-dimensional golden-section search on log λ (no closed form, but no optimisation package needed).
# Checked against R (optimize on the same log-likelihood, using dgamma/pgamma; reference_checks_R/ex_3_2_3_erlang_lung_check.R).
#
# Output: figures/ex_3_2_3_erlang_lung.pdf   (vector);  output/ex_3_2_3_erlang_lung.csv

using Plots, Printf, RDatasets, DataFrames, CSV

const BLUE   = RGB(0.18, 0.37, 0.52)
const ORANGE = RGB(0.88, 0.47, 0.22)
const GREEN  = RGB(0.18, 0.42, 0.18)
const FIGDIR = normpath(joinpath(@__DIR__, "..", "..", "figures"))
const OUTDIR = normpath(joinpath(@__DIR__, "..", "..", "output"))
mkpath(FIGDIR); mkpath(OUTDIR)

raw = dataset("survival", "lung")
y = Vector{Float64}(raw.Time)
δ = raw.Status .== 2
n = length(y)
@assert n == 228 && count(δ) == 165

logfact(k) = sum(log, 1:k; init = 0.0)
logf(t, m, λ) = m * log(λ) + (m - 1) * log(t) - λ * t - logfact(m - 1)
S(t, m, λ) = exp(-λ * t) * sum((λ * t)^j / exp(logfact(j)) for j in 0:m-1)
loglik(m, λ) = sum(δ[i] ? logf(y[i], m, λ) : log(S(y[i], m, λ)) for i in 1:n)

function golden(f, a, b; tol = 1e-13)                # maximise f on [a, b]
    φ = (sqrt(5) - 1) / 2
    c = b - φ * (b - a); d = a + φ * (b - a)
    fc = f(c); fd = f(d)
    while b - a > tol
        if fc > fd
            b, d, fd = d, c, fc; c = b - φ * (b - a); fc = f(c)
        else
            a, c, fc = c, d, fd; d = a + φ * (b - a); fd = f(d)
        end
    end
    return (a + b) / 2
end

fits = Dict{Int, NamedTuple}()
for m in 1:4
    lh = golden(l -> loglik(m, exp(l)), log(1e-6), log(1.0))
    λ = exp(lh); ℓ = loglik(m, λ)
    # median: S(t) = 0.5 by bisection
    lo, hi = 0.0, 5000.0
    for _ in 1:200; mid = (lo + hi) / 2; S(mid, m, λ) > 0.5 ? (lo = mid) : (hi = mid); end
    fits[m] = (λ = λ, ℓ = ℓ, mean = m / λ, median = (lo + hi) / 2)
    @printf("m = %d: lambda_hat = %.7f, mean = %.1f, median = %.1f, loglik = %.3f\n", m, λ, m / λ, (lo + hi) / 2, ℓ)
end
# checks against R
ref = Dict(1 => (0.002370928, -1162.338, 421.7758, 292.3527), 2 => (0.005299698, -1160.196, 377.38, 316.6873),
           3 => (0.008301557, -1188.635, 361.378, 322.1155))
for (m, (λr, ℓr, mr, medr)) in ref
    @assert isapprox(fits[m].λ, λr; rtol = 1e-5) && isapprox(fits[m].ℓ, ℓr; atol = 2e-3)
    @assert isapprox(fits[m].mean, mr; rtol = 1e-4) && isapprox(fits[m].median, medr; atol = 0.05)
end
# the log-likelihood is lower on both sides of the maximum
for m in 2:3
    λ = fits[m].λ
    @assert all(loglik(m, l) < fits[m].ℓ for l in (0.95λ, 0.99λ, 1.01λ, 1.05λ))
end

CSV.write(joinpath(OUTDIR, "ex_3_2_3_erlang_lung.csv"),
          DataFrame(m = 1:4, lambda = [fits[m].λ for m in 1:4], mean = [fits[m].mean for m in 1:4],
                    median = [fits[m].median for m in 1:4], loglik = [fits[m].ℓ for m in 1:4]))

# figure: Kaplan-Meier estimate and the Erlang fits for m = 2 and m = 3
tk = sort(unique(y[δ]))
Skm = cumprod([1 - count(δ .& (y .== t)) / count(y .>= t) for t in tk])
KM(t) = (i = searchsortedlast(tk, t); i == 0 ? 1.0 : Skm[i])
cens = sort(y[.!δ])
tg = range(0, maximum(y); length = 800)
p = plot([0; tk; maximum(y)], [1.0; Skm; Skm[end]]; seriestype = :steppost, color = BLUE, lw = 2.2, label = "Kaplan-Meier",
         xlabel = "Time (days)", ylabel = "S(t)", title = "Lung cancer data (n = $n)",
         xlims = (0, 1050), ylims = (0, 1.02), size = (600, 420), margin = 4Plots.mm, legend = :topright)
scatter!(p, cens, KM.(cens); markershape = :cross, color = BLUE, markersize = 4, markerstrokewidth = 1,
         label = "censored (on Kaplan-Meier curve)")
plot!(p, tg, [S(t, 2, fits[2].λ) for t in tg]; color = ORANGE, lw = 2, ls = :dash, label = "Erlang fit, m = 2")
plot!(p, tg, [S(t, 3, fits[3].λ) for t in tg]; color = GREEN, lw = 2, ls = :dot, label = "Erlang fit, m = 3")
savefig(p, joinpath(FIGDIR, "ex_3_2_3_erlang_lung.pdf"))
println("Saved ", joinpath(FIGDIR, "ex_3_2_3_erlang_lung.pdf"))
