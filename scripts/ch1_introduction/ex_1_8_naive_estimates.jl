# Section 1.8 (Survival Analysis with Censoring) — the two simple estimates of the survival function for the lung data,
# drawn without the Kaplan-Meier curve (which is only introduced in Section 1.10).
#
# Run from the code_package folder:
#     julia --project=. scripts/ch1_introduction/ex_1_8_naive_estimates.jl
#
# Same data and estimates as in ex_1_8_lung_censoring.jl:
#   A  treat censoring times as event times   Ŝ_A(t) = (1/n) Σ 1{y_i > t}
#   B  discard the censored individuals       Ŝ_B(t) = Σ δ_i 1{y_i > t} / Σ δ_i
#
# Output: figures/ex_1_8_naive_estimates.pdf   (vector)

using Plots, Statistics, RDatasets, Printf

const ORANGE = RGB(0.88, 0.47, 0.22)
const GREEN  = RGB(0.18, 0.42, 0.18)

const FIGDIR = normpath(joinpath(@__DIR__, "..", "..", "figures"))

raw = dataset("survival", "lung")
y = Vector{Float64}(raw.Time)
status = Vector{Int}(raw.Status)          # 1 = censored, 2 = dead
n = length(y)
δ = status .== 2
@assert n == 228 && sum(raw.Time) == 69593 && sum(status) == 393 && count(!, δ) == 63

SA(t) = mean(y .> t)
SB(t) = sum(δ .& (y .> t)) / sum(δ)

median_time(S, grid) = minimum(grid[[S(t) <= 0.5 for t in grid]])
grid = sort(unique(y))
@printf("median: A = %g, B = %g;  t = 300: A = %.2f, B = %.2f;  A at the largest observed time %g: %.3f\n",
        median_time(SA, grid), median_time(SB, grid), SA(300), SB(300), maximum(y), SA(maximum(y)))
@assert median_time(SA, grid) == 252 && median_time(SB, grid) == 226
@assert round(SA(300); digits = 2) == 0.40 && round(SB(300); digits = 2) == 0.39

tg = range(0, maximum(y); length = 1500)
p = plot(tg, SA.(tg); seriestype = :steppost, color = ORANGE, lw = 1.5, ls = :dash,
         label = "A: censoring times as events", xlabel = "Time (days)", ylabel = "S(t)",
         title = "Lung cancer data (n = $n)", ylims = (0, 1.02), size = (600, 420),
         margin = 4Plots.mm, legend = :topright)
plot!(p, tg, SB.(tg); seriestype = :steppost, color = GREEN, lw = 1.5, ls = :dot,
      label = "B: censored individuals discarded")
savefig(p, joinpath(FIGDIR, "ex_1_8_naive_estimates.pdf"))
println("Saved ", joinpath(FIGDIR, "ex_1_8_naive_estimates.pdf"))
