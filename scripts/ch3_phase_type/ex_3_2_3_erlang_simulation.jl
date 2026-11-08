# Section 3.2.3 (Erlang Distribution) — simulation: the Erlang distribution as a phase-type distribution
#
# Run from the code_package folder:
#     julia --project=. scripts/ch3_phase_type/ex_3_2_3_erlang_simulation.jl
#
# Erlang E_m(λ): m phases in series, each left at rate λ; α = (1, 0, ..., 0),
# T has -λ on the diagonal and λ on the superdiagonal.
# One simulated sample (no censoring) is compared with the theoretical density
# (eq:erlang-density) and survival function (eq:erlang-survival).
#
# Output: figures/ex_3_2_3_erlang_simulation.pdf, output/ex_3_2_3_erlang_simulation_sample.csv

using CSV, DataFrames, LinearAlgebra, Plots, Printf, Random, Statistics

# colour palette shared by all figures of the report
const BLUE   = RGB(0.18, 0.37, 0.52)
const ORANGE = RGB(0.88, 0.47, 0.22)
const GREEN  = RGB(0.18, 0.42, 0.18)
const LBLUE  = RGB(0.74, 0.84, 0.92)

using PhaseTypeSurvival

const SEED = 2026
const N = 200
const M = 3
const λ = 2.0
const BINWIDTH = 0.25

const FIGDIR = normpath(joinpath(@__DIR__, "..", "..", "figures"))
const OUTDIR = normpath(joinpath(@__DIR__, "..", "..", "output"))
mkpath(FIGDIR); mkpath(OUTDIR)

# --- representation: PH(α, T) --------------------------------------------------
α = [1.0; zeros(M - 1)]
T = Matrix(Diagonal(fill(-λ, M))) + Matrix(diagm(1 => fill(λ, M - 1)))
d = PhaseType(α, T)

f_true(t) = λ^M * t^(M - 1) * exp(-λ * t) / factorial(M - 1)          # (eq:erlang-density)
S_true(t) = exp(-λ * t) * sum((λ * t)^k / factorial(k) for k in 0:(M - 1))   # (eq:erlang-survival)

# sanity check: matrix formulas reduce to the Erlang formulas
for t in (0.0, 0.5, 1.0, 2.5, 5.0)
    @assert isapprox(pdf(d, t), f_true(t); atol = 1e-10)
    @assert isapprox(ccdf(d, t), S_true(t); atol = 1e-10)
end
@assert exitrates(d) == [zeros(M - 1); λ]
@assert isapprox(mean(d), M / λ; atol = 1e-10)
@assert isapprox(var(d), M / λ^2; atol = 1e-10)

# --- simulate: absorption times of the Markov chain, no censoring --------------
rng = MersenneTwister(SEED)
y = rand(rng, d, N)
@assert length(y) == N && all(y .> 0)
CSV.write(joinpath(OUTDIR, "ex_3_2_3_erlang_simulation_sample.csv"), DataFrame(time = y))

@printf("n = %d, m = %d, λ = %.1f\n", N, M, λ)
@printf("sample mean = %.3f (theory %.3f)\n", mean(y), M / λ)

# --- Kaplan-Meier estimate (no censoring: 1 - empirical CDF) -------------------
ys = sort(y)
km_t = [0.0; ys]
km_S = [1.0; 1 .- (1:N) ./ N]
@assert km_S[end] == 0.0 && issorted(km_S; rev = true)
D = maximum(max(abs((1 - i / N) - S_true(ys[i])), abs((1 - (i - 1) / N) - S_true(ys[i]))) for i in 1:N)
@printf("max |KM - S| = %.3f\n", D)

# --- figure --------------------------------------------------------------------
tg = range(0, maximum(y); length = 301)
edges = 0:BINWIDTH:(ceil(maximum(y) / BINWIDTH) * BINWIDTH)

p1 = histogram(y; bins = edges, normalize = :pdf, color = LBLUE, linecolor = BLUE,
               label = "simulated data (n = $N)", xlabel = "t", ylabel = "f(t)",
               title = "Density")
plot!(p1, tg, f_true.(tg); color = BLUE, lw = 2.5, label = "theoretical")

p2 = plot(tg, S_true.(tg); color = ORANGE, lw = 2.5, label = "theoretical",
          xlabel = "t", ylabel = "S(t)", title = "Survival function")
plot!(p2, km_t, km_S; seriestype = :steppost, color = BLUE, lw = 1.5, ls = :dash,
      label = "Kaplan-Meier")

fig = plot(p1, p2; layout = (1, 2), size = (900, 350), margin = 4Plots.mm)
savefig(fig, joinpath(FIGDIR, "ex_3_2_3_erlang_simulation.pdf"))
println("Saved ", joinpath(FIGDIR, "ex_3_2_3_erlang_simulation.pdf"))
