# Section 3.2.4 (Coxian Distribution) — simulation: the Coxian distribution as a phase-type distribution
#
# Run from the code_package folder:
#     julia --project=. scripts/ch3_phase_type/ex_3_2_4_coxian_simulation.jl
#
# Coxian with m = 3 phases (report notation): from phase j the process moves to
# phase j+1 at rate λ_j (j = 1, ..., m-1) or is absorbed at rate µ_j (j = 1, ..., m).
# There is no closed form for f and S, so the theoretical curves are the matrix
# formulas f(t) = α exp(Tt) d and S(t) = α exp(Tt) 1, computed by
# PhaseTypeSurvival's `pdf` and `ccdf`.
#
# Output: figures/ex_3_2_4_coxian_simulation.pdf, output/ex_3_2_4_coxian_simulation_sample.csv

using CSV, DataFrames, LinearAlgebra, Plots, Printf, Random, Statistics

# colour palette shared by all figures of the report
const BLUE   = RGB(0.18, 0.37, 0.52)
const ORANGE = RGB(0.88, 0.47, 0.22)
const GREEN  = RGB(0.18, 0.42, 0.18)
const LBLUE  = RGB(0.74, 0.84, 0.92)

using PhaseTypeSurvival

const SEED = 2026
const N = 200
const λ = [3.0, 2.0]          # λ_1, ..., λ_{m-1}: rates to the next phase
const μ = [0.5, 0.5, 1.5]     # µ_1, ..., µ_m: rates to absorption
const M = length(μ)
const BINWIDTH = 0.25

const FIGDIR = normpath(joinpath(@__DIR__, "..", "..", "figures"))
const OUTDIR = normpath(joinpath(@__DIR__, "..", "..", "output"))
mkpath(FIGDIR); mkpath(OUTDIR)

# --- representation -------------------------------------------------------------
# The package's coxian(rates, q): phase i lasts Exp(rates[i]) and then moves on with
# probability q[i], otherwise is absorbed. Thesis rates (λ_j, µ_j) map to
#   rates[j] = λ_j + µ_j,  q[j] = λ_j / (λ_j + µ_j),  rates[m] = µ_m.
rates = [λ .+ μ[1:end-1]; μ[end]]
q = λ ./ (λ .+ μ[1:end-1])
d = coxian(rates, q)

# sanity check 1: the representation equals the Coxian matrices of the report
T_report = Matrix(Diagonal(-[λ .+ μ[1:end-1]; μ[end]])) + Matrix(diagm(1 => λ))
@assert isapprox(subintensity(d), T_report; atol = 1e-12)
@assert d.α == [1.0; zeros(M - 1)]
@assert isapprox(exitrates(d), μ; atol = 1e-12)                  # (eq:coxian-d)

# sanity check 2 (no closed form available): f integrates to 1, S(t) = ∫_t^∞ f(s) ds,
# and the mean matches α(-T)^{-1}1. Trapezoidal rule on a fine grid.
trapz(x, y) = sum((x[i+1] - x[i]) * (y[i] + y[i+1]) / 2 for i in 1:length(x)-1)
sg = range(0, 40; length = 40001)
fg = pdf.(Ref(d), sg)
@assert isapprox(trapz(sg, fg), 1.0; atol = 1e-6)
@assert isapprox(trapz(sg, sg .* fg), mean(d); atol = 1e-5)
for t in (0.5, 1.0, 2.0)
    k = sg .>= t
    @assert isapprox(trapz(sg[k], fg[k]), ccdf(d, t); atol = 1e-5)
end
@assert isapprox(ccdf(d, 0.0), 1.0; atol = 1e-12)

# --- simulate: absorption times of the Markov chain, no censoring --------------
rng = MersenneTwister(SEED)
y = rand(rng, d, N)
@assert length(y) == N && all(y .> 0)
CSV.write(joinpath(OUTDIR, "ex_3_2_4_coxian_simulation_sample.csv"), DataFrame(time = y))

@printf("n = %d, m = %d, λ = %s, µ = %s\n", N, M, λ, μ)
@printf("sample mean = %.3f (theory %.3f)\n", mean(y), mean(d))

# --- Kaplan-Meier estimate (no censoring: 1 - empirical CDF) -------------------
ys = sort(y)
km_t = [0.0; ys]
km_S = [1.0; 1 .- (1:N) ./ N]
@assert km_S[end] == 0.0 && issorted(km_S; rev = true)
D = maximum(max(abs((1 - i / N) - ccdf(d, ys[i])), abs((1 - (i - 1) / N) - ccdf(d, ys[i]))) for i in 1:N)
@printf("max |KM - S| = %.3f\n", D)

# --- figure --------------------------------------------------------------------
tg = range(0, maximum(y); length = 301)
edges = 0:BINWIDTH:(ceil(maximum(y) / BINWIDTH) * BINWIDTH)

p1 = histogram(y; bins = edges, normalize = :pdf, color = LBLUE, linecolor = BLUE,
               label = "simulated data (n = $N)", xlabel = "t", ylabel = "f(t)",
               title = "Density")
plot!(p1, tg, pdf.(Ref(d), tg); color = BLUE, lw = 2.5, label = "theoretical")

p2 = plot(tg, ccdf.(Ref(d), tg); color = ORANGE, lw = 2.5, label = "theoretical",
          xlabel = "t", ylabel = "S(t)", title = "Survival function")
plot!(p2, km_t, km_S; seriestype = :steppost, color = BLUE, lw = 1.5, ls = :dash,
      label = "Kaplan-Meier")

fig = plot(p1, p2; layout = (1, 2), size = (900, 350), margin = 4Plots.mm)
savefig(fig, joinpath(FIGDIR, "ex_3_2_4_coxian_simulation.pdf"))
println("Saved ", joinpath(FIGDIR, "ex_3_2_4_coxian_simulation.pdf"))
