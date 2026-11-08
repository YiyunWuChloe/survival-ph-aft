# Section 3.2.1 (Exponential Distribution) — simulation: the exponential distribution as a phase-type distribution
#
# Run from the code_package folder:
#     julia --project=. scripts/ch3_phase_type/ex_3_2_1_exponential_simulation.jl
#
# Exp(λ) is PH(α, T) with m = 1: α = 1, T = -λ, d = λ.
# One simulated sample (no censoring) is compared with the theoretical
# density and survival function.
#
# Output: figures/ex_3_2_1_exponential_simulation.pdf, output/ex_3_2_1_exponential_simulation_sample.csv

using CSV, DataFrames, Plots, Printf, Random, Statistics

# colour palette shared by all figures of the report
const BLUE   = RGB(0.18, 0.37, 0.52)
const ORANGE = RGB(0.88, 0.47, 0.22)
const GREEN  = RGB(0.18, 0.42, 0.18)
const LBLUE  = RGB(0.74, 0.84, 0.92)

using PhaseTypeSurvival

const SEED = 2026
const N = 200
const λ = 1.0
const BINWIDTH = 0.25

const FIGDIR = normpath(joinpath(@__DIR__, "..", "..", "figures"))
const OUTDIR = normpath(joinpath(@__DIR__, "..", "..", "output"))
mkpath(FIGDIR); mkpath(OUTDIR)

# --- representation: PH(α, T) with m = 1 --------------------------------------
d = PhaseType([1.0], reshape([-λ], 1, 1))

# sanity check: the matrix formulas (eq:ph-survival, eq:ph-density) reduce to the exponential formulas (eq:exp-density, eq:exp-survival-hazard)
for t in (0.0, 0.5, 1.0, 2.5, 5.0)
    @assert isapprox(pdf(d, t), λ * exp(-λ * t); atol = 1e-10)
    @assert isapprox(ccdf(d, t), exp(-λ * t); atol = 1e-10)
    @assert isapprox(cdf(d, t), 1 - exp(-λ * t); atol = 1e-10)
end
@assert exitrates(d) == [λ]

# --- simulate: absorption times of the Markov chain, no censoring --------------
rng = MersenneTwister(SEED)
y = rand(rng, d, N)
@assert length(y) == N && all(y .> 0)
CSV.write(joinpath(OUTDIR, "ex_3_2_1_exponential_simulation_sample.csv"), DataFrame(time = y))

@printf("n = %d, true λ = %.2f, sample mean = %.3f, λ̂ = 1/mean = %.3f\n",
        N, λ, mean(y), 1 / mean(y))

# --- Kaplan-Meier estimate (no censoring: 1 - empirical CDF) -------------------
ys = sort(y)
km_t = [0.0; ys]
km_S = [1.0; 1 .- (1:N) ./ N]
@assert km_S[end] == 0.0 && issorted(km_S; rev = true)

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
savefig(fig, joinpath(FIGDIR, "ex_3_2_1_exponential_simulation.pdf"))
println("Saved ", joinpath(FIGDIR, "ex_3_2_1_exponential_simulation.pdf"))
