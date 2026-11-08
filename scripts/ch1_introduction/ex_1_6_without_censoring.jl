# Section 1.6 (Survival Analysis Without Censoring) — real data example (aircondit7)
#
# Run from the code_package folder:
#     julia --project=. scripts/ch1_introduction/ex_1_6_without_censoring.jl
#
# Data: `aircondit7` from the R package `boot` (read through RDatasets.jl): 24 times
# (hours) between failures of air-conditioning equipment. All 24 times are observed,
# so there is no censoring. The survival function S(t) = P(T > t) is estimated by the
# empirical survival function Ŝ(t) = (1/n) Σ 1{t_i > t}.
#
# Output: output/aircondit7.csv                 the data (id, hours), as used here
#         output/ex_1_6_table.csv               ordered times and Ŝ
#         figures/ex_1_6_without_censoring.pdf 

using CSV, DataFrames, Plots, Printf, Statistics

# colour palette shared by all figures of the report
const BLUE   = RGB(0.18, 0.37, 0.52)
const ORANGE = RGB(0.88, 0.47, 0.22)
const GREEN  = RGB(0.18, 0.42, 0.18)
const LBLUE  = RGB(0.74, 0.84, 0.92)

using RDatasets

const FIGDIR = normpath(joinpath(@__DIR__, "..", "..", "figures"))
const OUTDIR = normpath(joinpath(@__DIR__, "..", "..", "output"))
mkpath(FIGDIR); mkpath(OUTDIR)

# --- data -----------------------------------------------------------------------
raw = dataset("boot", "aircondit7")
@assert names(raw) == ["Hours"]
t = Vector{Float64}(raw.Hours)
n = length(t)
@assert n == 24 && all(t .> 0)
# independent check against the values printed by R's boot::aircondit7$hours
@assert Int.(t) == [3, 5, 5, 13, 14, 15, 22, 22, 23, 30, 36, 39, 44, 46, 50, 72, 79, 88,
                    97, 102, 139, 188, 197, 210]
CSV.write(joinpath(OUTDIR, "aircondit7.csv"), DataFrame(id = 1:n, hours = Int.(t)))

# --- summary --------------------------------------------------------------------
ts = sort(t)
ties = [(Int(v), count(==(v), ts)) for v in unique(ts) if count(==(v), ts) > 1]
@printf("n = %d, min = %g, max = %g, mean = %.2f, median = %.1f, sd = %.2f\n",
        n, minimum(t), maximum(t), mean(t), median(t), std(t))
@printf("quartiles = %.2f, %.2f\n", quantile(t, 0.25), quantile(t, 0.75))
println("tied values (value, count): ", ties)

# --- empirical survival function --------------------------------------------------
Shat(x) = mean(t .> x)
tab = DataFrame(i = 1:n, hours = Int.(ts), Shat = [Shat(x) for x in ts])
CSV.write(joinpath(OUTDIR, "ex_1_6_table.csv"), tab)
@assert Shat(0.0) == 1.0 && Shat(maximum(t)) == 0.0
@assert issorted(tab.Shat; rev = true)
@assert isapprox(Shat(3.0), 23 / 24)          # one time at 3: drops by 1/n
@assert isapprox(Shat(5.0), 21 / 24)          # two tied times at 5: drops by 2/n

# median survival time: smallest t with Ŝ(t) ≤ 0.5; a reading at t = 100
t_med = minimum(ts[[Shat(x) <= 0.5 for x in ts]])
@printf("median survival time (smallest t with Ŝ(t) <= 0.5) = %g, Ŝ(%g) = %.3f\n",
        t_med, t_med, Shat(t_med))
@printf("Ŝ(100) = %.3f  (%d of %d times exceed 100)\n", Shat(100.0), count(>(100), t), n)

# --- figure -------------------------------------------------------------------------
# Ŝ is drawn as horizontal segments only (a function has no vertical lines). At each distinct event
# time an open circle marks the value just before it, Ŝ(t-), and a filled dot marks the value at and
# just after it, Ŝ(t) = Ŝ(t+).
ut = unique(ts)                                          # distinct ordered times
before = [mean(t .>= s) for s in ut]                     # Ŝ(s-)
after  = [Shat(s) for s in ut]                           # Ŝ(s) = Ŝ(s+)
edges  = [0.0; ut; 220.0]
levels = [1.0; after]
fig = plot(; xlabel = "Time between failures (hours)", ylabel = "S(t)",
           title = "Air-conditioning failures (n = $n)", xlims = (0, 220), ylims = (0, 1.02),
           size = (560, 400), margin = 4Plots.mm, legend = :topright)
for k in eachindex(levels)
    plot!(fig, [edges[k], edges[k + 1]], [levels[k], levels[k]]; color = BLUE, lw = 1.5,
          label = k == 1 ? "empirical Ŝ(t)" : "")
end
scatter!(fig, ut, before; markercolor = :white, markerstrokecolor = BLUE, markerstrokewidth = 1,
         markersize = 3, label = "")
scatter!(fig, ut, after; markercolor = BLUE, markerstrokecolor = BLUE, markersize = 3, label = "")
plot!(fig, [0, t_med, t_med], [0.5, 0.5, 0]; color = ORANGE, lw = 1.5, ls = :dash,
      label = "median: Ŝ(t) = 0.5")
plot!(fig, [0, 100, 100], [Shat(100.0), Shat(100.0), 0]; color = GREEN, lw = 1.5, ls = :dot,
      label = "Ŝ(100)")
savefig(fig, joinpath(FIGDIR, "ex_1_6_without_censoring.pdf"))
println("Saved ", joinpath(FIGDIR, "ex_1_6_without_censoring.pdf"))
