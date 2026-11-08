# Sections 1.8 and 1.10 (Survival Analysis with Censoring; Kaplan-Meier Estimator) — lung data: timeline (Section 1.8) and the two simple estimates with the Kaplan-Meier estimate (Section 1.10)
#
# Run from the code_package folder:
#     julia --project=. scripts/ch1_introduction/ex_1_8_lung_censoring.jl
#
# Data: `lung` from the R package `survival` (read through RDatasets.jl): survival
# times in days of 228 patients with advanced lung cancer (NCCTG). status = 2 means
# death (event observed), status = 1 means censored.
#
# Three estimates of the survival function are compared:
#   A  treat censoring times as event times          Ŝ_A(t) = (1/n) Σ 1{y_i > t}
#   B  discard the censored individuals              Ŝ_B(t) = Σ δ_i 1{y_i > t} / Σ δ_i
#   KM the Kaplan-Meier (product-limit) estimator, which uses the censored times correctly
#
# Output: output/lung.csv                       the data used (id, time, status, event)
#         output/ex_1_8_first_rows.csv          first rows shown in the text
#         figures/ex_1_8_timeline.pdf           event times and censoring times, 20 patients
#         figures/ex_1_10_km_vs_naive.pdf     the three estimates  

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
raw = dataset("survival", "lung")
y = Vector{Float64}(raw.Time)                 # observed time (days)
status = Vector{Int}(raw.Status)              # 1 = censored, 2 = dead
n = length(y)
δ = status .== 2                              # true if the event (death) was observed
# independent check against R's survival::lung (n, sum of time, sum of status, censored count)
@assert n == 228 && sum(raw.Time) == 69593 && sum(status) == 393 && count(!, δ) == 63
@assert all(y .> 0) && all(s -> s in (1, 2), status)
CSV.write(joinpath(OUTDIR, "lung.csv"),
          DataFrame(id = 1:n, time = Int.(y), status = status, event = Int.(δ)))

n_event, n_cens = count(δ), count(!, δ)
@printf("n = %d, deaths = %d, censored = %d (%.1f%%)\n", n, n_event, n_cens, 100 * n_cens / n)
@printf("observed time: min = %g, max = %g (status at max: %s)\n", minimum(y), maximum(y),
        δ[argmax(y)] ? "dead" : "censored")
@printf("mean observed time: deaths = %.1f, censored = %.1f; median: deaths = %.1f, censored = %.1f\n",
        mean(y[δ]), mean(y[.!δ]), median(y[δ]), median(y[.!δ]))

first_rows = DataFrame(id = 1:8, time = Int.(y[1:8]), status = status[1:8])
CSV.write(joinpath(OUTDIR, "ex_1_8_first_rows.csv"), first_rows)
println("first 8 rows (time, status): ", [(Int(y[i]), status[i]) for i in 1:8])

# --- the three estimates ------------------------------------------------------------
SA(t) = mean(y .> t)
SB(t) = sum(δ .& (y .> t)) / sum(δ)

# Kaplan-Meier: S(t) = Π_{t_k <= t} (1 - d_k / n_k) over the distinct event times t_k,
# where n_k = number of individuals with y_i >= t_k and d_k = number of events at t_k.
event_times = sort(unique(y[δ]))
d_k = [count(δ .& (y .== t)) for t in event_times]
n_k = [count(y .>= t) for t in event_times]
km_S = cumprod(1 .- d_k ./ n_k)
KM(t) = (i = searchsortedlast(event_times, t); i == 0 ? 1.0 : km_S[i])
@assert issorted(km_S; rev = true) && all(0 .<= km_S .<= 1)
# independent check against R's survfit(Surv(time, status == 2) ~ 1, data = lung)
@assert isapprox(KM(100), 0.86396897; atol = 1e-7)
@assert isapprox(KM(200), 0.68027286; atol = 1e-7)
@assert isapprox(KM(300), 0.53060812; atol = 1e-7)
@assert isapprox(KM(500), 0.29326919; atol = 1e-7)
@assert isapprox(KM(1000), 0.05034557; atol = 1e-7)

# when nobody is censored, KM reduces to the empirical survival function of Section 1.6
let yy = y, dd = trues(n)
    et = sort(unique(yy)); dk = [count(yy .== t) for t in et]; nk = [count(yy .>= t) for t in et]
    S0 = cumprod(1 .- dk ./ nk)
    @assert all(isapprox(S0[i], mean(yy .> et[i]); atol = 1e-12) for i in eachindex(et))
end

median_time(S, grid) = minimum(grid[[S(t) <= 0.5 for t in grid]])
grid = sort(unique(y))
@printf("median survival time (smallest t with S(t) <= 0.5): A = %g, B = %g, KM = %g\n",
        median_time(SA, grid), median_time(SB, grid), median_time(KM, grid))
for t in (100, 200, 300, 500, 800)
    @printf("t = %4d:  A = %.3f   B = %.3f   KM = %.3f\n", t, SA(t), SB(t), KM(t))
end
@printf("at the largest observed time %g: A = %.3f, KM = %.3f\n", maximum(y), SA(maximum(y)), KM(maximum(y)))

# --- figure 1: timeline of the first 20 patients -------------------------------------
m = 20
p1 = plot(; xlabel = "Time (days)", ylabel = "Patient", yticks = (1:m, string.(m:-1:1)),
          ylims = (0.4, m + 0.6), xlims = (0, 1100), legend = :outertop, legend_column = 2,
          size = (560, 460), margin = 4Plots.mm)
for i in 1:m
    yi = m - i + 1
    plot!(p1, [0, y[i]], [yi, yi]; color = BLUE, lw = 1.2, label = "")
end
scatter!(p1, y[1:m][δ[1:m]], [m - i + 1 for i in 1:m if δ[i]]; color = ORANGE, markerstrokecolor = ORANGE, markershape = :circle,
         markersize = 5, label = "death (event observed)")
scatter!(p1, y[1:m][.!δ[1:m]], [m - i + 1 for i in 1:m if !δ[i]]; markercolor = :white,
         markerstrokecolor = GREEN, markerstrokewidth = 1.5, markershape = :circle, markersize = 5, label = "censored")
savefig(p1, joinpath(FIGDIR, "ex_1_8_timeline.pdf"))

# --- figure 2: the three estimates ----------------------------------------------------
tg = range(0, maximum(y); length = 1500)
cens_t = sort(y[.!δ])
p2 = plot(tg, KM.(tg); seriestype = :steppost, color = BLUE, lw = 2.5, label = "Kaplan-Meier",
          xlabel = "Time (days)", ylabel = "S(t)", title = "Lung cancer data (n = $n)",
          ylims = (0, 1.02), size = (600, 420), margin = 4Plots.mm, legend = :topright)
plot!(p2, tg, SA.(tg); seriestype = :steppost, color = ORANGE, lw = 1.5, ls = :dash,
      label = "A: censoring times as events")
plot!(p2, tg, SB.(tg); seriestype = :steppost, color = GREEN, lw = 1.5, ls = :dot,
      label = "B: censored individuals discarded")
scatter!(p2, cens_t, KM.(cens_t); markershape = :cross, color = BLUE, markersize = 4,
         markerstrokewidth = 1, label = "censored (on Kaplan-Meier curve)")
savefig(p2, joinpath(FIGDIR, "ex_1_10_km_vs_naive.pdf"))
println("Saved figures to ", FIGDIR)
