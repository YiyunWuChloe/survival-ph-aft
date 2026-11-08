# Section 1.10 (The Kaplan-Meier Estimator) — two examples that follow the definition
#
# Run from the code_package folder:
#     julia --project=. scripts/ch1_introduction/ex_1_10_km_examples.jl
#
# Example 1  A constructed data set of 7 individuals (times in months: 2, 3+, 4, 4, 6+, 7, 9, where
#            + marks a censored time). It is used to show how the product-limit estimate is
#            computed step by step. The values are chosen only to illustrate the calculation.
# Example 2  The `lung` data (R package `survival`, read through RDatasets.jl), Kaplan-Meier
#            curves for the two sexes (1 = male, 2 = female) and the number of patients at risk.
#
# The Kaplan-Meier estimator is written here directly from its definition. Its values are checked
# against R (`survival::survfit`), whose output was copied into the assertions below.
#
# Output: output/ex_1_10_toy_table.csv        the step-by-step calculation of example 1
#         output/ex_1_10_at_risk.csv          number at risk, example 2
#         figures/ex_1_10_km_toy.pdf          Kaplan-Meier curve of example 1
#         figures/ex_1_10_km_by_sex.pdf       Kaplan-Meier curves by sex, example 2  

using CSV, DataFrames, Plots, Printf

# colour palette shared by all figures of the report
const BLUE   = RGB(0.18, 0.37, 0.52)
const ORANGE = RGB(0.88, 0.47, 0.22)
const GREEN  = RGB(0.18, 0.42, 0.18)
const LBLUE  = RGB(0.74, 0.84, 0.92)

using RDatasets

const FIGDIR = normpath(joinpath(@__DIR__, "..", "..", "figures"))
const OUTDIR = normpath(joinpath(@__DIR__, "..", "..", "output"))
mkpath(FIGDIR); mkpath(OUTDIR)

# --- Kaplan-Meier from the definition ------------------------------------------------
# t_1 < ... < t_K are the distinct times with at least one event, d_k the number of events at t_k,
# n_k the number of individuals at risk just before t_k, and S(t) = prod_{t_k <= t} (n_k - d_k) / n_k.
function kaplan_meier(y, δ)
    tk = sort(unique(y[δ]))
    dk = [count((y .== t) .& δ) for t in tk]
    nk = [count(y .>= t) for t in tk]
    S = cumprod((nk .- dk) ./ nk)
    return (; tk, dk, nk, S)
end
km_at(km, t) = (i = searchsortedlast(km.tk, t); i == 0 ? 1.0 : km.S[i])
median_time(km) = (i = findfirst(km.S .<= 0.5); i === nothing ? NaN : km.tk[i])

# --- example 1: constructed data ---------------------------------------------------------
y1 = [2.0, 3, 4, 4, 6, 7, 9]
δ1 = [true, false, true, true, false, true, true]
km1 = kaplan_meier(y1, δ1)
@assert km1.tk == [2, 4, 7, 9] && km1.dk == [1, 2, 1, 1] && km1.nk == [7, 5, 2, 1]
@assert km1.S ≈ [6/7, 6/7 * 3/5, 6/7 * 3/5 * 1/2, 0.0]
@assert median_time(km1) == 7
# values from R: survfit(Surv(time, status) ~ 1) on the same seven observations
@assert isapprox.(km1.S, [0.857, 0.514, 0.257, 0.0]; atol = 5e-4) |> all
@printf("example 1: t_k = %s, n_k = %s, d_k = %s\n", km1.tk, km1.nk, km1.dk)
@printf("           S(t_k) = %s, median = %g\n", round.(km1.S, digits = 3), median_time(km1))

CSV.write(joinpath(OUTDIR, "ex_1_10_toy_table.csv"),
          DataFrame(k = 1:4, t = km1.tk, n_at_risk = km1.nk, n_events = km1.dk,
                    factor = (km1.nk .- km1.dk) ./ km1.nk, S = km1.S))

cens1 = y1[.!δ1]
p1 = plot(; xlabel = "Time (months)", ylabel = "Estimated survival probability",
          xlims = (0, 9.5), ylims = (0, 1.02), xticks = 0:9, yticks = 0:0.25:1,
          legend = false, size = (520, 380), margin = 4Plots.mm)
plot!(p1, [0; km1.tk], [1.0; km1.S]; seriestype = :steppost, color = BLUE, lw = 2.5)
plot!(p1, [0, median_time(km1), median_time(km1)], [0.5, 0.5, 0]; color = ORANGE, ls = :dot, lw = 1.5)
scatter!(p1, cens1, km_at.(Ref(km1), cens1); markershape = :cross, markercolor = BLUE,
         markerstrokecolor = BLUE, markerstrokewidth = 2, markersize = 8)
savefig(p1, joinpath(FIGDIR, "ex_1_10_km_toy.pdf"))

# --- example 2: lung data by sex -----------------------------------------------------------
raw = dataset("survival", "lung")
y = Vector{Float64}(raw.Time); δ = raw.Status .== 2; sex = Vector{Int}(raw.Sex)
@assert length(y) == 228 && sum(δ) == 165
@assert count(sex .== 1) == 138 && count(sex .== 2) == 90
km_m = kaplan_meier(y[sex .== 1], δ[sex .== 1])
km_f = kaplan_meier(y[sex .== 2], δ[sex .== 2])
@assert sum(km_m.dk) == 112 && sum(km_f.dk) == 53
@assert median_time(km_m) == 270 && median_time(km_f) == 426           # R: median 270 and 426

grid = [0, 200, 400, 600, 800]
at_risk(sx, t) = count((sex .== sx) .& (y .>= t))
risk_m = at_risk.(1, grid); risk_f = at_risk.(2, grid)
# R: summary(survfit(Surv(time, status == 2) ~ sex, data = lung), times = grid, extend = TRUE)
@assert risk_m == [138, 78, 31, 13, 6] && risk_f == [90, 66, 26, 11, 2]
@assert isapprox(km_at.(Ref(km_m), grid), [1.0, 0.6073, 0.2977, 0.1451, 0.0670]; atol = 5e-5)
@assert isapprox(km_at.(Ref(km_f), grid), [1.0, 0.7946, 0.5089, 0.3433, 0.0832]; atol = 5e-5)
@printf("example 2: n = %d men, %d women; deaths %d and %d; medians %g and %g days\n",
        count(sex .== 1), count(sex .== 2), sum(km_m.dk), sum(km_f.dk), median_time(km_m), median_time(km_f))
@printf("           at risk men   %s\n           at risk women %s\n", risk_m, risk_f)
@printf("           S men   %s\n           S women %s\n", round.(km_at.(Ref(km_m), grid), digits = 3),
        round.(km_at.(Ref(km_f), grid), digits = 3))
CSV.write(joinpath(OUTDIR, "ex_1_10_at_risk.csv"),
          DataFrame(time = grid, men = risk_m, women = risk_f,
                    S_men = km_at.(Ref(km_m), grid), S_women = km_at.(Ref(km_f), grid)))

p2 = plot(; xlabel = "Time (days)", ylabel = "Estimated survival probability",
          xlims = (0, 1050), ylims = (0, 1.02), legend = :topright, size = (600, 420),
          margin = 4Plots.mm)
for (km, yy, dd, ls, lab, col) in ((km_m, y[sex .== 1], δ[sex .== 1], :solid, "Men", BLUE),
                                   (km_f, y[sex .== 2], δ[sex .== 2], :dash, "Women", ORANGE))
    tx = [0; km.tk; maximum(yy)]                      # extend the last step to the largest observed time
    plot!(p2, tx, [1.0; km.S; km.S[end]]; seriestype = :steppost, color = col, lw = 1.8, ls = ls, label = lab)
    c = yy[.!dd]
    scatter!(p2, c, km_at.(Ref(km), c); markershape = :cross, markercolor = col, markerstrokecolor = col,
             markersize = 4, label = "")
end
savefig(p2, joinpath(FIGDIR, "ex_1_10_km_by_sex.pdf"))
println("Saved figures to ", FIGDIR)
