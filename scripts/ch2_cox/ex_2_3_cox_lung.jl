# Section 2.3 (example after The Breslow Estimator) — the Cox proportional hazards model for the lung data
#
# Run from the code_package folder:
#     julia --project=. scripts/ch2_cox/ex_2_3_cox_lung.jl
#
# Data: `lung` from the R package `survival` (read through RDatasets.jl), n = 228, 165 deaths, 63 censored.
# Model: h(t | x) = h0(t) exp(beta_1 (age - 62) + beta_2 female), female = 1 for women, 0 for men.
#        Age is centred at 62 (about the mean age, 62.4), so that the baseline hazard h0 is the hazard of a man aged 62.
#        A centring constant only changes h0 (see the discussion of the intercept in Section 2.1), not beta.
# Ties: the lung data have tied event times (165 deaths at 139 distinct times), so the partial likelihood uses the
#       Breslow approximation  L(beta) = prod_k exp(sum_{i in D_k} x_i beta) / (sum_{i in R(t_k)} exp(x_i beta))^{d_k}.
#       The Breslow estimator of the cumulative baseline hazard increments by d_k / sum_{i in R(t_k)} exp(x_i beta_hat).
# Everything is computed from scratch with the formulas of Chapter 2 and checked against R's
# survival::coxph(..., ties = "breslow") and basehaz() (reference values hard-coded below).
#
# Output: output/ex_2_3_cox_lung.csv        coefficient table
#         figures/ex_2_3_cox_lung.pdf       figure (vector)

using Plots, LinearAlgebra, Printf, DataFrames, CSV, RDatasets

const BLUE   = RGB(0.18, 0.37, 0.52)
const ORANGE = RGB(0.88, 0.47, 0.22)
const FIGDIR = normpath(joinpath(@__DIR__, "..", "..", "figures"))
const OUTDIR = normpath(joinpath(@__DIR__, "..", "..", "output"))
mkpath(FIGDIR); mkpath(OUTDIR)

raw = dataset("survival", "lung")
y = Vector{Float64}(raw.Time)
δ = raw.Status .== 2
age = Vector{Float64}(raw.Age)
female = Float64.(raw.Sex .== 2)
n = length(y)
@assert n == 228 && count(δ) == 165 && count(female .== 1) == 90
X = hcat(age .- 62, female)                        # n x 2 covariate matrix

# ---- partial likelihood (Breslow approximation), score and information ---------------------------------
tk = sort(unique(y[δ]))
K = length(tk)
dk = [count(δ .& (y .== t)) for t in tk]
@assert K == 139 && sum(dk) == 165
risk = [findall(y .>= t) for t in tk]               # risk sets R(t_k)
xsum_events = [vec(sum(X[δ .& (y .== t), :]; dims = 1)) for t in tk]   # sum of x_i over the events at t_k

function loglik_score_info(β)
    ℓ = 0.0; U = zeros(2); I = zeros(2, 2)
    w = exp.(X * β)
    for k in 1:K
        r = risk[k]
        s0 = sum(w[r]); s1 = vec(sum(w[r] .* X[r, :]; dims = 1)); s2 = X[r, :]' * (w[r] .* X[r, :])
        ℓ += dot(xsum_events[k], β) - dk[k] * log(s0)
        U .+= xsum_events[k] .- dk[k] .* s1 ./ s0
        I .+= dk[k] .* (s2 ./ s0 .- (s1 ./ s0) * (s1 ./ s0)')
    end
    return ℓ, U, I
end

β = zeros(2)
ℓ0 = loglik_score_info(β)[1]
for it in 1:25                                      # Newton-Raphson
    ℓ, U, I = loglik_score_info(β)
    step = I \ U
    β .+= step
    norm(step) < 1e-12 && break
end
ℓhat, U, I = loglik_score_info(β)
@assert norm(U) < 1e-8
V = inv(I)
se = sqrt.(diag(V))
z = β ./ se
HR = exp.(β)
lo = exp.(β .- 1.96 .* se); hi = exp.(β .+ 1.96 .* se)
@printf("loglik at 0: %.4f, at the maximum: %.4f\n", ℓ0, ℓhat)
for (j, nm) in enumerate(("age (per year)", "female"))
    @printf("%-15s beta = %9.6f  se = %8.6f  z = %7.4f  HR = %.4f  [%.4f, %.4f]\n", nm, β[j], se[j], z[j], HR[j], lo[j], hi[j])
end

# checks against R (survival::coxph, ties = "breslow")
@assert isapprox(β[1], 0.01701289; atol = 1e-6) && isapprox(β[2], -0.51256479; atol = 1e-6)
@assert isapprox(se[1], 0.00922195; atol = 1e-6) && isapprox(se[2], 0.16746206; atol = 1e-6)
@assert isapprox(ℓ0, -750.122; atol = 1e-3) && isapprox(ℓhat, -743.0797; atol = 1e-3)
@assert isapprox(lo[2], 0.431372; atol = 1e-5) && isapprox(hi[2], 0.831649; atol = 1e-5)

CSV.write(joinpath(OUTDIR, "ex_2_3_cox_lung.csv"),
          DataFrame(term = ["age (per year)", "female"], beta = β, se = se, z = z, HR = HR, HR_lower = lo, HR_upper = hi))

# ---- Breslow estimator of the cumulative baseline hazard (x = 0: a man aged 62) -------------------------
w = exp.(X * β)
h_k = [dk[k] / sum(w[risk[k]]) for k in 1:K]
H0 = cumsum(h_k)
H0_at(t) = (i = searchsortedlast(tk, t); i == 0 ? 0.0 : H0[i])
@assert isapprox(H0_at(100), 0.1718173; atol = 1e-6) && isapprox(H0_at(300), 0.7548918; atol = 1e-6)
@assert isapprox(H0_at(1000), 3.463283; atol = 1e-5)
S_cox(t, x) = exp(-H0_at(t) * exp(dot(x, β)))
@assert isapprox(S_cox(300, [0.0, 0.0]), 0.4700615; atol = 1e-6) && isapprox(S_cox(300, [0.0, 1.0]), 0.63626; atol = 1e-5)
@assert isapprox(S_cox(500, [0.0, 1.0]), 0.4144298; atol = 1e-5)
median_cox(x) = tk[findfirst(t -> S_cox(t, x) <= 0.5, tk)]
@printf("model-based median survival (age 62): men %g, women %g days\n", median_cox([0.0, 0.0]), median_cox([0.0, 1.0]))
@assert median_cox([0.0, 0.0]) == 285 && median_cox([0.0, 1.0]) == 394

# ---- Kaplan-Meier by sex, for comparison (as in Section 1.10) -------------------------------------------
function km(yy, dd)
    t = sort(unique(yy[dd])); S = 1.0; out = Float64[]
    for u in t
        S *= 1 - count(dd .& (yy .== u)) / count(yy .>= u); push!(out, S)
    end
    return t, out
end

# ---- figure ---------------------------------------------------------------------------------------------
tmax = maximum(y)
pa = plot([0; tk; tmax], [0.0; H0; H0[end]]; seriestype = :steppost, color = BLUE, lw = 2, label = "",
          xlabel = "Time (days)", ylabel = "Cumulative baseline hazard", title = "(a) Breslow estimate",
          xlims = (0, 1050), legend = false)
pb = plot(; xlabel = "Time (days)", ylabel = "Estimated survival probability", title = "(b) Cox model and Kaplan-Meier",
          xlims = (0, 1050), ylims = (0, 1.02), legend = :topright)
tg = sort(unique([0.0; tk; tmax]))
for (x, col, ls, lab) in (([0.0, 0.0], BLUE, :solid, "Men (Cox model)"), ([0.0, 1.0], ORANGE, :dash, "Women (Cox model)"))
    plot!(pb, tg, [S_cox(t, x) for t in tg]; seriestype = :steppost, color = col, lw = 2.2, ls = ls, label = lab)
end
for (sx, col, lab) in ((1, BLUE, "Men (Kaplan-Meier)"), (2, ORANGE, "Women (Kaplan-Meier)"))
    idx = raw.Sex .== sx
    tt, SS = km(y[idx], δ[idx])
    plot!(pb, [0; tt; maximum(y[idx])], [1.0; SS; SS[end]]; seriestype = :steppost, color = col, lw = 1, alpha = 0.5, label = lab)
end
fig = plot(pa, pb; layout = (1, 2), size = (1100, 420), margin = 4Plots.mm)
savefig(fig, joinpath(FIGDIR, "ex_2_3_cox_lung.pdf"))
println("Saved ", joinpath(FIGDIR, "ex_2_3_cox_lung.pdf"))
