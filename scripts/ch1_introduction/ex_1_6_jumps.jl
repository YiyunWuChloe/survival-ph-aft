# Section 1.6 (Survival Analysis Without Censoring) — Jumps of the empirical survival function: schematic vector figure
#
# Run from the code_package folder:
#     julia --project=. scripts/ch1_introduction/ex_1_6_jumps.jl
#
# Illustrates Ŝ(t) = (1/n) Σ 1{t_i > t} for n = 5 event times 1, 2, 2, 3, 4 (arbitrary time
# units), so that t_(2) = t_(3). Ŝ is drawn as horizontal segments only. At each event time an open
# circle marks the value just before it, Ŝ(t-), and a filled dot marks the value at and just after
# it, Ŝ(t) = Ŝ(t+). The times are chosen for illustration (they are not data).
#
# Output: figures/ex_1_6_jumps.pdf   (vector)

using Plots

const FIGDIR = normpath(joinpath(@__DIR__, "..", "..", "figures"))
mkpath(FIGDIR)

const BLUE = RGB(0.18, 0.37, 0.52)

ts = [1.0, 2.0, 2.0, 3.0, 4.0]                    # t_(1) <= ... <= t_(5), with t_(2) = t_(3)
n = length(ts)
Shat(x) = count(>(x), ts) / n
@assert issorted(ts) && Shat(0.0) == 1.0 && Shat(4.0) == 0.0

ut = unique(ts)                                   # distinct event times 1, 2, 3, 4
before = [count(>=(s), ts) / n for s in ut]       # Ŝ(s-)
after  = [Shat(s) for s in ut]                    # Ŝ(s) = Ŝ(s+)
@assert before == [1.0, 0.8, 0.4, 0.2] && after == [0.8, 0.4, 0.2, 0.0]

# step function as horizontal segments only (no vertical lines);
# open circles = Ŝ(t-), filled dots = Ŝ(t) = Ŝ(t+)
edges  = [0.0; ut; 5.0]
levels = [1.0; after]
fig = plot(; xlabel = "t", ylabel = "Ŝ(t)", xlims = (0, 5), ylims = (-0.02, 1.05),
           xticks = 0:5, yticks = 0:0.2:1, size = (560, 400), margin = 4Plots.mm,
           legend = :topright, framestyle = :axes, grid = false)
for k in eachindex(levels)
    plot!(fig, [edges[k], edges[k + 1]], [levels[k], levels[k]]; color = BLUE, lw = 2, label = "")
end
scatter!(fig, ut, before; markercolor = :white, markerstrokecolor = BLUE, markerstrokewidth = 1.5,
         markersize = 5, label = "Ŝ(t−)")
scatter!(fig, ut, after; markercolor = BLUE, markerstrokecolor = BLUE, markersize = 5,
         label = "Ŝ(t) = Ŝ(t+)")

savefig(fig, joinpath(FIGDIR, "ex_1_6_jumps.pdf"))
println("Saved ", joinpath(FIGDIR, "ex_1_6_jumps.pdf"))
