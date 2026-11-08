# Section 1.2 (Survival Analysis) — Two-state survival model (Figure 1): vector diagram with circles
#
# Run from the code_package folder:
#     julia --project=. scripts/ch1_introduction/ex_1_2_two_state.jl
#
# State 0 (initial) moves to state 1 (absorbing) at hazard rate h(t). The arrow without a source
# state marks where every individual starts, as in the transition diagrams of Chapter 3.
#
# Output: figures/ex_1_2_two_state.pdf   (vector)

using Plots, LaTeXStrings

const FIGDIR = normpath(joinpath(@__DIR__, "..", "..", "figures"))
mkpath(FIGDIR)

r = 0.5                                   # radius of the circles
θ = range(0, 2π; length = 200)
circle(cx) = (cx .+ r .* cos.(θ), r .* sin.(θ))

fig = plot(; legend = false, framestyle = :none, aspect_ratio = :equal,
           xlims = (-1.5, 4.0), ylims = (-1.2, 0.9), size = (560, 220), margin = 0Plots.mm)

plot!(fig, circle(0.0); color = :black, lw = 1.5)
plot!(fig, circle(3.0); color = :black, lw = 1.5)

plot!(fig, [-1.2, -r], [0.0, 0.0]; color = :black, lw = 1.5, arrow = :closed)    # start arrow
plot!(fig, [r, 3.0 - r], [0.0, 0.0]; color = :black, lw = 1.5, arrow = :closed)  # 0 -> 1

annotate!(fig, 0.0, 0.0, text("0", 16, :black))
annotate!(fig, 3.0, 0.0, text("1", 16, :black))
annotate!(fig, 1.75, 0.22, text(L"h(t)", 14, :black))
annotate!(fig, 0.0, -0.85, text("Initial", 11, :black))
annotate!(fig, 3.0, -0.85, text("Absorbing", 11, :black))

savefig(fig, joinpath(FIGDIR, "ex_1_2_two_state.pdf"))
println("Saved ", joinpath(FIGDIR, "ex_1_2_two_state.pdf"))
