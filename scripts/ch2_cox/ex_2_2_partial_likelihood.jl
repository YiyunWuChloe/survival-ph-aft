# Section 2.2 (The Partial Likelihood) — vector version of the illustration of the five-patient example
#
# Run from the code_package folder:
#     julia --project=. scripts/ch2_cox/ex_2_2_partial_likelihood.jl
#
# The draft figure was a screenshot of a slide. This redraws it as a vector figure.
# The example in the text assumes y_1 < y_2 < y_3 < y_4 < y_5, with events for patients 1, 3, 5
# (so the risk sets are R(y_1) = {1,...,5}, R(y_3) = {3,4,5}, R(y_5) = {5}) and censored times for
# patients 2 and 4. The values used here are y = (2, 2.5, 3, 4, 5), which satisfy this ordering.
# The values y = (2.5, 4, 3, 5, 2) written in the draft do not satisfy it, so they are not used.
#
# The formulas are drawn with GR's mathtext (Computer Modern, the font of LaTeX math), because the machine has no
# LaTeX installation. The tick numbers and axis names use the normal chart font (Helvetica). The event dots are orange (solid) and
# the censored times green (open); all text is black.
#
# Output: figures/ex_2_2_partial_likelihood.pdf   (vector)

using Plots, LaTeXStrings

const FIGDIR = normpath(joinpath(@__DIR__, "..", "..", "figures"))
mkpath(FIGDIR)

y     = [2.0, 2.5, 3.0, 4.0, 5.0]               # observed times of patients 1-5
event = [true, false, true, false, true]        # patients 1, 3, 5: event; patients 2, 4: censored

# the data match the example in the text
@assert issorted(y) && allunique(y)
@assert findall(event) == [1, 3, 5] && findall(.!event) == [2, 4]
# risk sets at the three event times (individuals with y_i >= y_k)
@assert [findall(y .>= y[k]) for k in findall(event)] == [[1, 2, 3, 4, 5], [3, 4, 5], [5]]

const BAR    = RGB(0.18, 0.37, 0.52)
const GREEN  = RGB(0.18, 0.42, 0.18)
const BLACK  = RGB(0.0, 0.0, 0.0)
const ORANGE = RGB(0.88, 0.47, 0.22)

p = plot(; xlims = (0, 7.6), ylims = (0.4, 6.1), xlabel = "Time in years", ylabel = "Patient", xticks = 0:7, yticks = 1:5, framestyle = :grid, ygrid = false,
         gridalpha = 1, gridcolor = :gray88, gridlinewidth = 1, legend = false,
         size = (680, 400), margin = 3Plots.mm, fontfamily = "Helvetica")
for i in 1:5
    plot!(p, [0, y[i]], [i, i]; color = BAR, lw = 9, label = "")
end
scatter!(p, y[event], findall(event); color = ORANGE, markerstrokecolor = ORANGE,
         markershape = :circle, markersize = 9)
scatter!(p, y[.!event], findall(.!event); markercolor = :white, markerstrokecolor = GREEN,
         markerstrokewidth = 3, markershape = :circle, markersize = 9)

# contribution of each event, written next to the event time
labels = Dict(1 => L"p_1 = \frac{\exp(x_1\beta)}{\sum_{i=1}^{5}\exp(x_i\beta)}",
              3 => L"p_3 = \frac{\exp(x_3\beta)}{\sum_{i=3}^{5}\exp(x_i\beta)}",
              5 => L"p_5 = \frac{\exp(x_5\beta)}{\sum_{i=5}^{5}\exp(x_i\beta)}")
# GR's mathtext strokes are thin, so each formula is drawn twice with a tiny horizontal shift,
# which gives the heavier look of LaTeX text in a PDF
const SHIFTS = (0.0, 0.004)
for (k, lab) in labels, dx in SHIFTS
    annotate!(p, y[k] + 0.25 + dx, k, text(lab, 12, :left, BLACK))
end
for dx in SHIFTS
    annotate!(p, 0.2 + dx, 5.8, text(L"\tilde{L}(\beta) = p_1 \times p_3 \times p_5", 14, :left, BLACK))
end

savefig(p, joinpath(FIGDIR, "ex_2_2_partial_likelihood.pdf"))
println("Saved ", joinpath(FIGDIR, "ex_2_2_partial_likelihood.pdf"))
