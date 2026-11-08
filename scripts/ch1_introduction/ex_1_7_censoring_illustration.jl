# Section 1.7 (Survival and Censoring Times) — vector version of the illustration of censored data
#
# Run from the code_package folder:
#     julia --project=. scripts/ch1_introduction/ex_1_7_censoring_illustration.jl
#
# This redraws the figure that was a raster screenshot in the draft ("Illustration of censored
# survival data"), with the same values as described in the text of Section 1.7:
#   patient 1: event at 2.5 years        patient 2: dropped out, censored at 4 years
#   patient 3: still alive at the end of the 5-year study, censored at 5 years
#   patient 4: event at 3 years
#
# Output: figures/ex_1_7_censoring_illustration.pdf   (vector)

using Plots

const FIGDIR = normpath(joinpath(@__DIR__, "..", "..", "figures"))
mkpath(FIGDIR)

const END_OF_STUDY = 5.0
time_obs = [2.5, 4.0, 5.0, 3.0]            # observed time of patients 1-4
event    = [true, false, false, true]      # true = event observed, false = censored

# the data are as described in the text
@assert all(time_obs .<= END_OF_STUDY)
@assert findall(.!event) == [2, 3]

const BAR    = RGB(0.18, 0.37, 0.52)       # dark blue
const GREEN  = RGB(0.18, 0.42, 0.18)
const ORANGE = RGB(0.88, 0.47, 0.22)

p = plot(; xlabel = "Time in years", ylabel = "Patient", xlims = (0, 6.3), ylims = (0.4, 4.6),
         xticks = 0:6, yticks = 1:4, framestyle = :grid, ygrid = false,
         gridalpha = 1, gridcolor = :gray88, gridlinewidth = 1, legend = false, size = (600, 300),
         margin = 3Plots.mm, fontfamily = "Helvetica")
plot!(p, [END_OF_STUDY, END_OF_STUDY], [0.4, 4.6]; color = GREEN, lw = 1.5, ls = :dot)
for i in 1:4
    plot!(p, [0, time_obs[i]], [i, i]; color = BAR, lw = 9, label = "")
end
scatter!(p, time_obs[event], findall(event); color = ORANGE, markerstrokecolor = ORANGE,
         markershape = :circle, markersize = 9)
scatter!(p, time_obs[.!event], findall(.!event); markercolor = :white, markerstrokecolor = GREEN,
         markerstrokewidth = 3, markershape = :circle, markersize = 9)
savefig(p, joinpath(FIGDIR, "ex_1_7_censoring_illustration.pdf"))
println("Saved ", joinpath(FIGDIR, "ex_1_7_censoring_illustration.pdf"))
