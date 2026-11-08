# Runs all scripts of the code package in the order of the report and checks that each of them finishes.
#
# Run from the code_package folder:
#     julia --project=. run_all.jl
#
# The four scripts "*_simulation.jl" of Chapter 3 use the Julia package PhaseTypeSurvival.jl, which is not part of this
# package. They are skipped, with a message, if that package is not available in the active environment (see README.md).

const SCRIPTS = [
    "ch1_introduction/ex_1_2_two_state.jl",
    "ch1_introduction/ex_1_6_jumps.jl",
    "ch1_introduction/ex_1_6_without_censoring.jl",
    "ch1_introduction/ex_1_7_censoring_illustration.jl",
    "ch1_introduction/ex_1_8_lung_censoring.jl",
    "ch1_introduction/ex_1_8_naive_estimates.jl",
    "ch1_introduction/ex_1_10_km_examples.jl",
    "ch2_cox/ex_2_2_partial_likelihood.jl",
    "ch2_cox/ex_2_3_cox_lung.jl",
    "ch3_phase_type/ex_3_2_1_exponential_simulation.jl",
    "ch3_phase_type/ex_3_2_1_exponential_lung.jl",
    "ch3_phase_type/ex_3_2_2_hyperexponential_simulation.jl",
    "ch3_phase_type/ex_3_2_2_hyperexponential_lung.jl",
    "ch3_phase_type/ex_3_2_3_erlang_simulation.jl",
    "ch3_phase_type/ex_3_2_3_erlang_lung.jl",
    "ch3_phase_type/ex_3_2_4_coxian_simulation.jl",
    "ch3_phase_type/ex_3_2_4_coxian_lung.jl",
    "ch4_ph_aft/ex_4_4_ph_aft_lung.jl",
]

has_phasetype = Base.find_package("PhaseTypeSurvival") !== nothing
failed = String[]
for s in SCRIPTS
    if occursin("_simulation", s) && !has_phasetype
        println("SKIPPED (needs PhaseTypeSurvival.jl): ", s)
        continue
    end
    println("---- ", s)
    # every script runs in its own Julia process, so that the scripts do not share variables
    ok = success(pipeline(`$(Base.julia_cmd()) --project=$(@__DIR__) $(joinpath(@__DIR__, "scripts", s))`; stdout = devnull, stderr = stderr))
    ok || push!(failed, s)
end
println(isempty(failed) ? "Done: all scripts finished." : "These scripts failed: " * join(failed, ", "))
