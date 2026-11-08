# Code for "A Survey of Survival Analysis Algorithms Using Markov Chains"

Author: Yiyun Wu, The University of Queensland (MSc Statistics, Minor Project, 2026).
Supervisor: Yoni Nazarathy.

Copyright (c) 2026 Yiyun Wu. All rights reserved.

This package contains the Julia scripts that produce all figures, tables of numbers and fitted models of the report,
and the R scripts that check the fitted models. The numbering of the scripts and figures follows the chapters and
sections of the report.

## Layout

```
code_package/
  README.md
  Project.toml, Manifest.toml      Julia environment (without PhaseTypeSurvival.jl, see below)
  run_all.jl                       runs all scripts in the order of the report
  scripts/
    ch1_introduction/              Chapter 1 (Sections 1.2 to 1.10)
    ch2_cox/                       Chapter 2 (Sections 2.2 and 2.3)
    ch3_phase_type/                Chapter 3 (Sections 3.2.1 to 3.2.4)
    ch4_ph_aft/                    Chapter 4 (Section 4.4)
  figures/                         all figures (vector PDF), named as the scripts
  output/                          numbers written by the scripts (CSV)
  reference_checks_R/              R scripts that give the reference values used by the Julia scripts
```

A script `ex_<chapter>_<section>[_<subsection>]_<topic>.jl` produces the figure `ex_<chapter>_<section>[_<subsection>]_<topic>.pdf`
in `figures/`, unless stated otherwise in the table below.

## Requirements

- Julia 1.12 (earlier 1.x versions should work). Install the packages with
  `julia --project=. -e 'using Pkg; Pkg.instantiate()'` in this folder.
- The data sets `lung` (R package `survival`) and `aircondit7` (R package `boot`) are read through the Julia
  package `RDatasets.jl`. No R installation is needed to run the Julia scripts.
- The scripts in `reference_checks_R/` need R with the package `survival`.
- The four scripts `ch3_phase_type/*_simulation.jl` use the Julia package `PhaseTypeSurvival.jl`, developed by Zhihao Qiao,
  to sample from phase-type distributions and to evaluate their density and survival function. This package is not included
  in this repository. All other scripts run without it.

## Run

```
julia --project=. scripts/ch1_introduction/ex_1_8_naive_estimates.jl     # one script
julia --project=. run_all.jl                                              # all scripts
```

All random numbers use fixed seeds, so the figures are reproduced exactly.

## Scripts, figures and sections

| Script | Section of the report | Output |
|---|---|---|
| `ch1_introduction/ex_1_2_two_state.jl` | 1.2 Survival Analysis | `figures/ex_1_2_two_state.pdf` |
| `ch1_introduction/ex_1_6_jumps.jl` | 1.6 Survival Analysis Without Censoring | `figures/ex_1_6_jumps.pdf` |
| `ch1_introduction/ex_1_6_without_censoring.jl` | 1.6 | `figures/ex_1_6_without_censoring.pdf`, `output/aircondit7.csv`, `output/ex_1_6_table.csv` |
| `ch1_introduction/ex_1_7_censoring_illustration.jl` | 1.7 Survival and Censoring Times | `figures/ex_1_7_censoring_illustration.pdf` |
| `ch1_introduction/ex_1_8_lung_censoring.jl` | 1.8 Survival Analysis with Censoring and 1.10 | `figures/ex_1_8_timeline.pdf`, `figures/ex_1_10_km_vs_naive.pdf`, `output/lung.csv`, `output/ex_1_8_first_rows.csv` |
| `ch1_introduction/ex_1_8_naive_estimates.jl` | 1.8 | `figures/ex_1_8_naive_estimates.pdf` |
| `ch1_introduction/ex_1_10_km_examples.jl` | 1.10 The Kaplan-Meier Estimator | `figures/ex_1_10_km_toy.pdf`, `figures/ex_1_10_km_by_sex.pdf`, `output/ex_1_10_toy_table.csv`, `output/ex_1_10_at_risk.csv` |
| `ch2_cox/ex_2_2_partial_likelihood.jl` | 2.2 The Partial Likelihood | `figures/ex_2_2_partial_likelihood.pdf` |
| `ch2_cox/ex_2_3_cox_lung.jl` | 2.3, example with the lung data | `figures/ex_2_3_cox_lung.pdf`, `output/ex_2_3_cox_lung.csv` |
| `ch3_phase_type/ex_3_2_1_exponential_simulation.jl` | 3.2.1 Exponential Distribution | `figures/ex_3_2_1_exponential_simulation.pdf`, `output/ex_3_2_1_exponential_simulation_sample.csv` |
| `ch3_phase_type/ex_3_2_1_exponential_lung.jl` | 3.2.1, example with the lung data | `figures/ex_3_2_1_exponential_lung.pdf` |
| `ch3_phase_type/ex_3_2_2_hyperexponential_simulation.jl` | 3.2.2 Hyperexponential Distribution | `figures/ex_3_2_2_hyperexponential_simulation.pdf`, `output/ex_3_2_2_hyperexponential_simulation_sample.csv` |
| `ch3_phase_type/ex_3_2_2_hyperexponential_lung.jl` | 3.2.2, example with the lung data | `figures/ex_3_2_2_hyperexponential_lung.pdf` |
| `ch3_phase_type/ex_3_2_3_erlang_simulation.jl` | 3.2.3 Erlang Distribution | `figures/ex_3_2_3_erlang_simulation.pdf`, `output/ex_3_2_3_erlang_simulation_sample.csv` |
| `ch3_phase_type/ex_3_2_3_erlang_lung.jl` | 3.2.3, example with the lung data | `figures/ex_3_2_3_erlang_lung.pdf`, `output/ex_3_2_3_erlang_lung.csv` |
| `ch3_phase_type/ex_3_2_4_coxian_simulation.jl` | 3.2.4 Coxian Distribution | `figures/ex_3_2_4_coxian_simulation.pdf`, `output/ex_3_2_4_coxian_simulation_sample.csv` |
| `ch3_phase_type/ex_3_2_4_coxian_lung.jl` | 3.2.4, example with the lung data | `figures/ex_3_2_4_coxian_lung.pdf`, `output/ex_3_2_4_coxian_lung.csv` |
| `ch4_ph_aft/ex_4_4_ph_aft_lung.jl` | 4.4, example with the lung data | `figures/ex_4_4_ph_aft_lung.pdf`, `output/ex_4_4_ph_aft_lung.csv` |

## Checks

- Every script checks its results with `@assert` (for example, the Kaplan-Meier values and the tabulated numbers of the text).
- The scripts that fit models to the lung data compare their results with reference values from R. The R scripts that give these
  values are in `reference_checks_R/` and have the same name as the Julia script with the suffix `_check`.
  - The Cox model (Section 2.3) is compared with `survival::coxph(..., ties = "breslow")` and `basehaz()`.
  - The exponential distribution (Section 3.2.1) is compared with `survreg(..., dist = "exponential")`.
  - The hyperexponential, Erlang and Coxian fits (Sections 3.2.2 to 3.2.4) are compared with a multi-start optimisation of the same log-likelihood in R.
  - The accelerated failure time models (Section 4.4) are compared with `survreg` (exponential and Weibull baselines)
    and with a multi-start optimisation in R (phase-type baselines).

## Notes

- The sign of the regression coefficients in `ch4_ph_aft/ex_4_4_ph_aft_lung.jl` is such that a positive coefficient lengthens the survival time
  (T = exp(eta) * tau). It is the opposite sign of the option `effect = :aft` of the function `simulate_data` in `PhaseTypeSurvival.jl`.
- The phase-type distributions fitted to the lung data are fitted without covariates in Chapter 3 and with covariates in Chapter 4. The rates of a
  phase-type distribution are not unique (different rates can give the same distribution), so the scripts report the distribution (mean, median, survival
  function) and the log-likelihood.
