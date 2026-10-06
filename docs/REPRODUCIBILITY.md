---
output:
  pdf_document: default
  html_document: default
---
# Reproducing the manuscript results

Run every command below from the repository root. The scripts are deliberately
named by their purpose so that each manuscript figure has one visible entry
point.

## Software

The analysis was prepared with R 4.4.2. Install the R dependencies and the
local package once:

```sh
Rscript scripts/00_install_dependencies.R
```

The key non-base dependencies are `CVXR` with the `osqp` solver, `Matrix`,
`foreach`, `doParallel`, `ggplot2`, `gridExtra`, `readr`, `glmnet`, `mgcv`,
`VariableScreening`, and `personalized`. Solver installation can require
system compilers on some platforms.

## Recreate figures from archived results

The `results/` CSV files are the archived Monte Carlo outputs used by the
plotting code (not rounded summaries). The canonical figure names and commands
are listed in [`FIGURE_MAP.md`](FIGURE_MAP.md). For example, the following
commands remake main-text Figures 1–5 without rerunning simulations:

```sh
Rscript scripts/manuscript_figure_1.R
Rscript scripts/manuscript_figure_2.R
Rscript scripts/manuscript_figure_3.R
Rscript scripts/manuscript_figure_4.R
Rscript scripts/manuscript_figure_5.R
```

Outputs are written under `figures/generated/`. The scripts also write
machine-readable summary tables alongside each figure. The A-ERM/ERM/BO script
uses the corrected A-ERM/ERM outputs in `results/main/` and the archived BO
outputs in `results/bo/`. The plotting scripts retain their completeness checks
and issue a warning when an archived CSV is incomplete; in that event they use
the available replications, exactly as in the submitted figure-generation
workflow.

To run those same figure scripts and display every resulting PDF in one
compiled HTML document, use:

```sh
Rscript scripts/render_reproduction_report.R
```

The source is `reports/reproduce_figures.Rmd`; it is a fully executable,
literate record of the figure-reproduction code. Rendering requires Pandoc,
which is bundled with RStudio.

## Re-run the simulations

The full main and identifiable grids each contain 80 tasks; the additional
simulations under varying IV strength, covariate dimension, and unmeasured
confounding strength contain 60 tasks. The full settings use 200 Monte Carlo
replications and 1,000 evaluation observations per task. They are intended for
parallel/HPC execution. Each `run_*` script accepts the existing
`--task-id N` command-line option. For example:

```sh
IVITR_N_REP=2 IVITR_N_EVAL=200 Rscript scripts/02_run_main_simulations.R --task-id 1
IVITR_N_REP=2 IVITR_N_EVAL=200 Rscript scripts/03_run_wald_main_simulations.R --task-id 1
IVITR_N_REP=2 IVITR_N_EVAL=200 Rscript scripts/04_run_variant_simulations.R --task-id 1
IVITR_N_REP=2 IVITR_N_EVAL=200 Rscript scripts/06_run_identifiable_simulations.R --task-id 1
```

Use script `05` for the Wald comparator in the additional simulations. Scripts
`07` and `08` run the mild setting and its Wald comparator, respectively.
Newly generated files go to `generated/`; they do not overwrite archived
manuscript results.

The outcome and treatment models for each simulation are documented in
[`SIMULATION_MODELS.md`](SIMULATION_MODELS.md). The source scripts retain the
complete grid definitions and random seeds (`set.seed(i)`).

## Public application demonstration

The application scripts expect variables `age`, `famhx`, `psa_bl`, `smoke`,
`bmi`, `endPoint`, `trt`, and `adh`. The public synthetic dataset has exactly
those fields and is generated deterministically:

```sh
Rscript scripts/01_prepare_synthetic_application_data.R
Rscript scripts/09_run_synthetic_application.R
Rscript scripts/figure5_synthetic_application.R
```

For laptop feasibility, script `09` defaults to a 300-subject fit, three folds,
three bandwidths, and in-sample prediction. Set `IVITR_REAL_FIT_N`,
`IVITR_REAL_NFOLD`, `IVITR_REAL_BW_LEN`, and `IVITR_REAL_PREDICT_ALL` before
calling the script to change these settings for a full authorised-data run.
For a one-threshold smoke test, also set `IVITR_REAL_D_SEQ=0.5`.
If a numerical solver fails for a strict deferral setting, the public
application runner records `fallback_zero_link` in `fit_status` and continues;
this fallback is explicitly visible in the saved output.

To use an authorised copy of the original analytic data, set `PCPT_DATA_PATH`
to a CSV with the same fields before running script `09`. Do not commit that
file or the resulting participant-level outputs to a public repository.

## Authorised PCPT figure reproduction

The exact application figures have their own canonical entry points, matching
the manuscript and supplement figure numbers. They intentionally require saved
fit objects because the participant-level PCPT data are not redistributable:

```sh
PCPT_FIGURE6_RESULT_RDATA=/path/analysis_single_more.RData \
  Rscript scripts/manuscript_figure_6.R
PCPT_FIGURE9_INPUT_RDATA=/path/analysis_single_more.RData \
  Rscript scripts/supplementary_figure_9.R
PCPT_FIGURE10_RESULT_RDATA=/path/analysis_single_90.RData \
PCPT_COVARIATE_RDATA=/path/analysis_single_more.RData \
  Rscript scripts/supplementary_figure_10.R
PCPT_FIGURE11_RESULT_RDATA=/path/analysis_single_80.RData \
PCPT_COVARIATE_RDATA=/path/analysis_single_more.RData \
  Rscript scripts/supplementary_figure_11.R
```

`PCPT_DATA_PATH` can replace `PCPT_COVARIATE_RDATA` when an authorised fit
does not embed the covariate matrix. Supplementary Figure 9 additionally uses
the `personalized` implementation of the Zhao et al. (2012) rule to replace
deferred decisions. See `FIGURE_MAP.md` for the full one-to-one map.
