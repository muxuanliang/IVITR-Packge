#!/usr/bin/env Rscript
Sys.setenv(IVITR_CONTINUOUS_RESULT_DIR = "results/continuous",
           IVITR_CONTINUOUS_FIGURE = "figures/generated/supplementary_figure_1.pdf")
source("code/figures/summary_continuous.R", chdir = FALSE)
