#!/usr/bin/env Rscript
Sys.setenv(IVITR_FIG_RESULT_DIR = "results/main",
           IVITR_BO_RESULT_DIR = "results/bo",
           IVITR_BO_FIG_OUTPUT_DIR = "figures/generated/aerm_erm_bo")
source("code/figures/summary_bo_corrected.R", chdir = FALSE)
