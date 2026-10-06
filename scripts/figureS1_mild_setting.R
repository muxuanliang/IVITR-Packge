#!/usr/bin/env Rscript
Sys.setenv(IVITR_MILD_DIR = "results/mild",
           IVITR_MILD_FIG_OUTPUT_DIR = "figures/generated/mild")
source("code/figures/summary_mild.R", chdir = FALSE)
