#!/usr/bin/env Rscript
Sys.setenv(IVITR_IDENT_DIR = "results/identifiable",
           IVITR_IDENT_FIG_OUTPUT_DIR = "figures/generated/identifiable")
source("code/figures/summary_ident.R", chdir = FALSE)
