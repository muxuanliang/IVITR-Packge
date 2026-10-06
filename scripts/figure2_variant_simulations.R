#!/usr/bin/env Rscript
Sys.setenv(IVITR_VARIANT_DIR = "results/variants",
           IVITR_VARIANT_FIG_OUTPUT_DIR = "figures/generated/variants")
source("code/figures/summary_variant.R", chdir = FALSE)
