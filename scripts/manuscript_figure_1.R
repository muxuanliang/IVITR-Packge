#!/usr/bin/env Rscript
source("code/figures/manuscript_illustrations.R", chdir = FALSE)
dir.create("figures/generated", recursive = TRUE, showWarnings = FALSE)
save_manuscript_figure_1("figures/generated/manuscript_figure_1.pdf")
