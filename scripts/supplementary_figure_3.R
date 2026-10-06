#!/usr/bin/env Rscript
# Supplementary Figure 3: deferral metrics for Settings (2,1) and (2,2).
source("scripts/figure1_main_simulation.R", chdir = FALSE)
source("code/figures/selected_main_figures.R", chdir = FALSE)
output_dir <- "figures/generated"
write_selected_main_figure(rbind(c(2, 1), c(2, 2)), review_metrics,
                           "supplementary_figure_3.pdf")
