#!/usr/bin/env Rscript
# Supplementary Figure 4: overall performance for Settings (1,1) and (2,1).
source("scripts/figure1_main_simulation.R", chdir = FALSE)
source("code/figures/selected_main_figures.R", chdir = FALSE)
output_dir <- "figures/generated"
write_selected_main_figure(rbind(c(1, 1), c(2, 1)), performance_metrics,
                           "supplementary_figure_4.pdf")
