#!/usr/bin/env Rscript
# Main-text Figure 4: overall performance for Settings (1,2) and (2,2).
source("scripts/figure1_main_simulation.R", chdir = FALSE)
source("code/figures/selected_main_figures.R", chdir = FALSE)
output_dir <- "figures/generated"
write_selected_main_figure(
  settings = rbind(c(1, 2), c(2, 2)),
  metrics = performance_metrics,
  filename = "manuscript_figure_4.pdf"
)
