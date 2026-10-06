#!/usr/bin/env Rscript
# Supplementary Figure 8: identifiable / no-U-to-A setting.
source("scripts/figure3_identifiable_setting.R", chdir = FALSE)
from <- "figures/generated/identifiable/ident_selected_metrics_setting_2_2_combined.pdf"
to <- "figures/generated/supplementary_figure_8.pdf"
if (!file.copy(from, to, overwrite = TRUE)) stop("Could not create ", to)
