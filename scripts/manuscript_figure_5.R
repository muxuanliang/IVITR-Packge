#!/usr/bin/env Rscript
# Main-text Figure 5: A-ERM, ERM, and BO comparison.
source("scripts/figure4_aerm_erm_bo.R", chdir = FALSE)
from <- "figures/generated/aerm_erm_bo/aerm_erm_bo_corrected.pdf"
to <- "figures/generated/manuscript_figure_5.pdf"
dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
if (!file.copy(from, to, overwrite = TRUE)) stop("Could not create ", to)
