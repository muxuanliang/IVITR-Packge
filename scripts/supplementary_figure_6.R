#!/usr/bin/env Rscript
# Supplementary Figure 6: simulation under varying IV strength.
source("scripts/figure2_variant_simulations.R", chdir = FALSE)
from <- "figures/generated/variants/variant_iv_strength.pdf"
to <- "figures/generated/supplementary_figure_6.pdf"
if (!file.copy(from, to, overwrite = TRUE)) stop("Could not create ", to)
