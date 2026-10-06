#!/usr/bin/env Rscript
# Supplementary Figure 7: simulation under varying covariate dimension.
source("scripts/figure2_variant_simulations.R", chdir = FALSE)
from <- "figures/generated/variants/variant_covariates_p.pdf"
to <- "figures/generated/supplementary_figure_7.pdf"
if (!file.copy(from, to, overwrite = TRUE)) stop("Could not create ", to)
