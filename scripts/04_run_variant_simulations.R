#!/usr/bin/env Rscript

# Re-runs one additional-simulation task (or all 60 default tasks).
root <- normalizePath(".")
out <- file.path(root, "generated", "variants")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(IVITR_VARIANT_DIR = out)
source(file.path(root, "code", "simulations", "test_variant.R"))
