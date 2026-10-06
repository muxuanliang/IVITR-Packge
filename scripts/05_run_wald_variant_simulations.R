#!/usr/bin/env Rscript

# Re-runs the Wald comparator for an additional-simulation task.
root <- normalizePath(".")
out <- file.path(root, "generated", "variants")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(IVITR_VARIANT_DIR = out)
source(file.path(root, "code", "simulations", "test_qiu_variant.R"))
