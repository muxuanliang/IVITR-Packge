#!/usr/bin/env Rscript

# Re-runs the mild-setting A-ERM/ERM simulation. Use the matching Wald script
# below to complete the comparison.
root <- normalizePath(".")
out <- file.path(root, "generated", "mild")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(IVITR_MILD_DIR = out)
source(file.path(root, "code", "simulations", "test_mild.R"))
