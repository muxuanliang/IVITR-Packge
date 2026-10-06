#!/usr/bin/env Rscript

# Re-runs one main-simulation grid task (or the full 80-task grid without an
# argument). Results are written to generated/main, preserving the archived
# manuscript results in results/main.
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(".")
out <- file.path(root, "generated", "main")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
setwd(out)
source(file.path(root, "code", "simulations", "test.R"))
