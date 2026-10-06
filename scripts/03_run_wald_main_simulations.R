#!/usr/bin/env Rscript

# Re-runs the Wald comparator. Pass --task-id N to select an array task.
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(".")
out <- file.path(root, "generated", "main")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
setwd(out)
source(file.path(root, "code", "simulations", "test_qiu.R"))
