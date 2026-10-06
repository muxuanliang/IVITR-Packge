#!/usr/bin/env Rscript

root <- normalizePath(".")
out <- file.path(root, "generated", "mild")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(IVITR_MILD_DIR = out)
source(file.path(root, "code", "simulations", "test_qiu_mild.R"))
