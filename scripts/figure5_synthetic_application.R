#!/usr/bin/env Rscript

root <- normalizePath(".")
result_dir <- file.path(root, "generated", "application")
if (!dir.exists(result_dir) || !length(list.files(result_dir, pattern = "[.]RData$"))) {
  stop("Run scripts/09_run_synthetic_application.R first.")
}
out <- file.path(root, "figures", "generated", "synthetic_application")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(REAL_PLOT_RESULT_DIR = result_dir, REAL_PLOT_OUTPUT_DIR = out,
           REAL_PLOT_DATA_PATH = file.path(root, "data", "synthetic", "pcpt_like.csv"))
source(file.path(root, "code", "application", "plot.R"))
