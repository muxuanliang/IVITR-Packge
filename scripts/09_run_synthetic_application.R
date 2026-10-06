#!/usr/bin/env Rscript

# Run the application pipeline on the public synthetic PCPT-like data.
root <- normalizePath(".")
data_path <- file.path(root, "data", "synthetic", "pcpt_like.csv")
if (!file.exists(data_path)) source(file.path(root, "scripts", "01_prepare_synthetic_application_data.R"))
out <- file.path(root, "generated", "application", "analysis_single_85.RData")
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
Sys.setenv(PCPT_DATA_PATH = data_path, IVITR_REAL_OUTPUT = out,
           IVITR_REAL_SKIP_DIAGNOSTIC = "TRUE",
           IVITR_SOURCE_DIR = file.path(root, "code", "IVITR", "R"))
set_if_unset <- function(name, value) {
  if (!nzchar(Sys.getenv(name, unset = ""))) {
    do.call(Sys.setenv, stats::setNames(list(value), name))
  }
}
# These defaults make the public synthetic demonstration practical on a laptop.
# Set these variables before calling the script to use a larger fit.
set_if_unset("IVITR_REAL_FIT_N", "300")
set_if_unset("IVITR_REAL_NFOLD", "3")
set_if_unset("IVITR_REAL_BW_LEN", "3")
set_if_unset("IVITR_REAL_PREDICT_ALL", "FALSE")
set_if_unset("IVITR_REAL_NCORES", "1")
source(file.path(root, "code", "application", "analysis_single.R"))
