#!/usr/bin/env Rscript

# Install the packages needed to reproduce the analyses, then install the
# local IVITR package.  Run this once from the repository root:
# Rscript scripts/00_install_dependencies.R

repo_root <- normalizePath(".")
if (!dir.exists(file.path(repo_root, "code", "IVITR"))) {
  stop("Run this script from the repository root.")
}
required <- c("CVXR", "osqp", "Matrix", "foreach", "doParallel", "ggplot2",
              "gridExtra", "readr", "glmnet", "mgcv", "VariableScreening",
              "knitr", "rmarkdown", "personalized")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  install.packages(missing, repos = "https://cloud.r-project.org")
}

pkg_dir <- file.path(repo_root, "code", "IVITR")
if (!dir.exists(pkg_dir)) stop("Run this script from the repository root.")
install.packages(pkg_dir, repos = NULL, type = "source")
message("Dependencies and local IVITR package installed.")
