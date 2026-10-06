#!/usr/bin/env Rscript

# Render the literate figure-reproduction report from the repository root.
if (!requireNamespace("rmarkdown", quietly = TRUE) ||
    !requireNamespace("knitr", quietly = TRUE)) {
  stop("Install rmarkdown and knitr with scripts/00_install_dependencies.R.")
}
if (!rmarkdown::pandoc_available()) {
  stop(
    "Pandoc 1.12.3 or later is required to render HTML. Install Pandoc or ",
    "render from RStudio (which bundles Pandoc), then rerun this script."
  )
}

repo_root <- normalizePath(".")
input <- file.path(repo_root, "reports", "reproduce_figures.Rmd")
if (!file.exists(input)) stop("Run this script from the repository root.")
rmarkdown::render(
  input = input,
  output_file = "reproduce_figures.html",
  output_dir = repo_root,
  knit_root_dir = repo_root,
  quiet = FALSE
)
