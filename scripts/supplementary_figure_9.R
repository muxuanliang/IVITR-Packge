#!/usr/bin/env Rscript
# Supplementary Figure 9: post-deferral PCPT decisions across d.
root <- normalizePath(".")
input_rdata <- Sys.getenv("PCPT_FIGURE9_INPUT_RDATA", unset = "")
if (!nzchar(input_rdata) || !file.exists(input_rdata)) {
  stop("Set PCPT_FIGURE9_INPUT_RDATA to the authorised analysis_single_more.RData input.",
       call. = FALSE)
}
source(file.path(root, "code", "application", "manuscript_application_figures.R"))
source(file.path(root, "code", "application", "post_deferral_figure.R"))
save_post_deferral_grid(
  input_rdata,
  file.path(root, "figures", "generated", "supplementary_figure_9.pdf")
)
