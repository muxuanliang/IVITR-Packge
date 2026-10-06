#!/usr/bin/env Rscript
# Supplementary Figure 10: application decisions under a 90% adherence threshold.
root <- normalizePath(".")
result_rdata <- Sys.getenv("PCPT_FIGURE10_RESULT_RDATA", unset = "")
data_path <- Sys.getenv("PCPT_DATA_PATH", unset = "")
covariate_rdata <- Sys.getenv("PCPT_COVARIATE_RDATA", unset = "")
if (!nzchar(result_rdata) || !file.exists(result_rdata)) {
  stop("Set PCPT_FIGURE10_RESULT_RDATA to the authorised analysis_single_90.RData fit.",
       call. = FALSE)
}
source(file.path(root, "code", "application", "manuscript_application_figures.R"))
save_application_decision_grid(result_rdata,
  file.path(root, "figures", "generated", "supplementary_figure_10.pdf"), data_path, covariate_rdata)
