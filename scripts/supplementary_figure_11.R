#!/usr/bin/env Rscript
# Supplementary Figure 11: application decisions under an 80% adherence threshold.
root <- normalizePath(".")
result_rdata <- Sys.getenv("PCPT_FIGURE11_RESULT_RDATA", unset = "")
data_path <- Sys.getenv("PCPT_DATA_PATH", unset = "")
covariate_rdata <- Sys.getenv("PCPT_COVARIATE_RDATA", unset = "")
if (!nzchar(result_rdata) || !file.exists(result_rdata)) {
  stop("Set PCPT_FIGURE11_RESULT_RDATA to the authorised analysis_single_80.RData fit.",
       call. = FALSE)
}
source(file.path(root, "code", "application", "manuscript_application_figures.R"))
save_application_decision_grid(result_rdata,
  file.path(root, "figures", "generated", "supplementary_figure_11.pdf"), data_path, covariate_rdata)
