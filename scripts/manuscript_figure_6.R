#!/usr/bin/env Rscript
# Main-text Figure 6: PCPT application decisions across d.
root <- normalizePath(".")
result_rdata <- Sys.getenv("PCPT_FIGURE6_RESULT_RDATA", unset = "")
data_path <- Sys.getenv("PCPT_DATA_PATH", unset = "")
covariate_rdata <- Sys.getenv("PCPT_COVARIATE_RDATA", unset = "")
if (!nzchar(result_rdata) || !file.exists(result_rdata)) {
  stop(
    "Figure 6 requires an authorised saved application fit. Set PCPT_FIGURE6_RESULT_RDATA to ",
    "analysis_single_more.RData; set PCPT_DATA_PATH too if that object has no x matrix. ",
    "For the public synthetic demonstration, run scripts/figure5_synthetic_application.R instead."
  )
}
target <- file.path(root, "figures", "generated", "manuscript_figure_6.pdf")
source(file.path(root, "code", "application", "manuscript_application_figures.R"))
# The manuscript Figure 6 uses the original PCPT treatment-label convention.
save_application_decision_grid(result_rdata, target, data_path, covariate_rdata,
                               positive_link_treatment = -1L)
