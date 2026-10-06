# This file is sourced after summary_R2.R, which defines result_summary,
# make_plot(), metric_y_ranges(), and shared_legend_grob().
write_selected_main_figure <- function(settings, metrics, filename) {
  stopifnot(is.matrix(settings), ncol(settings) == 2)
  keep <- vapply(seq_len(nrow(settings)), function(i) {
    result_summary$setting1 == settings[i, 1] & result_summary$setting2 == settings[i, 2]
  }, logical(nrow(result_summary)))
  selected <- result_summary[rowSums(keep) > 0, ]
  ranges <- metric_y_ranges(selected, metrics)
  plots <- unlist(lapply(seq_len(nrow(settings)), function(i) {
    lapply(sort(unique(selected$n)), function(n_value) {
      make_plot(n_value, settings[i, 1], settings[i, 2], metrics, ranges)
    })
  }), recursive = FALSE)
  if (any(vapply(plots, is.null, logical(1)))) stop("Missing results for selected figure.")
  combined <- shared_legend_grob(plots, nrow = nrow(settings), ncol = 2, which.legend = 1)
  output <- file.path(output_dir, filename)
  ggsave(output, combined, width = 10, height = if (length(metrics) == 3) 8 else 7,
         device = grDevices::pdf, bg = "white", useDingbats = FALSE)
  message("Wrote ", normalizePath(output))
}
