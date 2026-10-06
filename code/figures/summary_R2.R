library(ggplot2)
library(grid)
library(gridExtra)

get_arg_value <- function(args, name, default = NULL){
  prefix <- paste0("--", name, "=")
  hit <- grep(paste0("^", prefix), args, value = TRUE)
  if (length(hit) > 0){
    return(sub(prefix, "", hit[1], fixed = TRUE))
  }
  idx <- which(args == paste0("--", name))
  if (length(idx) > 0 && idx[1] < length(args)){
    return(args[idx[1] + 1])
  }
  default
}

resolve_path <- function(path){
  if (grepl("^(/|~)", path)){
    return(path.expand(path))
  }
  file.path(getwd(), path)
}

args <- commandArgs(trailingOnly = TRUE)
result_dir <- resolve_path(get_arg_value(
  args, "result-dir", Sys.getenv("IVITR_FIG_RESULT_DIR", unset = "simulation_R2")
))
output_dir <- resolve_path(get_arg_value(
  args, "output-dir", Sys.getenv("IVITR_FIG_OUTPUT_DIR", unset = result_dir)
))
expected_reps <- as.integer(get_arg_value(
  args, "expected-reps", Sys.getenv("IVITR_N_REP", unset = "200")
))

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

result_files <- sort(list.files(result_dir, pattern = "^res_.*[.]csv$",
                                full.names = TRUE))
if (length(result_files) == 0){
  stop("No result CSV files found in ", result_dir, ".")
}

read_result <- function(path){
  x <- read.csv(path, stringsAsFactors = FALSE)
  x$source_file <- basename(path)
  x$source_is_qiu <- grepl("_qiu[.]csv$", x$source_file)
  x
}

raw_res <- do.call(rbind, lapply(result_files, read_result))

apply_qiu_zero_defer_correction <- function(raw_res){
  target_rows <- raw_res$source_is_qiu &
    raw_res$method == "Debiased" &
    raw_res$metric == "porp rejected" &
    !is.na(raw_res$value) &
    raw_res$value <= .Machine$double.eps
  if (!any(target_rows)){
    return(raw_res)
  }

  zero_defer_keys <- paste(raw_res$source_file[target_rows],
                           raw_res$index[target_rows], sep = "\r")
  row_keys <- paste(raw_res$source_file, raw_res$index, sep = "\r")
  qiu_zero_rows <- raw_res$source_is_qiu &
    raw_res$method == "Qiu" &
    row_keys %in% zero_defer_keys

  raw_res$value[qiu_zero_rows & raw_res$metric == "porp rejected"] <- 0
  raw_res$value[
    qiu_zero_rows &
      raw_res$metric == "prop deferred with bounds crossing zero"
  ] <- NA_real_

  message("Corrected ", length(unique(zero_defer_keys)),
          " Qiu replications with zero target deferral.")
  raw_res
}

raw_res <- apply_qiu_zero_defer_correction(raw_res)

file_completeness <- do.call(rbind, lapply(split(raw_res, raw_res$source_file),
                                           function(x){
  data.frame(
    source_file = x$source_file[1],
    source_type = if (x$source_is_qiu[1]) "Qiu" else "Main",
    rows = nrow(x),
    reps = length(unique(x$index)),
    expected_rows = expected_reps * 10,
    expected_reps = expected_reps,
    complete = nrow(x) == expected_reps * 10 &&
      length(unique(x$index)) == expected_reps,
    n = x$n[1],
    delta = x$delta[1],
    setting1 = x$setting1[1],
    setting2 = x$setting2[1],
    d = x$d[1]
  )
}))
file_completeness <- file_completeness[order(file_completeness$n,
                                             file_completeness$delta,
                                             file_completeness$setting1,
                                             file_completeness$setting2,
                                             file_completeness$d,
                                             file_completeness$source_type),]
write.csv(file_completeness,
          file.path(output_dir, "result_file_completeness.csv"),
          row.names = FALSE)

if (any(!file_completeness$complete)){
  warning(sum(!file_completeness$complete),
          " result files are incomplete; plots use available non-NA values.")
}

res <- raw_res[
  (raw_res$source_is_qiu & raw_res$method == "Qiu") |
    (!raw_res$source_is_qiu & raw_res$method %in% c("Debiased", "No-debias")),
]
if (!any(res$method == "Qiu")){
  warning("No Qiu rows found; Wald curves will be missing from the plots.")
}

summary_one <- function(value){
  value <- value[!is.na(value)]
  n_obs <- length(value)
  mean_value <- if (n_obs == 0) NA_real_ else mean(value)
  sd_value <- if (n_obs <= 1) NA_real_ else sd(value)
  se <- if (n_obs <= 1) NA_real_ else sd_value / sqrt(n_obs)
  data.frame(
    mean = mean_value,
    sd = sd_value,
    se = se,
    n_obs = n_obs,
    lower = mean_value - 1.96 * se,
    upper = mean_value + 1.96 * se
  )
}

split_vars <- c("n", "delta", "setting1", "setting2", "method", "d", "metric")
result_summary <- do.call(rbind, lapply(split(res, res[split_vars], drop = TRUE),
                                        function(x){
  cbind(x[1, split_vars], summary_one(x$value))
}))
rownames(result_summary) <- NULL

metric_labels <- c(
  "loss in value" = "Loss-val",
  "value on non-rejected" = "NonDef-val",
  "porp rejected" = "Def-prop",
  "prop deferred with bounds crossing zero" = "PropCorrectDef",
  "prop wrong treatment among non-deferred" = "PropIncorrectTrt"
)
method_labels <- c(
  "Debiased" = "Augmented ERM",
  "No-debias" = "ERM",
  "Qiu" = "Wald"
)

result_summary$metric_label <- unname(metric_labels[result_summary$metric])
result_summary$method_label <- unname(method_labels[result_summary$method])
result_summary$metric_label[is.na(result_summary$metric_label)] <-
  result_summary$metric[is.na(result_summary$metric_label)]
result_summary$method_label[is.na(result_summary$method_label)] <-
  result_summary$method[is.na(result_summary$method_label)]
result_summary$metric_label <- factor(
  result_summary$metric_label,
  levels = c("Loss-val", "NonDef-val", "Def-prop",
             "PropCorrectDef", "PropIncorrectTrt")
)
method_levels <- intersect(c("Augmented ERM", "ERM", "Wald"),
                           unique(result_summary$method_label))
result_summary$method_label <- factor(result_summary$method_label,
                                      levels = method_levels)
write.csv(result_summary, file.path(output_dir, "result_summary.csv"),
          row.names = FALSE)

method_colors <- c(
  "Augmented ERM" = "#F8766D",
  "ERM" = "#00BA38",
  "Wald" = "#619CFF"
)
method_linetypes <- c(
  "Augmented ERM" = "solid",
  "ERM" = "dashed",
  "Wald" = "longdash"
)

theme_sim <- function(){
  theme_grey(base_size = 9, base_family = "Helvetica") +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      plot.background = element_rect(fill = "white", colour = NA),
      legend.background = element_rect(fill = "white", colour = NA),
      legend.box.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "grey92", colour = NA),
      panel.grid.major = element_line(colour = "white", linewidth = 0.35),
      panel.grid.minor = element_line(colour = "white", linewidth = 0.18),
      strip.background = element_rect(fill = "grey85", colour = NA),
      strip.text = element_text(size = 9),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 8),
      plot.title = element_blank()
    )
}

metric_y_ranges <- function(data, metrics){
  data <- data[data$metric %in% metrics & !is.na(data$mean),]
  if (nrow(data) == 0){
    return(NULL)
  }

  ranges <- do.call(rbind, lapply(split(data, data$metric_label, drop = TRUE),
                                  function(x){
    y_range <- range(x$mean, na.rm = TRUE)
    if (!all(is.finite(y_range))){
      return(NULL)
    }
    if (diff(y_range) == 0){
      pad <- max(abs(y_range[1]) * 0.04, 0.02)
      y_range <- y_range + c(-pad, pad)
    }
    data.frame(metric_label = x$metric_label[1],
               ymin = y_range[1], ymax = y_range[2])
  }))
  rownames(ranges) <- NULL
  ranges
}

y_range_blanks <- function(data, y_ranges){
  if (is.null(y_ranges) || nrow(y_ranges) == 0){
    return(NULL)
  }

  deltas <- sort(unique(data$delta))
  x_value <- min(data$d, na.rm = TRUE)
  blanks <- do.call(rbind, lapply(seq_len(nrow(y_ranges)), function(i){
    expand.grid(
      metric_label = y_ranges$metric_label[i],
      delta = deltas,
      d = x_value,
      mean = c(y_ranges$ymin[i], y_ranges$ymax[i]),
      KEEP.OUT.ATTRS = FALSE
    )
  }))
  blanks$metric_label <- factor(blanks$metric_label,
                                levels = levels(result_summary$metric_label))
  blanks
}

plot_metric_group <- function(data, metrics, x_label, y_label,
                              y_ranges = NULL){
  data <- data[data$metric %in% metrics,]
  blanks <- y_range_blanks(data, y_ranges)
  p <- ggplot(data, aes(x = d, y = mean, color = method_label,
                        linetype = method_label, group = method_label)) +
    geom_point(size = 1.5, na.rm = TRUE) +
    geom_line(linewidth = 0.5, na.rm = TRUE)
  if (!is.null(blanks)){
    p <- p + geom_blank(data = blanks, aes(x = d, y = mean),
                        inherit.aes = FALSE)
  }
  p +
    facet_grid(metric_label ~ delta, scales = "free_y") +
    scale_x_continuous(breaks = sort(unique(data$d))) +
    scale_color_manual(values = method_colors[method_levels], drop = FALSE) +
    scale_linetype_manual(values = method_linetypes[method_levels],
                          drop = FALSE) +
    labs(x = x_label, y = y_label) +
    theme_sim()
}

make_plot <- function(n_value, setting1_value, setting2_value, metrics,
                      y_ranges = NULL){
  data_plot <- result_summary[
    result_summary$n == n_value &
      result_summary$setting1 == setting1_value &
      result_summary$setting2 == setting2_value,
  ]
  if (nrow(data_plot) == 0){
    return(NULL)
  }

  plot_metric_group(
    data_plot, metrics,
    x_label = paste("Sample Size", n_value),
    y_label = paste0("Setting (", setting1_value, ",", setting2_value, ")"),
    y_ranges = y_ranges
  )
}

shared_legend_grob <- function(plots, nrow = 1, ncol = length(plots),
                               position = c("bottom", "right"),
                               which.legend = 1){
  position <- match.arg(position)
  legend_plot <- plots[[which.legend]] +
    theme(legend.title = element_blank(), legend.position = position) +
    guides(color = guide_legend(nrow = 1, byrow = TRUE),
           linetype = guide_legend(nrow = 1, byrow = TRUE))
  grobs <- ggplotGrob(legend_plot)$grobs
  legend_index <- which(vapply(grobs, function(x) x$name, character(1)) ==
                          "guide-box")
  if (length(legend_index) == 0){
    stop("No legend found in the selected plot.")
  }
  legend <- grobs[[legend_index[1]]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)

  plot_grobs <- lapply(plots, function(x) x + theme(legend.position = "none"))
  arranged <- do.call(gridExtra::arrangeGrob,
                      c(plot_grobs, nrow = nrow, ncol = ncol))

  switch(position,
         "bottom" = gridExtra::arrangeGrob(
           arranged, legend, ncol = 1,
           heights = unit.c(unit(1, "npc") - lheight, lheight)
         ),
         "right" = gridExtra::arrangeGrob(
           arranged, legend, ncol = 2,
           widths = unit.c(unit(1, "npc") - lwidth, lwidth)
         ))
}

save_combined <- function(filename, plots, width, height){
  combined <- shared_legend_grob(plots, nrow = 2, ncol = 2,
                                 which.legend = 1)
  ggsave(
    filename = file.path(output_dir, filename),
    plot = combined,
    width = width, height = height,
    device = grDevices::pdf,
    bg = "white",
    useDingbats = FALSE
  )
}

performance_metrics <- c("loss in value", "value on non-rejected",
                         "porp rejected")
review_metrics <- c("prop deferred with bounds crossing zero",
                    "prop wrong treatment among non-deferred")

for (n_value in sort(unique(result_summary$n))){
  for (setting1_value in sort(unique(result_summary$setting1))){
    for (setting2_value in sort(unique(result_summary$setting2))){
      data_plot <- result_summary[
        result_summary$n == n_value &
          result_summary$setting1 == setting1_value &
          result_summary$setting2 == setting2_value,
      ]
      if (nrow(data_plot) == 0){
        next
      }

      performance_plot <- plot_metric_group(
        data_plot, performance_metrics,
        x_label = paste("Sample Size", n_value),
        y_label = paste0("Setting (", setting1_value, ",", setting2_value, ")"),
        y_ranges = metric_y_ranges(data_plot, performance_metrics)
      )
      review_plot <- plot_metric_group(
        data_plot, review_metrics,
        x_label = paste("Sample Size", n_value),
        y_label = paste0("Setting (", setting1_value, ",", setting2_value, ")"),
        y_ranges = metric_y_ranges(data_plot, review_metrics)
      )

      ggsave(
        filename = file.path(output_dir,
                             paste0("performance_metrics_Setting_",
                                    setting1_value, "_", setting2_value,
                                    "_n_", n_value, ".pdf")),
        plot = performance_plot,
        width = 6, height = 6,
        device = grDevices::pdf,
        bg = "white",
        useDingbats = FALSE
      )
      ggsave(
        filename = file.path(output_dir,
                             paste0("review_metrics_Setting_",
                                    setting1_value, "_", setting2_value,
                                    "_n_", n_value, ".pdf")),
        plot = review_plot,
        width = 6, height = 4,
        device = grDevices::pdf,
        bg = "white",
        useDingbats = FALSE
      )
    }
  }
}

for (setting2_value in sort(unique(result_summary$setting2))){
  combined_data <- result_summary[result_summary$setting2 == setting2_value,]
  perf_y_ranges <- metric_y_ranges(combined_data, performance_metrics)
  review_y_ranges <- metric_y_ranges(combined_data, review_metrics)

  perf_plots <- list(
    make_plot(300, 1, setting2_value, performance_metrics, perf_y_ranges),
    make_plot(500, 1, setting2_value, performance_metrics, perf_y_ranges),
    make_plot(300, 2, setting2_value, performance_metrics, perf_y_ranges),
    make_plot(500, 2, setting2_value, performance_metrics, perf_y_ranges)
  )
  review_plots <- list(
    make_plot(300, 1, setting2_value, review_metrics, review_y_ranges),
    make_plot(500, 1, setting2_value, review_metrics, review_y_ranges),
    make_plot(300, 2, setting2_value, review_metrics, review_y_ranges),
    make_plot(500, 2, setting2_value, review_metrics, review_y_ranges)
  )

  if (!any(vapply(perf_plots, is.null, logical(1)))){
    save_combined(paste0("performance_metrics_setting2_",
                         setting2_value, "_combined.pdf"),
                  perf_plots, width = 10, height = 8)
  }
  if (!any(vapply(review_plots, is.null, logical(1)))){
    save_combined(paste0("review_metrics_setting2_",
                         setting2_value, "_combined.pdf"),
                  review_plots, width = 10, height = 7)
  }
}

message("Read results from ", normalizePath(result_dir))
message("Wrote plots and summaries to ", normalizePath(output_dir))
