library(ggplot2)

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
  args, "result-dir", Sys.getenv("IVITR_MILD_DIR", unset = "simulation_R2_mild")
))
output_dir <- resolve_path(get_arg_value(
  args, "output-dir", Sys.getenv("IVITR_MILD_FIG_OUTPUT_DIR",
                                unset = result_dir)
))
expected_reps <- as.integer(get_arg_value(
  args, "expected-reps", Sys.getenv("IVITR_N_REP", unset = "200")
))

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

result_files <- sort(list.files(result_dir, pattern = "^res_mild_.*[.]csv$",
                                full.names = TRUE))
if (length(result_files) == 0){
  stop("No mild result CSV files found in ", result_dir, ".")
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
    eta = x$eta[1],
    delta = x$delta[1],
    setting1 = x$setting1[1],
    setting2 = x$setting2[1],
    d = x$d[1],
    p_covariates = x$p_covariates[1],
    za_coef = x$za_coef[1],
    x2_coef = x$x2_coef[1],
    u_a_multiplier = x$u_a_multiplier[1]
  )
}))
file_completeness <- file_completeness[order(file_completeness$n,
                                             file_completeness$d,
                                             file_completeness$source_type),]
write.csv(file_completeness,
          file.path(output_dir, "mild_result_file_completeness.csv"),
          row.names = FALSE)

if (any(!file_completeness$complete)){
  warning(sum(!file_completeness$complete),
          " mild result files are incomplete; plots use available non-NA values.")
}

res <- raw_res[
  (raw_res$source_is_qiu & raw_res$method == "Qiu") |
    (!raw_res$source_is_qiu & raw_res$method %in% c("Debiased", "No-debias")),
]
if (!any(res$method == "Qiu")){
  warning("No Qiu rows found; Ad-hoc curves will be missing from the plots.")
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

split_vars <- c("n", "eta", "delta", "setting1", "setting2", "method", "d",
                "metric", "p_covariates", "za_coef", "x2_coef",
                "u_a_multiplier")
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
  "Qiu" = "Ad-hoc"
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
method_levels <- intersect(c("Augmented ERM", "ERM", "Ad-hoc"),
                           unique(result_summary$method_label))
result_summary$method_label <- factor(result_summary$method_label,
                                      levels = method_levels)
write.csv(result_summary, file.path(output_dir, "mild_result_summary.csv"),
          row.names = FALSE)

method_colors <- c(
  "Augmented ERM" = "#F8766D",
  "ERM" = "#00BA38",
  "Ad-hoc" = "#619CFF"
)
method_linetypes <- c(
  "Augmented ERM" = "solid",
  "ERM" = "dashed",
  "Ad-hoc" = "longdash"
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

metric_y_ranges <- function(data){
  data <- data[!is.na(data$mean),]
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

  x_value <- min(data$d, na.rm = TRUE)
  blanks <- do.call(rbind, lapply(seq_len(nrow(y_ranges)), function(i){
    expand.grid(
      metric_label = y_ranges$metric_label[i],
      n_label = unique(data$n_label),
      d = x_value,
      mean = c(y_ranges$ymin[i], y_ranges$ymax[i]),
      KEEP.OUT.ATTRS = FALSE
    )
  }))
  blanks$metric_label <- factor(blanks$metric_label,
                                levels = levels(result_summary$metric_label))
  blanks$n_label <- factor(blanks$n_label, levels = levels(data$n_label))
  blanks
}

plot_metric_group <- function(data, metrics){
  data <- data[data$metric %in% metrics,]
  data <- data[order(data$metric_label, data$n, data$method_label, data$d),]
  data$n_label <- factor(paste0("n = ", data$n),
                         levels = paste0("n = ", sort(unique(data$n))))

  y_ranges <- metric_y_ranges(data)
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
    facet_grid(metric_label ~ n_label, scales = "free_y") +
    scale_x_continuous(breaks = sort(unique(data$d))) +
    scale_color_manual(values = method_colors[method_levels], drop = FALSE) +
    scale_linetype_manual(values = method_linetypes[method_levels],
                          drop = FALSE) +
    labs(x = "d", y = "Metric value") +
    theme_sim()
}

performance_metrics <- c("loss in value", "value on non-rejected",
                         "porp rejected")
review_metrics <- c("prop deferred with bounds crossing zero",
                    "prop wrong treatment among non-deferred")

performance_plot <- plot_metric_group(result_summary, performance_metrics)
review_plot <- plot_metric_group(result_summary, review_metrics)

ggsave(
  filename = file.path(output_dir, "mild_performance_metrics.pdf"),
  plot = performance_plot,
  width = 7, height = 6,
  device = grDevices::pdf,
  bg = "white",
  useDingbats = FALSE
)
ggsave(
  filename = file.path(output_dir, "mild_review_metrics.pdf"),
  plot = review_plot,
  width = 7, height = 4.5,
  device = grDevices::pdf,
  bg = "white",
  useDingbats = FALSE
)

message("Read mild results from ", normalizePath(result_dir))
message("Wrote mild plots and summaries to ", normalizePath(output_dir))
