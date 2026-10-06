library(ggplot2)

resolve_path <- function(path){
  if (grepl("^(/|~)", path)){
    return(path.expand(path))
  }
  file.path(getwd(), path)
}

default_result_dir <- function(){
  if (dir.exists("real_data/result")){
    return("real_data/result")
  }
  "result"
}

read_data_file <- function(path){
  if (!file.exists(path)){
    stop("Data file not found: ", path, call. = FALSE)
  }
  if (requireNamespace("readr", quietly = TRUE)){
    return(readr::read_csv(path, show_col_types = FALSE))
  }
  read.csv(path, stringsAsFactors = FALSE)
}

load_covariates <- function(){
  rdata_path <- Sys.getenv("REAL_PLOT_COVARIATE_RDATA", unset = "")
  if (!nzchar(rdata_path)){
    candidates <- c("real_data/analysis_single_more.RData",
                    "analysis_single_more.RData")
    candidates <- candidates[file.exists(candidates)]
    if (length(candidates) > 0){
      rdata_path <- candidates[1]
    }
  }

  if (nzchar(rdata_path)){
    env <- new.env()
    load(resolve_path(rdata_path), envir = env)
    if (exists("x", env)){
      return(as.matrix(env$x))
    }
  }

  data_path <- Sys.getenv("REAL_PLOT_DATA_PATH",
                          unset = Sys.getenv("PCPT_DATA_PATH", unset = ""))
  if (!nzchar(data_path)){
    stop(
      "Set REAL_PLOT_COVARIATE_RDATA, REAL_PLOT_DATA_PATH, or PCPT_DATA_PATH ",
      "to an authorised input.",
      call. = FALSE
    )
  }
  data <- read_data_file(path.expand(data_path))
  data <- data[stats::complete.cases(data),]
  x <- model.matrix(~age + famhx + psa_bl + smoke + bmi, data)
  x <- x[, -1, drop = FALSE]
  x[, 1] <- scale(x[, 1])
  x[, 6] <- scale(x[, 6])
  x
}

threshold_from_file <- function(path){
  value <- sub("^.*analysis_single_([0-9]+)[.]RData$", "\\1", basename(path))
  if (identical(value, basename(path))){
    return(NA_integer_)
  }
  as.integer(value)
}

decision_from_link <- function(link, d){
  threshold <- 0.5-d
  decision <- rep(0L, length(link))
  # Preserve the original real-data plotting convention.
  decision[link >= threshold] <- 1L
  decision[link <= -threshold] <- -1L
  decision
}

decision_levels <- c("-1", "0", "1")

application_d_label <- function(d){
  value <- d / (1 - d)
  label <- sub("\\.?0+$", "", sprintf("%.2f", value))
  paste0("d=", label)
}

load_result_decisions <- function(path, pc_scores){
  env <- new.env()
  load(path, envir = env)
  required <- c("fit_pred", "d_seq")
  missing_required <- required[!vapply(required, exists, logical(1), envir = env)]
  if (length(missing_required) > 0){
    stop("Missing variables in ", path, ": ",
         paste(missing_required, collapse = ", "), call. = FALSE)
  }

  threshold <- threshold_from_file(path)
  predict_idx <- if (exists("predict_idx", env)){
    env$predict_idx
  } else {
    seq_along(env$fit_pred[[1]])
  }
  if (max(predict_idx) > nrow(pc_scores)){
    stop("predict_idx in ", basename(path), " exceeds available covariate rows.",
         call. = FALSE)
  }

  rows <- lapply(seq_along(env$d_seq), function(d_index){
    d_value <- env$d_seq[d_index]
    decision <- decision_from_link(env$fit_pred[[d_index]], d_value)
    data.frame(
      threshold = threshold,
      threshold_label = paste0(threshold, "%"),
      d = d_value,
      d_label = application_d_label(d_value),
      subject_index = predict_idx,
      PC1 = pc_scores[predict_idx, 1],
      PC2 = pc_scores[predict_idx, 2],
      decision = factor(as.character(decision), levels = decision_levels),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

decision_colors <- c("-1" = "#9EC9E2", "0" = "#E69F00", "1" = "#FF1F5B")
decision_shapes <- c("-1" = 1, "0" = 0, "1" = 2)
agreement_colors <- c(Agree = "grey35", Mixed = "purple")
decision_guides <- guides(
  color = guide_legend(
    override.aes = list(
      shape = unname(decision_shapes[decision_levels]),
      size = 2.2,
      alpha = 1
    )
  ),
  shape = "none"
)
application_width <- 587 / 72
application_height <- 344 / 72
application_point_size <- 0.5
application_point_alpha <- 0.3

application_theme <- function(base_size = 9){
  theme_gray(base_size = base_size) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.key.size = grid::unit(0.35, "lines"),
      legend.spacing.x = grid::unit(0.15, "lines"),
      legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
      plot.title = element_text(hjust = 0.5, size = base_size + 1),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size, face = "bold"),
      panel.spacing = grid::unit(0.8, "lines")
    )
}

plot_decisions <- function(data, title = NULL){
  ggplot(data, aes(PC1, PC2, color = decision, shape = decision)) +
    geom_point(alpha = application_point_alpha, size = application_point_size,
               stroke = 0.3) +
    facet_wrap(~ d_label, ncol = 3) +
    scale_color_manual(name = "Decision", values = decision_colors,
                       limits = decision_levels, breaks = decision_levels,
                       drop = FALSE) +
    scale_shape_manual(name = "Decision", values = decision_shapes,
                       limits = decision_levels, breaks = decision_levels,
                       drop = FALSE) +
    decision_guides +
    labs(title = title, x = "PC1", y = "PC2") +
    application_theme(base_size = 9)
}

result_dir <- resolve_path(Sys.getenv("REAL_PLOT_RESULT_DIR",
                                      unset = default_result_dir()))
output_dir <- resolve_path(Sys.getenv("REAL_PLOT_OUTPUT_DIR",
                                      unset = file.path(result_dir, "figures")))
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

result_files <- sort(list.files(result_dir, pattern = "^analysis_single_[0-9]+[.]RData$",
                                full.names = TRUE))
if (length(result_files) == 0){
  stop("No threshold result files found in ", result_dir, call. = FALSE)
}

x <- load_covariates()
pca <- prcomp(x)
pc_scores <- pca$x[, 1:2, drop = FALSE]

plot_data <- do.call(rbind, lapply(result_files, load_result_decisions,
                                   pc_scores = pc_scores))
plot_data$threshold_label <- factor(
  plot_data$threshold_label,
  levels = paste0(sort(unique(plot_data$threshold)), "%")
)
plot_data$d_label <- factor(plot_data$d_label,
                            levels = application_d_label(unique(plot_data$d)))

write.csv(plot_data[, c("threshold", "d", "subject_index", "PC1", "PC2",
                        "decision")],
          file.path(output_dir, "decision_plot_data.csv"),
          row.names = FALSE)

decision_counts <- as.data.frame(
  table(threshold = plot_data$threshold, d = plot_data$d,
        decision = plot_data$decision),
  stringsAsFactors = FALSE
)
names(decision_counts)[names(decision_counts) == "Freq"] <- "n"
decision_counts$total <- ave(decision_counts$n,
                             decision_counts$threshold, decision_counts$d,
                             FUN = sum)
decision_counts$prop <- decision_counts$n / decision_counts$total
write.csv(decision_counts,
          file.path(output_dir, "decision_counts_by_threshold.csv"),
          row.names = FALSE)

for (threshold in sort(unique(plot_data$threshold))){
  threshold_data <- plot_data[plot_data$threshold == threshold,]
  threshold_plot <- plot_decisions(
    threshold_data,
    title = paste0("Adherence threshold ", threshold, "%")
  )
  ggsave(
    filename = file.path(output_dir,
                         paste0("decision_pca_threshold_", threshold, ".pdf")),
    plot = threshold_plot,
    width = application_width,
    height = application_height,
    bg = "white",
    useDingbats = FALSE
  )
}

combined_plot <- ggplot(plot_data, aes(PC1, PC2, color = decision,
                                       shape = decision)) +
  geom_point(alpha = application_point_alpha, size = application_point_size,
             stroke = 0.3) +
  facet_grid(threshold_label ~ d_label) +
  scale_color_manual(name = "Decision", values = decision_colors,
                     limits = decision_levels, breaks = decision_levels,
                     drop = FALSE) +
  scale_shape_manual(name = "Decision", values = decision_shapes,
                     limits = decision_levels, breaks = decision_levels,
                     drop = FALSE) +
  decision_guides +
  labs(x = "PC1", y = "PC2") +
  application_theme(base_size = 7)
ggsave(
  filename = file.path(output_dir, "decision_pca_threshold_comparison.pdf"),
  plot = combined_plot,
  width = application_width,
  height = application_height * 1.45,
  bg = "white",
  useDingbats = FALSE
)

agreement_rows <- lapply(split(plot_data, plot_data$d), function(x_d){
  wide <- reshape(
    x_d[, c("subject_index", "threshold", "decision")],
    idvar = "subject_index",
    timevar = "threshold",
    direction = "wide"
  )
  decision_cols <- grep("^decision[.]", names(wide), value = TRUE)
  decision_mat <- as.matrix(wide[, decision_cols, drop = FALSE])
  all_same <- apply(decision_mat, 1, function(z){
    length(unique(as.character(z))) == 1
  })
  robust_decision <- ifelse(all_same, as.character(decision_mat[, 1]), "Mixed")
  data.frame(
    d = unique(x_d$d),
    d_label = unique(as.character(x_d$d_label)),
    subject_index = wide$subject_index,
    all_same = all_same,
    robust_decision = robust_decision,
    PC1 = pc_scores[wide$subject_index, 1],
    PC2 = pc_scores[wide$subject_index, 2],
    stringsAsFactors = FALSE
  )
})
agreement_data <- do.call(rbind, agreement_rows)
agreement_data$d_label <- factor(agreement_data$d_label,
                                 levels = levels(plot_data$d_label))

agreement_summary <- aggregate(all_same ~ d + d_label, agreement_data, mean)
names(agreement_summary)[names(agreement_summary) == "all_same"] <- "prop_all_thresholds_agree"
write.csv(agreement_summary,
          file.path(output_dir, "decision_agreement_by_d.csv"),
          row.names = FALSE)

pairwise_agreement <- do.call(rbind, lapply(split(plot_data, plot_data$d),
                                            function(x_d){
  thresholds <- sort(unique(x_d$threshold))
  pairs <- combn(thresholds, 2, simplify = FALSE)
  do.call(rbind, lapply(pairs, function(pair){
    x1 <- x_d[x_d$threshold == pair[1], c("subject_index", "decision")]
    x2 <- x_d[x_d$threshold == pair[2], c("subject_index", "decision")]
    merged <- merge(x1, x2, by = "subject_index", suffixes = c("_1", "_2"))
    data.frame(
      d = unique(x_d$d),
      d_label = unique(as.character(x_d$d_label)),
      threshold_pair = paste0(pair[1], "% vs ", pair[2], "%"),
      agreement = mean(merged$decision_1 == merged$decision_2),
      stringsAsFactors = FALSE
    )
  }))
}))
write.csv(pairwise_agreement,
          file.path(output_dir, "decision_pairwise_agreement_by_d.csv"),
          row.names = FALSE)

agreement_plot <- ggplot(agreement_summary,
                         aes(x = d, y = prop_all_thresholds_agree)) +
  geom_line(linewidth = 0.5, color = "black") +
  geom_point(size = 1.7, color = "black") +
  scale_x_continuous(breaks = sort(unique(agreement_summary$d))) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "d", y = "Proportion with same decision across thresholds") +
  application_theme(base_size = 9)
ggsave(
  filename = file.path(output_dir, "decision_threshold_agreement.pdf"),
  plot = agreement_plot,
  width = application_width,
  height = application_height,
  bg = "white",
  useDingbats = FALSE
)

robustness_plot <- ggplot(
  agreement_data,
  aes(PC1, PC2, color = ifelse(all_same, "Agree", "Mixed"))
) +
  geom_point(alpha = application_point_alpha, size = application_point_size,
             stroke = 0.3) +
  facet_wrap(~ d_label, ncol = 3) +
  scale_color_manual(name = "Threshold robustness",
                     values = agreement_colors) +
  labs(x = "PC1", y = "PC2") +
  application_theme(base_size = 9)
ggsave(
  filename = file.path(output_dir, "decision_threshold_robustness_pca.pdf"),
  plot = robustness_plot,
  width = application_width,
  height = application_height,
  bg = "white",
  useDingbats = FALSE
)

message("Wrote real-data plots and summaries to ", normalizePath(output_dir))
