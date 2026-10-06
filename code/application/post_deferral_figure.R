# Reproduce Supplementary Figure 9: replace deferred decisions with the
# Zhao et al. (2012) outcome-weighted-learning rule, as in the supplement.

save_post_deferral_grid <- function(input_rdata, output_file) {
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("personalized", quietly = TRUE) ||
      !requireNamespace("IVITR", quietly = TRUE)) {
    stop("Supplementary Figure 9 requires ggplot2, personalized, and IVITR.",
         call. = FALSE)
  }
  if (!file.exists(input_rdata)) stop("Input not found: ", input_rdata, call. = FALSE)
  input <- new.env(parent = globalenv())
  load(input_rdata, envir = input)
  required <- c("data", "x", "sampled.idx", "index", "fit", "d_seq", "z_tfd")
  missing <- required[!vapply(required, exists, logical(1), envir = input)]
  if (length(missing)) stop("Input is missing: ", paste(missing, collapse = ", "), call. = FALSE)

  data <- get("data", envir = input)
  x <- as.matrix(get("x", envir = input))
  sampled_idx <- get("sampled.idx", envir = input)
  split_index <- get("index", envir = input)
  train <- sampled_idx == split_index
  test <- !train
  y <- ifelse(data$endPoint == "Healthy", 1, -1)
  fit_rej <- personalized::fit.subgroup(
    x = x[train, , drop = FALSE], y = y[train],
    trt = get("z_tfd", envir = input)[train],
    propensity.func = function(x, trt) 0.5,
    loss = "owl_hinge_flip_loss", method = "weighting"
  )

  pca <- stats::prcomp(x)
  scores <- pca$x[test, 1:2, drop = FALSE]
  d_seq <- get("d_seq", envir = input)
  fit <- get("fit", envir = input)
  rows <- lapply(seq_along(d_seq), function(j) {
    d <- d_seq[j]
    link <- as.numeric(IVITR::predict.iv_itr(fit[[j]], x[test, , drop = FALSE], type = "link"))
    threshold <- 0.5 - d
    decision <- ifelse(link >= threshold, -1L,
                       ifelse(link <= -threshold, 1L, 0L))
    deferred <- which(decision == 0L)
    if (length(deferred)) {
      decision[deferred] <- as.integer(
        stats::predict(fit_rej, newx = x[test, , drop = FALSE][deferred, , drop = FALSE],
                       type = "trt.group")
      )
    }
    data.frame(PC1 = scores[, 1], PC2 = scores[, 2],
               decision = factor(decision, levels = c(-1, 1)),
               d_label = application_d_label(d))
  })
  plot_data <- do.call(rbind, rows)
  plot_data$d_label <- factor(plot_data$d_label, levels = application_d_label(d_seq))
  plot <- ggplot2::ggplot(plot_data,
                          ggplot2::aes(PC1, PC2, color = decision, shape = decision)) +
    ggplot2::geom_point(alpha = 0.3, size = 0.5, stroke = 0.3) +
    ggplot2::facet_wrap(~ d_label, ncol = 3) +
    ggplot2::scale_color_manual(name = "Decision", values = c("-1" = "#9EC9E2", "1" = "#FF1F5B")) +
    ggplot2::scale_shape_manual(name = "Decision", values = c("-1" = 3, "1" = 1)) +
    ggplot2::guides(shape = "none") + ggplot2::labs(x = "PC1", y = "PC2") +
    ggplot2::theme_gray(base_size = 12) +
    ggplot2::theme(legend.position = "bottom", legend.title = ggplot2::element_blank())
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(output_file, plot, width = 587 / 72, height = 344 / 72,
                  bg = "white", useDingbats = FALSE)
}
