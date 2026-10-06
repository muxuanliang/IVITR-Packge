# Helpers for PCPT decision figures. They use an authorised saved fit object
# and either its saved covariate matrix or the authorised analytic CSV.

application_d_label <- function(d) {
  value <- d / (1 - d)
  paste0("d=", sub("\\.?0+$", "", sprintf("%.2f", value)))
}

load_application_fit <- function(result_rdata, data_path = "", covariate_rdata = "") {
  if (!file.exists(result_rdata)) {
    stop("Saved application fit not found: ", result_rdata, call. = FALSE)
  }
  fit_env <- new.env(parent = emptyenv())
  load(result_rdata, envir = fit_env)
  required <- c("fit_pred", "d_seq")
  missing <- required[!vapply(required, exists, logical(1), envir = fit_env)]
  if (length(missing)) {
    stop("Saved fit is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (exists("x", envir = fit_env)) {
    x <- as.matrix(get("x", envir = fit_env))
  } else if (nzchar(covariate_rdata) && file.exists(covariate_rdata)) {
    covariate_env <- new.env(parent = emptyenv())
    load(covariate_rdata, envir = covariate_env)
    if (!exists("x", envir = covariate_env)) {
      stop("PCPT_COVARIATE_RDATA does not contain x.", call. = FALSE)
    }
    x <- as.matrix(get("x", envir = covariate_env))
  } else {
    if (!nzchar(data_path) || !file.exists(data_path)) {
      stop("The saved fit does not contain x. Set PCPT_DATA_PATH to the authorised analytic CSV.",
           call. = FALSE)
    }
    data <- utils::read.csv(data_path, stringsAsFactors = FALSE)
    data <- data[stats::complete.cases(data), , drop = FALSE]
    x <- stats::model.matrix(~ age + famhx + psa_bl + smoke + bmi, data)
    x <- x[, -1, drop = FALSE]
    x[, 1] <- scale(x[, 1])
    x[, 6] <- scale(x[, 6])
  }
  predict_index <- if (exists("predict_idx", envir = fit_env)) {
    get("predict_idx", envir = fit_env)
  } else if (exists("sampled.idx", envir = fit_env) && exists("index", envir = fit_env)) {
    which(get("sampled.idx", envir = fit_env) != get("index", envir = fit_env))
  } else {
    seq_along(get("fit_pred", envir = fit_env)[[1]])
  }
  if (length(predict_index) != length(get("fit_pred", envir = fit_env)[[1]]) ||
      max(predict_index) > nrow(x)) {
    stop("Saved predictions and covariate rows are not aligned.", call. = FALSE)
  }
  list(fit_pred = get("fit_pred", envir = fit_env), d_seq = get("d_seq", envir = fit_env),
       x = x, predict_index = predict_index)
}

save_application_decision_grid <- function(result_rdata, output_file, data_path = "",
                                           covariate_rdata = "",
                                           positive_link_treatment = 1L) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.", call. = FALSE)
  }
  input <- load_application_fit(result_rdata, data_path, covariate_rdata)
  scores <- stats::prcomp(input$x)$x[, 1:2, drop = FALSE]
  rows <- lapply(seq_along(input$d_seq), function(j) {
    d <- input$d_seq[j]
    link <- as.numeric(input$fit_pred[[j]])
    threshold <- 0.5 - d
    decision <- ifelse(link >= threshold, positive_link_treatment,
                       ifelse(link <= -threshold, -positive_link_treatment, 0L))
    data.frame(PC1 = scores[input$predict_index, 1],
               PC2 = scores[input$predict_index, 2],
               decision = factor(decision, levels = c(-1, 0, 1)),
               d_label = application_d_label(d))
  })
  plot_data <- do.call(rbind, rows)
  plot_data$d_label <- factor(plot_data$d_label, levels = application_d_label(input$d_seq))
  plot <- ggplot2::ggplot(plot_data,
                          ggplot2::aes(PC1, PC2, color = decision, shape = decision)) +
    ggplot2::geom_point(alpha = 0.3, size = 0.5, stroke = 0.3) +
    ggplot2::facet_wrap(~ d_label, ncol = 3) +
    ggplot2::scale_color_manual(name = "Decision",
      values = c("-1" = "#9EC9E2", "0" = "#E69F00", "1" = "#FF1F5B")) +
    ggplot2::scale_shape_manual(name = "Decision", values = c("-1" = 3, "0" = 2, "1" = 1)) +
    ggplot2::guides(shape = "none") + ggplot2::labs(x = "PC1", y = "PC2") +
    ggplot2::theme_gray(base_size = 12) +
    ggplot2::theme(legend.position = "bottom", legend.title = ggplot2::element_blank())
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(output_file, plot, width = 587 / 72, height = 344 / 72,
                  bg = "white", useDingbats = FALSE)
  invisible(plot_data)
}
