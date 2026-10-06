library(ggplot2)

theme_illustration <- function() {
  theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_text(face = "bold"),
      plot.margin = margin(12, 20, 12, 12)
    )
}

save_manuscript_figure_1 <- function(path) {
  dat <- data.frame(
    patient = factor(paste("Patient", 1:5), levels = paste("Patient", 5:1)),
    lower = c(-2, 0.20, -0.25, -1.50, -0.75),
    upper = c(-1, 1.20, 1.50, 0.50, 1.00)
  )
  p <- ggplot(dat, aes(y = patient)) +
    geom_vline(xintercept = 0, colour = "darkorange", linetype = "dashed") +
    geom_segment(aes(x = lower, xend = upper, yend = patient),
                 colour = "#5B9BD5", linewidth = 0.8) +
    geom_point(aes(x = lower), colour = "#5B9BD5", size = 2) +
    geom_point(aes(x = upper), colour = "#5B9BD5", size = 2) +
    geom_text(aes(x = lower, label = "L(X)"), vjust = 1.7, parse = TRUE) +
    geom_text(aes(x = upper, label = "U(X)"), vjust = 1.7, parse = TRUE) +
    annotate("text", x = -1.45, y = 5.25, label = "D*(X) is identified as -1") +
    annotate("text", x = 0.88, y = 4.25, label = "D*(X) is identified as 1") +
    annotate("text", x = 0.45, y = 3.25, label = "D*(X) is not identifiable") +
    annotate("text", x = -0.45, y = 2.25, label = "D*(X) is not identifiable") +
    annotate("text", x = 0.20, y = 1.25, label = "D*(X) is not identifiable") +
    scale_x_continuous(limits = c(-2.5, 2), breaks = seq(-2.5, 2, 0.5)) +
    labs(x = expression("Possible values of " * Delta(X))) +
    theme_illustration()
  ggsave(path, p, width = 8, height = 5, device = grDevices::pdf,
         bg = "white", useDingbats = FALSE)
}

save_manuscript_figure_2 <- function(path) {
  dat <- data.frame(
    patient = factor(paste("Patient", 1:5), levels = paste("Patient", 5:1)),
    lower = c(-2, 0.25, -0.25, -1.50, -0.75),
    upper = c(-1, 1.20, 1.50, 0.50, 1.00)
  )
  p <- ggplot(dat, aes(y = patient)) +
    geom_vline(xintercept = 0, colour = "darkorange", linetype = "dashed") +
    geom_segment(aes(x = lower, xend = upper, yend = patient),
                 colour = "#5B9BD5", linewidth = 0.8) +
    geom_point(aes(x = lower), colour = "#5B9BD5", size = 2) +
    geom_point(aes(x = upper), colour = "#5B9BD5", size = 2) +
    annotate("segment", x = -0.25, xend = 1.5, y = 2.82, yend = 2.82,
             colour = "orange", linewidth = 0.9) +
    annotate("segment", x = -0.25, xend = 0, y = 2.62, yend = 2.62,
             colour = "#7FBF3F", linewidth = 0.9) +
    annotate("segment", x = -1.5, xend = 0, y = 1.82, yend = 1.82,
             colour = "#7FBF3F", linewidth = 0.9) +
    annotate("segment", x = 0, xend = 0.5, y = 1.82, yend = 1.82,
             colour = "orange", linewidth = 0.9) +
    annotate("segment", x = -0.75, xend = 0, y = 0.82, yend = 0.82,
             colour = "#7FBF3F", linewidth = 0.9) +
    annotate("segment", x = 0, xend = 1, y = 0.82, yend = 0.82,
             colour = "orange", linewidth = 0.9) +
    annotate("text", x = -1.25, y = 5.25,
             label = "D*(X) is identified as -1; D*_{1/3}(X) is always -1") +
    annotate("text", x = 0.78, y = 4.25,
             label = "D*(X) is identified as 1; D*_{1/3}(X) is always 1") +
    annotate("text", x = 0.62, y = 3.27,
             label = "|U(X)| is much larger than |L(X)|; D*_{1/3}(X) is 1") +
    annotate("text", x = -0.53, y = 2.27,
             label = "|L(X)| is much larger than |U(X)|; D*_{1/3}(X) is -1") +
    annotate("text", x = 0.15, y = 1.27,
             label = "|U(X)| and |L(X)| are comparable; D*_{1/3}(X) is 0") +
    scale_x_continuous(limits = c(-2.5, 2), breaks = seq(-2.5, 2, 0.5)) +
    labs(x = expression("Possible values of " * Delta(X))) +
    theme_illustration()
  ggsave(path, p, width = 8, height = 5, device = grDevices::pdf,
         bg = "white", useDingbats = FALSE)
}

save_manuscript_figure_3 <- function(path) {
  z <- seq(-1, 2, length.out = 600)
  loss <- data.frame(
    z = rep(z, 2),
    value = c(ifelse(z < 0, 1 - z / 0.5, ifelse(z < 1, 1 - z, 0)),
              pmax(1 - z, 0)),
    type = rep(c("Generalized hinge (d = 0.5)", "Hinge"), each = length(z))
  )
  p <- ggplot(loss, aes(z, value, colour = type, linetype = type)) +
    geom_vline(xintercept = c(0, 1), colour = "darkorange", linetype = "dashed") +
    geom_line(linewidth = 0.8) +
    scale_colour_manual(values = c("Generalized hinge (d = 0.5)" = "red",
                                   "Hinge" = "grey55")) +
    scale_linetype_manual(values = c("Generalized hinge (d = 0.5)" = "solid",
                                     "Hinge" = "dashed")) +
    labs(x = "z", y = "Loss(z)", colour = "Loss.Type", linetype = "Loss.Type") +
    theme_minimal(base_size = 12) + theme(legend.position = "bottom")
  ggsave(path, p, width = 7, height = 4.4, device = grDevices::pdf,
         bg = "white", useDingbats = FALSE)
}
