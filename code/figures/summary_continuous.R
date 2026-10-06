library(ggplot2)
library(grid)
library(gridExtra)

continuous_dir <- Sys.getenv("IVITR_CONTINUOUS_RESULT_DIR", unset = "results/continuous")
continuous_output <- Sys.getenv("IVITR_CONTINUOUS_FIGURE", unset = "figures/generated/supplementary_figure_1.pdf")
files <- sort(list.files(continuous_dir, pattern = "_cont(_qiu)?[.]csv$", full.names = TRUE))
if (!length(files)) stop("No continuous-outcome result files found in ", continuous_dir)
res <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
res <- res[res$setting1 == 1 & res$setting2 == 2 &
             res$metric %in% c("loss in value", "value on non-rejected", "porp rejected"), ]
res$metric <- factor(res$metric,
                     levels = c("loss in value", "value on non-rejected", "porp rejected"),
                     labels = c("Loss-val", "NonDef-val", "Def-prop"))
res$method <- factor(res$method,
                     levels = c("Debiased", "No-debias", "Qiu"),
                     labels = c("Augmented ERM", "ERM", "Ad-hoc"))
summary_res <- aggregate(value ~ n + delta + d + metric + method, res, mean, na.rm = TRUE)
summary_res$eta_label <- factor(paste0("eta = ", summary_res$delta),
                                levels = c("eta = 0.5", "eta = 1.5"))

make_panel <- function(n_value) {
  data <- summary_res[summary_res$n == n_value, ]
  ggplot(data, aes(d, value, colour = method, linetype = method, group = method)) +
    geom_line(linewidth = 0.5) + geom_point(size = 1.4) +
    facet_grid(metric ~ eta_label, scales = "free_y") +
    scale_x_continuous(breaks = sort(unique(data$d))) +
    scale_colour_manual(values = c("Augmented ERM" = "#F8766D", "ERM" = "#00BA38", "Ad-hoc" = "#619CFF")) +
    scale_linetype_manual(values = c("Augmented ERM" = "solid", "ERM" = "dashed", "Ad-hoc" = "longdash")) +
    labs(x = paste0("d (n = ", n_value, ")"), y = "Metric value") +
    theme_grey(base_size = 9) + theme(legend.position = "bottom", legend.title = element_blank())
}

plots <- lapply(c(300, 500), make_panel)
legend_grob <- {
  grobs <- ggplotGrob(plots[[1]])$grobs
  grobs[[which(vapply(grobs, function(x) x$name, character(1)) == "guide-box")[1]]]
}
combined <- arrangeGrob(
  arrangeGrob(plots[[1]] + theme(legend.position = "none"),
              plots[[2]] + theme(legend.position = "none"), ncol = 2),
  legend_grob, ncol = 1, heights = unit.c(unit(1, "npc") - sum(legend_grob$height), sum(legend_grob$height))
)
dir.create(dirname(continuous_output), recursive = TRUE, showWarnings = FALSE)
ggsave(continuous_output, combined, width = 10, height = 5, device = grDevices::pdf,
       bg = "white", useDingbats = FALSE)
message("Wrote ", normalizePath(continuous_output))
