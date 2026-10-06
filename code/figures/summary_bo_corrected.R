library(ggplot2)
library(grid)
library(gridExtra)

get_arg_value <- function(args, name, default=NULL){
  prefix <- paste0("--", name, "=")
  hit <- grep(paste0("^", prefix), args, value=TRUE)
  if (length(hit) > 0){
    return(sub(prefix, "", hit[1], fixed=TRUE))
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

result_file <- function(result_dir, n, eta, setting1, setting2, suffix=""){
  file.path(result_dir,
            paste0("res_", n, "_", eta, "_", eta, "_", eta, "_",
                   setting1, "_", setting2, suffix, ".csv"))
}

theme_sim <- function(){
  theme_grey(base_size=9, base_family="Helvetica") +
    theme(
      legend.position="bottom",
      legend.title=element_blank(),
      plot.background=element_rect(fill="white", colour=NA),
      legend.background=element_rect(fill="white", colour=NA),
      legend.box.background=element_rect(fill="white", colour=NA),
      panel.background=element_rect(fill="grey92", colour=NA),
      panel.grid.major=element_line(colour="white", linewidth=0.35),
      panel.grid.minor=element_line(colour="white", linewidth=0.18),
      strip.background=element_rect(fill="grey85", colour=NA),
      strip.text=element_text(size=9),
      axis.title=element_text(size=9),
      axis.text=element_text(size=8),
      plot.title=element_blank()
    )
}

args <- commandArgs(trailingOnly=TRUE)
main_dir <- resolve_path(get_arg_value(
  args, "main-dir", Sys.getenv("IVITR_FIG_RESULT_DIR", unset="simulation_R2")
))
bo_dir <- resolve_path(get_arg_value(
  args, "bo-dir", Sys.getenv("IVITR_BO_RESULT_DIR", unset="simulation")
))
output_dir <- resolve_path(get_arg_value(
  args, "output-dir", Sys.getenv("IVITR_BO_FIG_OUTPUT_DIR",
                                  unset="figures/generated/aerm_erm_bo")
))
expected_reps <- as.integer(get_arg_value(
  args, "expected-reps", Sys.getenv("IVITR_N_REP", unset="200")
))

dir.create(output_dir, showWarnings=FALSE, recursive=TRUE)

# The main simulation results contain the corrected A-ERM and ERM estimates.
# The BO column is read separately because BO is not part of the main run.
settings <- expand.grid(setting1=c(1, 2), setting2=2,
                        n=c(300, 500), eta=c(0.5, 1.5),
                        KEEP.OUT.ATTRS=FALSE)
metrics <- c("loss in value", "value on non-rejected")

read_one_setting <- function(spec){
  main_path <- result_file(main_dir, spec$n, spec$eta,
                           spec$setting1, spec$setting2, "_0.5")
  bo_path <- result_file(bo_dir, spec$n, spec$eta,
                         spec$setting1, spec$setting2, "_bo_compare")
  if (!file.exists(main_path) || !file.exists(bo_path)){
    stop("Missing comparison input for n=", spec$n, ", eta=", spec$eta,
         ", setting= (", spec$setting1, ",", spec$setting2, "): ",
         paste(c(main_path, bo_path)[!file.exists(c(main_path, bo_path))],
               collapse="; "))
  }

  main <- read.csv(main_path, stringsAsFactors=FALSE)
  main <- main[main$d == 0.5 & main$metric %in% metrics &
                 main$method %in% c("Debiased", "No-debias"),
               c("value", "metric", "method", "index")]
  main$method <- ifelse(main$method == "Debiased", "A-ERM", "ERM")

  bo <- read.csv(bo_path, stringsAsFactors=FALSE)
  bo <- bo[bo$metric %in% metrics & bo$method == "Bo",
           c("value", "metric", "method", "index")]
  bo$method <- "Bo"

  out <- rbind(main, bo)
  out$n <- spec$n
  out$eta <- spec$eta
  out$setting1 <- spec$setting1
  out$setting2 <- spec$setting2
  out
}

res <- do.call(rbind, lapply(seq_len(nrow(settings)), function(i){
  read_one_setting(settings[i,,drop=FALSE])
}))

expected_groups <- expand.grid(n=c(300, 500), eta=c(0.5, 1.5),
                               setting1=c(1, 2), setting2=2,
                               metric=metrics, method=c("A-ERM", "ERM", "Bo"),
                               KEEP.OUT.ATTRS=FALSE)
observed <- aggregate(index ~ n + eta + setting1 + setting2 + metric + method,
                      res, function(x) length(unique(x)))
names(observed)[names(observed) == "index"] <- "n_reps"
data_check <- merge(expected_groups, observed,
                    by=c("n", "eta", "setting1", "setting2", "metric", "method"),
                    all.x=TRUE, sort=FALSE)
data_check$n_reps[is.na(data_check$n_reps)] <- 0L
data_check$complete <- data_check$n_reps == expected_reps
write.csv(data_check, file.path(output_dir, "aerm_erm_bo_data_check.csv"),
          row.names=FALSE)
if (any(!data_check$complete)){
  warning(sum(!data_check$complete), " A-ERM/ERM/BO groups are incomplete; plotting available replications.")
}

metric_labels <- c("loss in value"="Loss-val",
                   "value on non-rejected"="NonDef-val")
res$metric_label <- factor(unname(metric_labels[res$metric]),
                           levels=unname(metric_labels[metrics]))
res$method <- factor(res$method, levels=c("A-ERM", "ERM", "Bo"))
res$eta_label <- factor(as.character(res$eta), levels=c("0.5", "1.5"))
write.csv(res, file.path(output_dir, "aerm_erm_bo_corrected_data.csv"),
          row.names=FALSE)

method_colors <- c("A-ERM"="red", "ERM"="orange", "Bo"="grey55")

plot_panel <- function(setting1_value, n_value){
  data_plot <- res[res$setting1 == setting1_value & res$n == n_value,]
  ggplot(data_plot, aes(x=method, y=value, colour=method)) +
    geom_boxplot(outlier.size=1, linewidth=0.35, na.rm=TRUE) +
    facet_grid(metric_label ~ eta_label, scales="free_y") +
    scale_colour_manual(values=method_colors, drop=FALSE) +
    labs(x=paste0("Sample Size ", n_value),
         y=paste0("Setting (", setting1_value, ",2)")) +
    theme_sim()
}

shared_legend_grob <- function(plots){
  legend_plot <- plots[[1]] + theme(legend.position="bottom")
  grobs <- ggplotGrob(legend_plot)$grobs
  legend_index <- which(vapply(grobs, function(x) x$name, character(1)) ==
                          "guide-box")
  if (length(legend_index) == 0){
    stop("No legend found in A-ERM/ERM/BO plot.")
  }
  legend <- grobs[[legend_index[1]]]
  legend_height <- sum(legend$height)
  panels <- gridExtra::arrangeGrob(
    grobs=lapply(plots, function(x) x + theme(legend.position="none")),
    nrow=2, ncol=2
  )
  gridExtra::arrangeGrob(
    panels, legend, ncol=1,
    heights=unit.c(unit(1, "npc") - legend_height, legend_height)
  )
}

plots <- list(
  plot_panel(1, 300), plot_panel(1, 500),
  plot_panel(2, 300), plot_panel(2, 500)
)
ggsave(
  filename=file.path(output_dir, "aerm_erm_bo_corrected.pdf"),
  plot=shared_legend_grob(plots), width=10, height=8,
  device=grDevices::pdf, bg="white", useDingbats=FALSE
)

message("Read corrected A-ERM/ERM results from ", normalizePath(main_dir))
message("Read BO results from ", normalizePath(bo_dir))
message("Wrote corrected comparison figure to ", normalizePath(output_dir))
