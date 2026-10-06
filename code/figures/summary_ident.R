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

summary_one <- function(value){
  value <- value[!is.na(value)]
  n_obs <- length(value)
  mean_value <- if (n_obs == 0) NA_real_ else mean(value)
  sd_value <- if (n_obs <= 1) NA_real_ else sd(value)
  se <- if (n_obs <= 1) NA_real_ else sd_value / sqrt(n_obs)
  data.frame(
    mean=mean_value,
    sd=sd_value,
    se=se,
    n_obs=n_obs,
    lower=mean_value - 1.96 * se,
    upper=mean_value + 1.96 * se
  )
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
result_dir <- resolve_path(get_arg_value(
  args, "result-dir", Sys.getenv("IVITR_IDENT_DIR", unset="simulation_R2_ident")
))
output_dir <- resolve_path(get_arg_value(
  args, "output-dir", Sys.getenv("IVITR_IDENT_FIG_OUTPUT_DIR", unset=result_dir)
))
expected_reps <- as.integer(get_arg_value(
  args, "expected-reps", Sys.getenv("IVITR_N_REP", unset="200")
))

dir.create(output_dir, showWarnings=FALSE, recursive=TRUE)
result_files <- sort(list.files(result_dir, pattern="^res_ident_.*[.]csv$",
                                full.names=TRUE))
if (length(result_files) == 0){
  stop("No identifiable result CSV files found in ", result_dir, ".")
}

raw_res <- do.call(rbind, lapply(result_files, function(path){
  x <- read.csv(path, stringsAsFactors=FALSE)
  x$source_file <- basename(path)
  x
}))

file_completeness <- do.call(rbind, lapply(split(raw_res, raw_res$source_file),
                                           function(x){
  data.frame(
    source_file=x$source_file[1], rows=nrow(x),
    reps=length(unique(x$index)), expected_rows=expected_reps * 30,
    expected_reps=expected_reps,
    complete=nrow(x) == expected_reps * 30 &&
      length(unique(x$index)) == expected_reps,
    n=x$n[1], eta=x$eta[1], setting1=x$setting1[1],
    setting2=x$setting2[1], d=x$d[1]
  )
}))
file_completeness <- file_completeness[order(file_completeness$setting1,
                                             file_completeness$setting2,
                                             file_completeness$n,
                                             file_completeness$eta,
                                             file_completeness$d),]
write.csv(file_completeness,
          file.path(output_dir, "ident_result_file_completeness.csv"),
          row.names=FALSE)
if (any(!file_completeness$complete)){
  warning(sum(!file_completeness$complete),
          " identifiable result files are incomplete; plots use available values.")
}

metrics <- c("worst-case loss in value",
             "prop wrong treatment among non-deferred")
methods <- c("Augmented ERM", "ERM", "Ad-hoc")
res <- raw_res[raw_res$metric %in% metrics & raw_res$method %in% methods,]
if (nrow(res) == 0){
  stop("None of the requested identifiable metrics were found.")
}

# The Ad-hoc rule never defers.  Use its d=0.5 estimates as one reference line
# over every d value, rather than mixing in separately simulated Ad-hoc fits.
proposed <- res[res$method != "Ad-hoc",]
adhoc <- res[res$method == "Ad-hoc" & abs(res$d - 0.5) < 1e-12,]
if (nrow(adhoc) == 0){
  stop("No Ad-hoc rows at d=0.5 were found.")
}
target_d <- unique(proposed[c("n", "eta", "setting1", "setting2", "d")])
adhoc_reference <- merge(
  adhoc[setdiff(names(adhoc), "d")], target_d,
  by=c("n", "eta", "setting1", "setting2"), all=FALSE, sort=FALSE
)
res_plot <- rbind(proposed, adhoc_reference[, names(proposed)])

split_vars <- c("n", "eta", "setting1", "setting2", "method", "d", "metric")
result_summary <- do.call(rbind, lapply(split(res_plot, res_plot[split_vars],
                                              drop=TRUE), function(x){
  cbind(x[1, split_vars], summary_one(x$value))
}))
rownames(result_summary) <- NULL

metric_labels <- c(
  "worst-case loss in value"="Loss-val",
  "prop wrong treatment among non-deferred"="PropIncorrectTrt"
)
result_summary$metric_label <- factor(unname(metric_labels[result_summary$metric]),
                                      levels=unname(metric_labels[metrics]))
method_labels <- c("Augmented ERM"="Augmented ERM", "ERM"="ERM", "Ad-hoc"="Wald")
result_summary$method_label <- factor(unname(method_labels[result_summary$method]),
                                      levels=unname(method_labels[methods]))
write.csv(result_summary, file.path(output_dir, "ident_result_summary.csv"),
          row.names=FALSE)

method_colors <- c("Augmented ERM"="#F8766D", "ERM"="#00BA38", "Wald"="#619CFF")
method_linetypes <- c("Augmented ERM"="solid", "ERM"="dashed", "Wald"="longdash")

metric_y_ranges <- function(data){
  do.call(rbind, lapply(split(data[!is.na(data$mean),],
                            data$metric_label[!is.na(data$mean)], drop=TRUE),
                         function(x){
    y_range <- range(x$mean)
    if (diff(y_range) == 0){
      y_range <- y_range + c(-0.02, 0.02)
    }
    data.frame(metric_label=x$metric_label[1], ymin=y_range[1], ymax=y_range[2])
  }))
}

y_range_blanks <- function(data, y_ranges){
  x_value <- min(data$d, na.rm=TRUE)
  blanks <- do.call(rbind, lapply(seq_len(nrow(y_ranges)), function(i){
    expand.grid(
      metric_label=y_ranges$metric_label[i],
      eta_label=levels(data$eta_label),
      d=x_value,
      mean=c(y_ranges$ymin[i], y_ranges$ymax[i]),
      KEEP.OUT.ATTRS=FALSE
    )
  }))
  blanks$metric_label <- factor(blanks$metric_label,
                                levels=levels(data$metric_label))
  blanks$eta_label <- factor(blanks$eta_label, levels=levels(data$eta_label))
  blanks
}

plot_n <- function(data, n_value, y_ranges){
  data <- data[data$n == n_value,]
  data$eta_label <- factor(paste0("eta = ", data$eta),
                           levels=paste0("eta = ", sort(unique(data$eta))))
  blanks <- y_range_blanks(data, y_ranges)
  ggplot(data, aes(x=d, y=mean, colour=method_label,
                   linetype=method_label, group=method_label)) +
    geom_point(size=1.5, na.rm=TRUE) +
    geom_line(linewidth=0.5, na.rm=TRUE) +
    geom_blank(data=blanks, aes(x=d, y=mean), inherit.aes=FALSE) +
    facet_grid(metric_label ~ eta_label, scales="free_y") +
    scale_x_continuous(breaks=sort(unique(data$d))) +
    scale_colour_manual(values=method_colors, drop=FALSE) +
    scale_linetype_manual(values=method_linetypes, drop=FALSE) +
    labs(x=paste0("d   (n = ", n_value, ")"), y="Metric value") +
    theme_sim()
}

shared_legend_grob <- function(plots){
  legend_plot <- plots[[1]] + theme(legend.position="bottom")
  grobs <- ggplotGrob(legend_plot)$grobs
  legend <- grobs[[which(vapply(grobs, function(x) x$name, character(1)) ==
                         "guide-box")[1]]]
  legend_height <- sum(legend$height)
  panels <- do.call(gridExtra::arrangeGrob,
                    c(lapply(plots, function(x) x + theme(legend.position="none")),
                      ncol=1))
  gridExtra::arrangeGrob(
    panels, legend, ncol=1,
    heights=unit.c(unit(1, "npc") - legend_height, legend_height)
  )
}

for (setting1_value in sort(unique(result_summary$setting1))){
  for (setting2_value in sort(unique(result_summary$setting2))){
    data_plot <- result_summary[
      result_summary$setting1 == setting1_value &
        result_summary$setting2 == setting2_value,
    ]
    n_values <- sort(unique(data_plot$n))
    y_ranges <- metric_y_ranges(data_plot)
    plots <- lapply(n_values, function(n_value) plot_n(data_plot, n_value, y_ranges))
    if (length(plots) == 0){
      next
    }
    ggsave(
      filename=file.path(output_dir,
                         paste0("ident_selected_metrics_setting_", setting1_value,
                                "_", setting2_value, "_combined.pdf")),
      plot=shared_legend_grob(plots), width=7.5, height=9,
      device=grDevices::pdf, bg="white", useDingbats=FALSE
    )
  }
}

message("Read identifiable results from ", normalizePath(result_dir))
message("Wrote plots and summaries to ", normalizePath(output_dir))
