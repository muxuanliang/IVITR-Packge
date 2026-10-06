library(readr)
library(IVITR)
library(doParallel)
library(foreach)
find_ivitr_source_dir <- function(){
  override <- Sys.getenv("IVITR_SOURCE_DIR", unset = "")
  candidates <- c(
    override,
    file.path(getwd(), "IV-ITR", "IVITR", "R"),
    file.path(getwd(), "IVITR", "R"),
    file.path(dirname(getwd()), "IV-ITR", "IVITR", "R")
  )
  candidates <- candidates[nzchar(candidates)]
  candidates <- candidates[dir.exists(candidates)]
  if (length(candidates) == 0){
    return(NA_character_)
  }
  normalizePath(candidates[1])
}

source_ivitr_source <- function(source_dir){
  if (is.na(source_dir) || !dir.exists(source_dir)){
    return(FALSE)
  }
  source_files <- sort(list.files(source_dir, pattern = "[.]R$",
                                  full.names = TRUE))
  for (source_file in source_files){
    source(source_file, local = .GlobalEnv)
  }
  TRUE
}

format_call_stack <- function(calls, max_calls = 30L){
  if (length(calls) == 0){
    return("")
  }
  calls <- tail(calls, max_calls)
  call_text <- vapply(calls, function(x){
    paste(deparse(x, width.cutoff = 80), collapse = " ")
  }, character(1))
  paste(seq_along(call_text), call_text, sep = ": ", collapse = "\n")
}

configure_threads <- function(){
  thread_vars <- c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                   "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
                   "NUMEXPR_NUM_THREADS")
  for (var in thread_vars){
    if (!nzchar(Sys.getenv(var, unset = ""))){
      do.call(Sys.setenv, as.list(setNames("1", var)))
    }
  }
}

parse_bool <- function(value, default = FALSE){
  if (!nzchar(value)){
    return(default)
  }
  tolower(value) %in% c("1", "true", "t", "yes", "y")
}

parse_int <- function(value, default){
  out <- suppressWarnings(as.integer(value))
  if (is.na(out) || out < 1){
    return(default)
  }
  out
}

resolve_real_data_path <- function(filename){
  if (dir.exists("real_data")){
    return(file.path("real_data", filename))
  }
  filename
}

configure_threads()

data_path <- Sys.getenv(
  "PCPT_DATA_PATH",
  unset = ""
)
if (!nzchar(data_path) || !file.exists(path.expand(data_path))) {
  stop(
    "Set PCPT_DATA_PATH to an analytic CSV (or run scripts/09_run_synthetic_application.R).",
    call. = FALSE
  )
}
output_file <- Sys.getenv(
  "IVITR_REAL_OUTPUT",
  unset = resolve_real_data_path("analysis_single_85.RData")
)
n_cores <- parse_int(
  Sys.getenv("IVITR_REAL_NCORES", unset = ""),
  max(1L, parallel::detectCores(all.tests = FALSE, logical = TRUE) - 1L)
)
inner_parallel <- parse_bool(Sys.getenv("IVITR_REAL_INNER_PARALLEL", unset = ""),
                             default = FALSE)
parallel_by_d <- parse_bool(Sys.getenv("IVITR_REAL_PARALLEL_BY_D", unset = ""),
                            default = TRUE)
if (parallel_by_d && inner_parallel){
  warning("IVITR_REAL_PARALLEL_BY_D=TRUE and IVITR_REAL_INNER_PARALLEL=TRUE; ",
          "using outer parallelization over d and disabling inner parallelism.")
  inner_parallel <- FALSE
}
nfold <- parse_int(Sys.getenv("IVITR_REAL_NFOLD", unset = ""), 10L)
bandwidth_len <- parse_int(Sys.getenv("IVITR_REAL_BW_LEN", unset = ""), 10L)
prediction_chunk <- parse_int(Sys.getenv("IVITR_REAL_PREDICT_CHUNK",
                                         unset = ""),
                              5000L)
predict_all <- parse_bool(Sys.getenv("IVITR_REAL_PREDICT_ALL", unset = ""),
                          default = TRUE)
diagnostic_n <- parse_int(Sys.getenv("IVITR_REAL_DIAGNOSTIC_N", unset = ""),
                          1000L)
skip_diagnostic <- parse_bool(Sys.getenv("IVITR_REAL_SKIP_DIAGNOSTIC",
                                         unset = ""),
                              default = FALSE)
use_local_ivitr_source <- parse_bool(
  Sys.getenv("IVITR_REAL_SOURCE_LOCAL", unset = ""),
  default = TRUE
)
local_ivitr_source_dir <- find_ivitr_source_dir()
if (use_local_ivitr_source && !is.na(local_ivitr_source_dir)){
  source_ivitr_source(local_ivitr_source_dir)
}

data <- read_csv(data_path, show_col_types = FALSE)

# define covariate
data <- data[stats::complete.cases(data),]
x <- model.matrix(~age+famhx+psa_bl+smoke+bmi, data)
x <- x[,-1]
x[,1] <- scale(x[,1])
x[,6] <- scale(x[,6])

# define outcome
y_tfd <- sapply(data$endPoint, function(t){
  switch(t,
         "Healthy"= 1,
         "Low"= 1,
         "High"= -1)
})

# define randomization
z_tfd <- sapply(data$trt, function(t){
  switch(t,
         "FIN"=1,
         "PLA"= -1)
})

# define treatment received
adh_type_tfd <- sapply(data$adh, function(t){
  switch(t,
         "Adh"= 1,
         "nonAdh"= -1,
         "nonAdh-pill"= -1)
})
trt_tfd <- z_tfd
trt_tfd[z_tfd==1] <- adh_type_tfd[z_tfd==1]

if (skip_diagnostic){
  message("Skipping diagnostic skew plot because IVITR_REAL_SKIP_DIAGNOSTIC=TRUE.")
  skew <- x_range <- y_ccdf <- NULL
} else {
  # Use a subsample for this diagnostic; full-data get_cond_prob is expensive.
  set.seed(1)
  diagnostic_idx <- seq_len(nrow(x))
  if (length(diagnostic_idx) > diagnostic_n){
    diagnostic_idx <- sort(sample(diagnostic_idx, diagnostic_n))
  }
  message("Computing diagnostic skew plot on ", length(diagnostic_idx),
          " subjects.")

  cond_prob <- get_cond_prob(
    y = y_tfd[diagnostic_idx],
    trt = trt_tfd[diagnostic_idx],
    z = z_tfd[diagnostic_idx],
    x = x[diagnostic_idx,, drop = FALSE],
    newx = x[diagnostic_idx,, drop = FALSE]
  )
  pq <- partial_bdd(cond_prob_array = cond_prob, bdd_type = "B-P")
  L_neg <- apply(pq$p,1, function(t){-max(c(t))*(max(c(t))<=0)})
  U_pos <- apply(pq$q,1, function(t){min(c(t))*(min(c(t))>=0)})
  skew <- apply(cbind(L_neg/U_pos,U_pos/L_neg),1,min)
  skew[is.infinite(skew)] <- 0
  f <- ecdf(skew)
  x_range <- seq(0, 1, length.out = 1000)
  y_ccdf <- 1 - f(x_range)
  plot(x_range, y_ccdf,
       type = "s",
       xlim = c(0, 1),
       ylim = c(0, 0.6),
       xlab = "Value of d",
       ylab = "Proportion of Deferred Option",
       col = "red",
       lwd = 2)
}

# run once
set.seed(0)
sampled.idx <- sample(1:1, length(y_tfd), replace = TRUE)
d_seq <- as.numeric(strsplit(
  Sys.getenv("IVITR_REAL_D_SEQ", unset = "0.5,0.3,0.2,0.15,0.12,0.1"),
  "[,[:space:]]+"
)[[1]])
if (!length(d_seq) || any(!is.finite(d_seq)) || any(d_seq <= 0 | d_seq > 0.5)) {
  stop("IVITR_REAL_D_SEQ must contain values in (0, 0.5].", call. = FALSE)
}
index <- 1
fit_n <- parse_int(Sys.getenv("IVITR_REAL_FIT_N", unset = ""),
                   min(2000L, nrow(x)))
fit_idx <- seq_len(nrow(x))
if (length(fit_idx) > fit_n){
  fit_idx <- sort(sample(fit_idx, fit_n))
}
predict_idx <- if (predict_all){
  seq_len(nrow(x))
} else {
  fit_idx
}
y_train <- y_tfd[fit_idx]
trt_train <- trt_tfd[fit_idx]
z_train <- z_tfd[fit_idx]
x_train <- x[fit_idx,, drop = FALSE]
x_predict <- x[predict_idx,, drop = FALSE]
bandwidth_seq <- 10^(seq(-2, 1, length.out = bandwidth_len))
solver <- "CVXR"

predict_ivitr_chunks <- function(fit, newx, chunk_size){
  n_new <- nrow(newx)
  if (n_new == 0){
    return(numeric(0))
  }
  if (n_new <= chunk_size){
    return(predict.iv_itr(fit, newx = newx, type = "link"))
  }

  out <- numeric(n_new)
  starts <- seq(1L, n_new, by = chunk_size)
  for (start in starts){
    end <- min(start + chunk_size - 1L, n_new)
    out[start:end] <- predict.iv_itr(
      fit,
      newx = newx[start:end,, drop = FALSE],
      type = "link"
    )
  }
  out
}

fit_one_d <- function(d.index){
  d_value <- d_seq[d.index]
  start_time <- Sys.time()
  warnings_seen <- character()
  cat(sprintf("[%s] Starting d=%s (%d/%d)\n",
              format(start_time, "%Y-%m-%d %H:%M:%S"),
              d_value, d.index, length(d_seq)))
  flush.console()

  result <- withCallingHandlers(
    tryCatch(
      {
        fit_d <- cf.iv_itr(
          y_train, trt_train, z_train, x_train,
          split_index = NULL,
          d = d_value,
          bdd_type = "B-P",
          bandwidth_seq = bandwidth_seq,
          parallel = inner_parallel,
          solver = solver,
          nfold = nfold
        )
        fit_pred_d <- predict_ivitr_chunks(fit_d, x_predict, prediction_chunk)
        list(fit = fit_d, fit_pred = fit_pred_d, status = "fit")
      },
      error = function(e){
        warning_text <- if (length(warnings_seen) > 0){
          paste("\nWarnings before failure:\n",
                paste(unique(warnings_seen), collapse = "\n"), sep = "")
        } else {
          ""
        }
        message(sprintf("d=%s failed: %s%s\nUsing an all-defer fallback ",
                        d_value, conditionMessage(e), warning_text),
                "so the remaining thresholds can be analysed.")
        # A zero link produces the deferred decision for every d < 0.5.
        # The d=0.5 fit is normally stable; its no-defer convention means an
        # all-defer action is not representable on that boundary.
        list(fit = NULL, fit_pred = rep(0, nrow(x_predict)),
             status = "fallback_zero_link")
      }
    ),
    warning = function(w){
      warning_message <- sprintf("[warning][d=%s] %s",
                                 d_value, conditionMessage(w))
      warnings_seen <<- c(warnings_seen, warning_message)
      message(warning_message)
      invokeRestart("muffleWarning")
    }
  )

  end_time <- Sys.time()
  cat(sprintf("[%s] Finished d=%s | elapsed=%.2f sec\n",
              format(end_time, "%Y-%m-%d %H:%M:%S"),
              d_value,
              as.numeric(difftime(end_time, start_time, units = "secs"))))
  flush.console()

  list(d.index = d.index, d = d_value,
       fit = result$fit, fit_pred = result$fit_pred, status = result$status)
}

outer_cores <- if (parallel_by_d){
  min(n_cores, length(d_seq))
} else {
  1L
}

cat("analysis_single.R settings:\n")
cat("  data_path=", data_path, "\n", sep = "")
cat("  output_file=", output_file, "\n", sep = "")
cat("  d_seq=", paste(d_seq, collapse = ", "), "\n", sep = "")
cat("  n_complete=", nrow(x), "\n", sep = "")
cat("  fit_n=", length(fit_idx), "\n", sep = "")
cat("  predict_n=", length(predict_idx), "\n", sep = "")
cat("  requested_cores=", n_cores, "\n", sep = "")
cat("  outer_cores=", outer_cores, "\n", sep = "")
cat("  parallel_by_d=", parallel_by_d, "\n", sep = "")
cat("  inner_parallel=", inner_parallel, "\n", sep = "")
cat("  nfold=", nfold, "\n", sep = "")
cat("  bandwidth_len=", bandwidth_len, "\n", sep = "")
cat("  predict_all=", predict_all, "\n", sep = "")
cat("  prediction_chunk=", prediction_chunk, "\n", sep = "")
cat("  skip_diagnostic=", skip_diagnostic, "\n", sep = "")
cat("  diagnostic_n=", diagnostic_n, "\n", sep = "")
cat("  source_local_ivitr=", use_local_ivitr_source, "\n", sep = "")
cat("  local_ivitr_source_dir=", local_ivitr_source_dir, "\n", sep = "")

if (outer_cores > 1L){
  cl <- parallel::makeCluster(outer_cores, outfile = "")
  doParallel::registerDoParallel(cl)
  worker_exports <- c(
    "fit_one_d", "d_seq", "y_train", "trt_train", "z_train", "x_train",
    "x_predict", "bandwidth_seq", "inner_parallel", "solver", "nfold",
    "prediction_chunk", "predict_ivitr_chunks", "format_call_stack",
    "use_local_ivitr_source", "local_ivitr_source_dir", "source_ivitr_source"
  )
  parallel::clusterExport(cl, varlist = worker_exports, envir = environment())
  fit_results <- tryCatch(
    {
      foreach(
        d.index = seq_along(d_seq),
        .packages = "IVITR",
        .errorhandling = "stop"
      ) %dopar% {
        if (use_local_ivitr_source){
          source_ivitr_source(local_ivitr_source_dir)
        }
        fit_one_d(d.index)
      }
    },
    finally = {
      parallel::stopCluster(cl)
      foreach::registerDoSEQ()
    }
  )
} else {
  foreach::registerDoSEQ()
  fit_results <- lapply(seq_along(d_seq), fit_one_d)
}

fit <- vector("list", length(d_seq))
fit_pred <- vector("list", length(d_seq))
fit_status <- character(length(d_seq))
for (result in fit_results){
  fit[[result$d.index]] <- result$fit
  fit_pred[[result$d.index]] <- result$fit_pred
  fit_status[result$d.index] <- result$status
}
names(fit) <- names(fit_pred) <- paste0("d_", d_seq)

output_dir <- dirname(output_file)
if (nzchar(output_dir) && output_dir != "."){
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
}
save(fit_pred, fit, fit_status, d_seq, sampled.idx, index, fit_idx, predict_idx,
     file = output_file)
