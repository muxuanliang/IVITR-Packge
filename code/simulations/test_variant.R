sim_cores <- function(){
  n_cores <- as.integer(Sys.getenv("IVITR_NCORES", unset=NA_character_))
  if (is.na(n_cores)){
    n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset=NA_character_))
  }
  if (is.na(n_cores)){
    n_cores <- as.integer(Sys.getenv("SLURM_CPUS_ON_NODE", unset=NA_character_))
  }
  if (is.na(n_cores)){
    n_cores <- as.integer(Sys.getenv("LSB_DJOB_NUMPROC", unset=NA_character_))
  }
  if (is.na(n_cores)){
    n_cores <- parallel::detectCores(all.tests = FALSE, logical = TRUE)
  }
  if (is.na(n_cores) || n_cores < 1){
    n_cores <- 1
  }
  n_cores
}

configure_hpc_threads <- function(){
  thread_vars <- c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                   "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
                   "NUMEXPR_NUM_THREADS")
  for (var in thread_vars){
    if (!nzchar(Sys.getenv(var, unset=""))){
      do.call(Sys.setenv, as.list(setNames("1", var)))
    }
  }
}

parse_num_seq <- function(value, default){
  if (!nzchar(value)){
    return(default)
  }
  as.numeric(strsplit(value, "[,[:space:]]+")[[1]])
}

parse_chr_seq <- function(value, default){
  if (!nzchar(value)){
    return(default)
  }
  out <- strsplit(value, "[,[:space:]]+")[[1]]
  out[nzchar(out)]
}

safe_label <- function(value){
  gsub("[^A-Za-z0-9]+", "-", as.character(value))
}

resolve_result_dir <- function(result_dir){
  if (grepl("^(/|~)", result_dir)){
    return(path.expand(result_dir))
  }
  file.path(getwd(), result_dir)
}

numeric_or_named <- function(level, named_values){
  numeric_value <- suppressWarnings(as.numeric(level))
  if (!is.na(numeric_value)){
    return(numeric_value)
  }
  if (!level %in% names(named_values)){
    stop("Unknown level: ", level, ". Allowed levels are ",
         paste(names(named_values), collapse=", "),
         " or a numeric value.")
  }
  unname(named_values[level])
}

za_coef_from_level <- function(level){
  numeric_or_named(level, c(strong=12, existing=8, weak=4))
}

p_from_level <- function(level){
  as.integer(numeric_or_named(level, c(existing=5, more=10)))
}

u_a_multiplier_from_level <- function(level){
  numeric_or_named(level, c(existing=1, weak=0.5, none=0))
}

variant_specs <- function(){
  za_levels <- parse_chr_seq(Sys.getenv("IVITR_ZA_LEVEL_SEQ", unset=""),
                             c("strong", "existing", "weak"))
  p_levels <- parse_chr_seq(Sys.getenv("IVITR_P_LEVEL_SEQ", unset=""),
                            c("existing", "more"))
  u_a_levels <- parse_chr_seq(Sys.getenv("IVITR_UA_LEVEL_SEQ", unset=""),
                              c("existing", "weak", "none"))

  za_levels <- za_levels[za_levels != "existing"]
  p_levels <- p_levels[p_levels != "existing"]
  u_a_levels <- u_a_levels[u_a_levels != "existing"]

  rbind(
    data.frame(variant_family="baseline", variant_level="existing",
               za_level="existing", p_level="existing",
               u_a_level="existing", stringsAsFactors=FALSE),
    data.frame(variant_family="z_a", variant_level=za_levels,
               za_level=za_levels, p_level="existing",
               u_a_level="existing", stringsAsFactors=FALSE),
    data.frame(variant_family="covariates", variant_level=p_levels,
               za_level="existing", p_level=p_levels,
               u_a_level="existing", stringsAsFactors=FALSE),
    data.frame(variant_family="u_a", variant_level=u_a_levels,
               za_level="existing", p_level="existing",
               u_a_level=u_a_levels, stringsAsFactors=FALSE)
  )
}

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

simulation_grid <- function(){
  base_grid <- expand.grid(
    setting1=parse_num_seq(Sys.getenv("IVITR_SETTING1_SEQ", unset=""), 2),
    setting2=parse_num_seq(Sys.getenv("IVITR_SETTING2_SEQ", unset=""), 2),
    n=parse_num_seq(Sys.getenv("IVITR_N_SEQ", unset=""), c(300, 500)),
    eta=parse_num_seq(Sys.getenv("IVITR_ETA_SEQ", unset=""), 1.5),
    d=parse_num_seq(Sys.getenv("IVITR_D_SEQ", unset=""), c(0.5, 0.4, 0.3, 0.2, 0.1)),
    KEEP.OUT.ATTRS=FALSE
  )
  specs <- variant_specs()
  do.call(rbind, lapply(seq_len(nrow(specs)), function(i){
    out <- cbind(base_grid, specs[rep(i, nrow(base_grid)),, drop=FALSE])
    rownames(out) <- NULL
    out
  }))
}

selected_task_id <- function(args=commandArgs(trailingOnly=TRUE)){
  task_id <- get_arg_value(args, "task-id", Sys.getenv("SLURM_ARRAY_TASK_ID", unset=""))
  if (!nzchar(task_id)){
    task_id <- Sys.getenv("LSB_JOBINDEX", unset="")
  }
  if (!nzchar(task_id)){
    return(NA_integer_)
  }
  as.integer(task_id)
}

draw_treatment <- function(x, u, z, lambda, za_coef){
  rbinom(length(u), 1, prob=plogis(za_coef * z + x[,1] - 7 * x[,2] +
                                     lambda * (1 + x[,1]) * u))
}

outcome_linear_predictor <- function(x, u, trt, xi, delta, setting1, setting2){
  x1 <- x[,1]
  x2 <- x[,2]
  if (setting1 == 1){
    g1 <- 1 - x1 - x2 + xi * u
  } else {
    g1 <- 1 - x1^2 + x2^2 + xi * x1 * x2 * u
  }
  if (setting2 == 1){
    g2 <- 0.442 * (1 - x1 - x2 + delta * u) * trt
  } else {
    g2 <- (x2 - delta * 0.25 * x1^2 + 4 * (x2 > 0) - 6 * (x2 <= 0)) * trt
  }
  g1 + g2
}

draw_outcome <- function(x, u, trt, xi, delta, setting1, setting2){
  rbinom(length(u), 1,
         prob=plogis(outcome_linear_predictor(x, u, trt, xi, delta,
                                              setting1, setting2)))
}

mean_outcome <- function(x, u, trt, xi, delta, setting1, setting2){
  plogis(outcome_linear_predictor(x, u, trt, xi, delta, setting1, setting2))
}

cond_mean_outcome <- function(x_mat, a, xi, delta, setting1, setting2){
  x1 <- x_mat[,1]
  x2 <- x_mat[,2]
  if (setting1 == 1){
    g1_int <- 1 - x1 - x2
    g1_slope <- rep(xi, length(x1))
  } else {
    g1_int <- 1 - x1^2 + x2^2
    g1_slope <- xi * x1 * x2
  }
  if (setting2 == 1){
    g2_int <- 0.442 * (1 - x1 - x2) * a
    g2_slope <- 0.442 * delta * a
  } else {
    g2_int <- (x2 - delta * 0.25 * x1^2 + 4 * (x2 > 0) - 6 * (x2 <= 0)) * a
    g2_slope <- 0
  }

  intercept <- g1_int + g2_int
  slope <- g1_slope + g2_slope
  log1pexp <- function(v){
    out <- numeric(length(v))
    pos <- v > 0
    out[pos] <- v[pos] + log1p(exp(-v[pos]))
    out[!pos] <- log1p(exp(v[!pos]))
    out
  }
  cond_mean <- (log1pexp(intercept + slope) - log1pexp(intercept - slope)) / (2 * slope)
  near_zero <- abs(slope) < 1e-8
  cond_mean[near_zero] <- plogis(intercept[near_zero])
  cond_mean
}

value_nonrejected <- function(y1, y0, decision){
  nonrejected <- decision != 0
  if (!any(nonrejected)){
    return(NA_real_)
  }
  sum(y1 * (decision == 1) + y0 * (decision == -1)) / sum(nonrejected)
}

bdd_loss_from_decision <- function(decision, p, q, d){
  p_max <- if (is.null(dim(p))) as.numeric(p) else apply(p, 1, max)
  q_min <- if (is.null(dim(q))) as.numeric(q) else apply(q, 1, min)
  mean(d * (q_min - p_max) * (q_min > 0 & p_max < 0) * (decision == 0) +
         (decision == 1) * (p_max < 0) * abs(p_max) +
         (decision == -1) * (q_min > 0) * abs(q_min))
}

sim <- function(n=500, eta=0.5, setting1=1, setting2=2, d=0.5,
                variant_family="z_a", variant_level="existing",
                za_level="existing", p_level="existing",
                u_a_level="existing",
                n_rep=as.integer(Sys.getenv("IVITR_N_REP", unset="200")),
                n_eval=as.integer(Sys.getenv("IVITR_N_EVAL", unset="1000")),
                result_dir=Sys.getenv("IVITR_VARIANT_DIR",
                                      unset="simulation_R2_variant")){
  configure_hpc_threads()
  library(doParallel)
  n_cores <- sim_cores()
  cl <- NULL
  if (n_cores > 1){
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    on.exit(stopCluster(cl), add = TRUE)
  } else {
    foreach::registerDoSEQ()
  }

  p <- p_from_level(p_level)
  za_coef <- za_coef_from_level(za_level)
  u_a_multiplier <- u_a_multiplier_from_level(u_a_level)
  lambda <- eta * u_a_multiplier
  xi <- delta <- eta

  helper_exports <- c("draw_treatment", "outcome_linear_predictor", "draw_outcome",
                      "mean_outcome", "cond_mean_outcome", "value_nonrejected",
                      "bdd_loss_from_decision")
  print(system.time(res <- foreach(i=1:n_rep, .packages="IVITR",
                                   .export=helper_exports,
                                   .errorhandling = "stop",.combine = rbind)%dopar%{
    set.seed(i)

    x <- array(runif(n*p, min=-1, max=1), c(n, p))
    u <- runif(n, min=-1, max=1)
    z <- (rnorm(n)>0)
    trt <- draw_treatment(x, u, z, lambda, za_coef)
    y <- draw_outcome(x, u, trt, xi, delta, setting1, setting2)

    # flip 0,1 to -1, 1
    z_tfd <- 2*(z-0.5)
    trt_tfd <- 2*(trt-0.5)
    y_tfd <- 2*(y-0.5)

    # fit
    fit <- cf.iv_itr(y_tfd, trt_tfd,z_tfd,x, split_index = NULL, d=d, bdd_type="B-P")
    fit0 <- cf.iv_itr(y_tfd, trt_tfd,z_tfd,x, split_index = NULL, d=d, bdd_type="B-P", debias = FALSE)

    # new data
    newx <- array(runif(n_eval*p, min=-1, max=1), c(n_eval, p))
    newu <- runif(n_eval, min=-1, max=1)
    y1 <- mean_outcome(newx, newu, rep(1, n_eval), xi, delta, setting1, setting2)
    y0 <- mean_outcome(newx, newu, rep(0, n_eval), xi, delta, setting1, setting2)
    true_effect_x <- cond_mean_outcome(newx, 1, xi, delta, setting1, setting2) -
      cond_mean_outcome(newx, 0, xi, delta, setting1, setting2)

    fit_pred <- predict.iv_itr(fit, newx=newx, type="link")
    fit0_pred <- predict.iv_itr(fit0, newx=newx, type="link")

    fit_decision <- defer_decision_from_link(fit_pred, d=d)
    fit0_decision <- defer_decision_from_link(fit0_pred, d=d)

    # value on non-rejected
    value_fit_nonreject <- value_nonrejected(y1, y0, fit_decision)
    value_fit0_nonreject <- value_nonrejected(y1, y0, fit0_decision)

    # proportion rejected
    pr_fit <- mean(fit_decision == 0)
    pr_fit0 <- mean(fit0_decision == 0)

    # loss
    newz <- (rnorm(n_eval)>0)
    newtrt <- draw_treatment(newx, newu, newz, lambda, za_coef)
    newy <- draw_outcome(newx, newu, newtrt, xi, delta, setting1, setting2)
    newz_tfd <- 2*(newz-0.5)
    newtrt_tfd <- 2*(newtrt-0.5)
    newy_tfd <- 2*(newy-0.5)

    cond_prob <- get_cond_prob(y=newy_tfd, trt=newtrt_tfd, z=newz_tfd, x=newx)
    prob_bdd <- partial_bdd(cond_prob, bdd_type = "B-P")

    review_fit <- defer_review_proportions(
      estimated = fit_decision, true_effect = true_effect_x,
      p = prob_bdd$p, q = prob_bdd$q,
      d = d, estimated_type = "decision"
    )
    review_fit0 <- defer_review_proportions(
      estimated = fit0_decision, true_effect = true_effect_x,
      p = prob_bdd$p, q = prob_bdd$q,
      d = d, estimated_type = "decision"
    )

    loss_fit <- bdd_loss_from_decision(fit_decision, prob_bdd$p, prob_bdd$q, d)
    loss_fit0 <- bdd_loss_from_decision(fit0_decision, prob_bdd$p, prob_bdd$q, d)

    res <- data.frame(value=c(value_fit_nonreject, value_fit0_nonreject,
                              pr_fit, pr_fit0,
                              loss_fit, loss_fit0,
                              review_fit["prop_bounds_cross_zero_among_deferred"],
                              review_fit0["prop_bounds_cross_zero_among_deferred"],
                              review_fit["prop_wrong_treatment_among_nondeferred"],
                              review_fit0["prop_wrong_treatment_among_nondeferred"]),
                      d=d, metric=rep(c("value on non-rejected", "porp rejected", "loss in value",
                                        "prop deferred with bounds crossing zero",
                                        "prop wrong treatment among non-deferred"), each=2),
                      method=c("Debiased", "No-debias"), index=i, n=n, lambda=lambda, xi=xi, delta=delta,
                      setting1=setting1, setting2=setting2,
                      variant_family=variant_family,
                      variant_level=variant_level,
                      za_level=za_level,
                      za_coef=za_coef,
                      p_level=p_level,
                      p_covariates=p,
                      u_a_level=u_a_level,
                      u_a_multiplier=u_a_multiplier)
    res
  }))
  expected_rows <- n_rep * 10
  if (is.null(res) || nrow(res) != expected_rows){
    stop("Expected ", expected_rows, " result rows, but got ",
         if (is.null(res)) 0 else nrow(res), ".")
  }
  result_path <- resolve_result_dir(result_dir)
  dir.create(result_path, showWarnings = FALSE, recursive = TRUE)
  output_name <- paste0("res_variant_",
                        safe_label(variant_family), "_",
                        safe_label(variant_level),
                        "_n", n,
                        "_eta", eta,
                        "_s", setting1, "_", setting2,
                        "_d", d, ".csv")
  write.csv(res, file=file.path(result_path, output_name),
            row.names=FALSE)
}

time_sim <- function(n=500, eta=0.5, setting1=1, setting2=2, d=0.5,
                     variant_family="z_a", variant_level="existing",
                     za_level="existing", p_level="existing",
                     u_a_level="existing",
                     result_dir=Sys.getenv("IVITR_VARIANT_DIR",
                                           unset="simulation_R2_variant"),
                     timing_file=file.path(resolve_result_dir(result_dir),
                                           "runtime_log_variant.csv")){
  dir.create(dirname(timing_file), showWarnings = FALSE, recursive = TRUE)
  za_coef <- za_coef_from_level(za_level)
  p_covariates <- p_from_level(p_level)
  u_a_multiplier <- u_a_multiplier_from_level(u_a_level)
  lambda <- eta * u_a_multiplier

  start_time <- Sys.time()
  setting_label <- paste0("n=", n, ", eta=", eta,
                          ", setting=(", setting1, ",", setting2, ")",
                          ", d=", d,
                          ", variant=", variant_family, ":", variant_level,
                          ", za=", za_level, "(", za_coef, ")",
                          ", p=", p_covariates,
                          ", u_a=", u_a_level, "(", lambda, ")")
  cat(sprintf("[%s] Starting variant sim: %s\n",
              format(start_time, "%Y-%m-%d %H:%M:%S %Z"), setting_label))
  flush.console()

  error_message <- NA_character_
  timing <- system.time(
    status <- tryCatch(
      {
        sim(n=n, eta=eta, setting1=setting1, setting2=setting2, d=d,
            variant_family=variant_family, variant_level=variant_level,
            za_level=za_level, p_level=p_level, u_a_level=u_a_level,
            result_dir=result_dir)
        "success"
      },
      error=function(e){
        error_message <<- conditionMessage(e)
        "error"
      }
    )
  )
  end_time <- Sys.time()

  cat(sprintf("[%s] Finished variant sim: %s | status=%s | elapsed=%.2f sec\n",
              format(end_time, "%Y-%m-%d %H:%M:%S %Z"), setting_label,
              status, unname(timing[["elapsed"]])))
  flush.console()

  runtime <- data.frame(
    start_time=format(start_time, "%Y-%m-%d %H:%M:%S %Z"),
    end_time=format(end_time, "%Y-%m-%d %H:%M:%S %Z"),
    elapsed_sec=unname(timing[["elapsed"]]),
    user_sec=unname(timing[["user.self"]]),
    system_sec=unname(timing[["sys.self"]]),
    status=status,
    error=error_message,
    n=n, eta=eta, lambda=lambda, xi=eta, delta=eta,
    setting1=setting1, setting2=setting2, d=d,
    variant_family=variant_family,
    variant_level=variant_level,
    za_level=za_level,
    za_coef=za_coef,
    p_level=p_level,
    p_covariates=p_covariates,
    u_a_level=u_a_level,
    u_a_multiplier=u_a_multiplier
  )

  write.table(runtime, file=timing_file, sep=",", row.names=FALSE,
              col.names=!file.exists(timing_file), append=file.exists(timing_file))
  if (status != "success"){
    stop(error_message, call. = FALSE)
  }
  invisible(runtime)
}

run_driver <- function(args=commandArgs(trailingOnly=TRUE)){
  grid <- simulation_grid()
  task_id <- selected_task_id(args)

  cat("Variant simulation grid has", nrow(grid), "tasks.\n")
  if (!is.na(task_id)){
    if (task_id < 1 || task_id > nrow(grid)){
      stop("task-id must be between 1 and ", nrow(grid), ".")
    }
    grid <- grid[task_id,, drop=FALSE]
    cat("Running task", task_id, "only.\n")
  }
  cat("IVITR_NCORES=", sim_cores(), "\n", sep="")
  cat("IVITR_N_REP=", Sys.getenv("IVITR_N_REP", unset="200"), "\n", sep="")
  cat("IVITR_N_EVAL=", Sys.getenv("IVITR_N_EVAL", unset="1000"), "\n", sep="")
  cat("IVITR_VARIANT_DIR=", Sys.getenv("IVITR_VARIANT_DIR",
                                       unset="simulation_R2_variant"),
      "\n", sep="")
  cat("Variant design=marginal_only\n")
  cat("Variant families:",
      paste(unique(grid$variant_family), collapse=", "), "\n")

  for (row_idx in seq_len(nrow(grid))){
    time_sim(
      n = grid$n[row_idx],
      eta = grid$eta[row_idx],
      setting1 = grid$setting1[row_idx],
      setting2 = grid$setting2[row_idx],
      d = grid$d[row_idx],
      variant_family = grid$variant_family[row_idx],
      variant_level = grid$variant_level[row_idx],
      za_level = grid$za_level[row_idx],
      p_level = grid$p_level[row_idx],
      u_a_level = grid$u_a_level[row_idx]
    )
  }
}

if (!interactive()){
  run_driver()
}
