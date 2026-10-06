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
  expand.grid(
    setting1=parse_num_seq(Sys.getenv("IVITR_SETTING1_SEQ", unset=""), c(1, 2)),
    setting2=parse_num_seq(Sys.getenv("IVITR_SETTING2_SEQ", unset=""), c(1, 2)),
    n=parse_num_seq(Sys.getenv("IVITR_N_SEQ", unset=""), c(300, 500)),
    eta=parse_num_seq(Sys.getenv("IVITR_ETA_SEQ", unset=""), c(0.5, 1.5)),
    d=parse_num_seq(Sys.getenv("IVITR_D_SEQ", unset=""), c(0.5, 0.4, 0.3, 0.2, 0.1)),
    KEEP.OUT.ATTRS=FALSE
  )
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

draw_treatment <- function(x, u, z, lambda){
  rbinom(length(u), 1, prob=plogis(8 * z + x[,1] - 7 * x[,2] +
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

sim <- function(n=500, eta=0.5, setting1=1, setting2=1, d=0.5,
                n_rep=as.integer(Sys.getenv("IVITR_N_REP", unset="200")),
                n_eval=as.integer(Sys.getenv("IVITR_N_EVAL", unset="1000"))){
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

  p <- 5
  lambda <- xi <- delta <- eta

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
    trt <- draw_treatment(x, u, z, lambda)
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
    newtrt <- draw_treatment(newx, newu, newz, lambda)
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
                      setting1=setting1, setting2=setting2)
    res
  }))
  expected_rows <- n_rep * 10
  if (is.null(res) || nrow(res) != expected_rows){
    stop("Expected ", expected_rows, " result rows, but got ",
         if (is.null(res)) 0 else nrow(res), ".")
  }
  dir.create(file.path(getwd(), "simulation_R2"), showWarnings = FALSE, recursive = TRUE)
  write.csv(res, file=paste0(getwd(),"/simulation_R2/res_",n,"_",lambda,"_", xi, "_", delta, "_",setting1, "_", setting2, "_", d,".csv"))
}

time_sim <- function(n=500, eta=0.5, setting1=1, setting2=1, d=0.5,
                     timing_file=file.path(getwd(), "simulation_R2", "runtime_log.csv")){
  dir.create(dirname(timing_file), showWarnings = FALSE, recursive = TRUE)

  start_time <- Sys.time()
  setting_label <- paste0("n=", n, ", eta=", eta, ", setting=(", setting1, ",", setting2, "), d=", d)
  cat(sprintf("[%s] Starting sim: %s\n",
              format(start_time, "%Y-%m-%d %H:%M:%S %Z"), setting_label))
  flush.console()

  error_message <- NA_character_
  timing <- system.time(
    status <- tryCatch(
      {
        sim(n=n, eta=eta, setting1=setting1, setting2=setting2, d=d)
        "success"
      },
      error=function(e){
        error_message <<- conditionMessage(e)
        "error"
      }
    )
  )
  end_time <- Sys.time()

  cat(sprintf("[%s] Finished sim: %s | status=%s | elapsed=%.2f sec\n",
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
    n=n, eta=eta, lambda=eta, xi=eta, delta=eta,
    setting1=setting1, setting2=setting2, d=d
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

  cat("Main simulation grid has", nrow(grid), "tasks.\n")
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

  for (row_idx in seq_len(nrow(grid))){
    time_sim(
      n = grid$n[row_idx],
      eta = grid$eta[row_idx],
      setting1 = grid$setting1[row_idx],
      setting2 = grid$setting2[row_idx],
      d = grid$d[row_idx]
    )
  }
}

if (!interactive()){
  run_driver()
}
