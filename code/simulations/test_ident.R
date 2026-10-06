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

resolve_result_dir <- function(result_dir){
  if (grepl("^(/|~)", result_dir)){
    return(path.expand(result_dir))
  }
  file.path(getwd(), result_dir)
}

safe_label <- function(value){
  gsub("[^A-Za-z0-9]+", "-", as.character(value))
}

ident_result_name <- function(n, eta, setting1, setting2, p_covariates,
                              za_coef, x1_coef, x2_coef, d){
  paste0("res_ident_n", n,
         "_eta", eta,
         "_s", setting1, "_", setting2,
         "_p", p_covariates,
         "_z", safe_label(za_coef),
         "_x1", safe_label(x1_coef),
         "_x2", safe_label(x2_coef),
         "_d", d,
         ".csv")
}

simulation_grid <- function(){
  expand.grid(
    setting1=parse_num_seq(Sys.getenv("IVITR_SETTING1_SEQ", unset=""), c(1, 2)),
    setting2=parse_num_seq(Sys.getenv("IVITR_SETTING2_SEQ", unset=""), c(1, 2)),
    n=parse_num_seq(Sys.getenv("IVITR_N_SEQ", unset=""), c(300, 500)),
    eta=parse_num_seq(Sys.getenv("IVITR_ETA_SEQ", unset=""), c(0.5, 1.5)),
    d=parse_num_seq(Sys.getenv("IVITR_D_SEQ", unset=""),
                    c(0.5, 0.4, 0.3, 0.2, 0.1)),
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

first_stage_shift <- function(x, x1_coef, x2_coef){
  x1_coef * x[,1] - x2_coef * x[,2]
}

treatment_prob <- function(x, z_tfd, za_coef, x1_coef, x2_coef){
  plogis(za_coef * z_tfd + first_stage_shift(x, x1_coef, x2_coef))
}

first_stage_delta <- function(x, za_coef, x1_coef, x2_coef){
  treatment_prob(x, rep(1, nrow(x)), za_coef, x1_coef, x2_coef) -
    treatment_prob(x, rep(-1, nrow(x)), za_coef, x1_coef, x2_coef)
}

first_stage_bounds <- function(za_coef, x1_coef, x2_coef){
  max_shift <- abs(x1_coef) + abs(x2_coef)
  shifts <- c(-max_shift, 0, max_shift)
  delta <- plogis(za_coef + shifts) - plogis(-za_coef + shifts)
  c(min=min(delta), max=max(delta))
}

draw_treatment <- function(x, z_tfd, za_coef, x1_coef, x2_coef){
  rbinom(length(z_tfd), 1,
         prob=treatment_prob(x, z_tfd, za_coef, x1_coef, x2_coef))
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
    g2 <- (x2 - delta * 0.25 * x1^2 +
             4 * (x2 > 0) - 6 * (x2 <= 0)) * trt
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
    g2_int <- (x2 - delta * 0.25 * x1^2 +
                 4 * (x2 > 0) - 6 * (x2 <= 0)) * a
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
  cond_mean <- (log1pexp(intercept + slope) -
                  log1pexp(intercept - slope)) / (2 * slope)
  near_zero <- abs(slope) < 1e-8
  cond_mean[near_zero] <- plogis(intercept[near_zero])
  cond_mean
}

decision_from_effect <- function(effect){
  decision <- ifelse(effect >= 0, 1L, -1L)
  decision[!is.finite(effect)] <- NA_integer_
  as.integer(decision)
}

value_nonrejected <- function(y1, y0, decision){
  observed <- !is.na(decision)
  nonrejected <- observed & decision != 0
  if (!any(nonrejected)){
    return(NA_real_)
  }
  mean(y1[nonrejected] * (decision[nonrejected] == 1) +
         y0[nonrejected] * (decision[nonrejected] == -1))
}

value_population <- function(y1, y0, decision, defer_value=0){
  value <- rep(defer_value, length(decision))
  treat <- !is.na(decision) & decision == 1
  control <- !is.na(decision) & decision == -1
  value[treat] <- y1[treat]
  value[control] <- y0[control]
  mean(value[!is.na(value)])
}

wrong_treatment_prop <- function(decision, true_effect){
  observed <- !is.na(decision) & decision != 0 & is.finite(true_effect)
  if (!any(observed)){
    return(NA_real_)
  }
  mean(decision[observed] * true_effect[observed] < 0)
}

effect_metric <- function(estimated, truth, type=c("bias", "mae", "mse")){
  type <- match.arg(type)
  observed <- is.finite(estimated) & is.finite(truth)
  if (!any(observed)){
    return(NA_real_)
  }
  error <- estimated[observed] - truth[observed]
  switch(type,
         bias=mean(error),
         mae=mean(abs(error)),
         mse=mean(error^2))
}

# These are the exact conditional probabilities under the identifiable DGP.
# They are used only to evaluate the population bound/risk, not to fit a method.
true_cond_prob <- function(x, xi, delta, setting1, setting2,
                           za_coef, x1_coef, x2_coef){
  n_eval <- nrow(x)
  mu0 <- cond_mean_outcome(x, 0, xi, delta, setting1, setting2)
  mu1 <- cond_mean_outcome(x, 1, xi, delta, setting1, setting2)
  cond_prob <- array(0, c(2, 2, 2, n_eval))

  for (z_index in 1:2){
    z_tfd <- if (z_index == 1) -1 else 1
    prob_a1 <- treatment_prob(x, rep(z_tfd, n_eval), za_coef, x1_coef, x2_coef)
    for (a_index in 1:2){
      prob_a <- if (a_index == 1) 1 - prob_a1 else prob_a1
      mu_a <- if (a_index == 1) mu0 else mu1
      cond_prob[1, a_index, z_index,] <- prob_a * (1 - mu_a)
      cond_prob[2, a_index, z_index,] <- prob_a * mu_a
    }
  }
  cond_prob
}

loss_components <- function(decision, lower, upper, true_effect, d){
  observed <- !is.na(decision) & is.finite(lower) & is.finite(upper) &
    is.finite(true_effect)
  if (!any(observed)){
    return(c(worst_case=NA_real_, true_cate=NA_real_))
  }

  decision <- decision[observed]
  lower <- lower[observed]
  upper <- upper[observed]
  true_effect <- true_effect[observed]
  crosses_zero <- lower < 0 & upper > 0
  deferred <- decision == 0

  # The implementation's d is the manuscript rejection-cost parameter after
  # reparameterization.  The defer term is evaluated on the same ambiguous
  # (bound-crossing) stratum as value_bdd, with the true CATE replacing the
  # worst-case interval width.
  worst_case <- d * (upper - lower) * crosses_zero * deferred +
    (decision == 1) * pmax(-lower, 0) +
    (decision == -1) * pmax(upper, 0)
  true_cate <- d * abs(true_effect) * crosses_zero * deferred +
    (decision == 1) * pmax(-true_effect, 0) +
    (decision == -1) * pmax(true_effect, 0)

  # With population conditional probabilities, tau(X) belongs to [lower, upper].
  # This verifies that the reported bound is genuinely a worst-case loss.
  if (any(true_cate > worst_case + 1e-10)){
    stop("The worst-case loss is smaller than the true-CATE loss; check the d reparameterization or probability bounds.")
  }
  c(worst_case=mean(worst_case), true_cate=mean(true_cate))
}

safe_proposed_decision <- function(fit_fn, y, trt, z, x, newx, d,
                                   debias=TRUE, predict_fn=predict.iv_itr){
  tryCatch({
    fit <- fit_fn(y, trt, z, x, split_index=NULL, d=d, bdd_type="B-P",
                  debias=debias)
    link <- predict_fn(fit, newx=newx, type="link")
    decision <- defer_decision_from_link(link, d=d)
    if (length(decision) != nrow(newx) || anyNA(decision) ||
        !all(decision %in% c(-1L, 0L, 1L))){
      stop("The optimizer returned an invalid prediction.")
    }
    list(decision=as.integer(decision), status="fit", error=NA_character_)
  }, error=function(e){
    list(decision=rep.int(0L, nrow(newx)), status="fallback_all_defer",
         error=gsub("[\r\n]+", " ", conditionMessage(e)))
  })
}

metric_rows <- function(method_values, d, i, n, eta, xi, delta,
                        setting1, setting2, p_covariates, za_coef, x1_coef, x2_coef,
                        first_stage_mean, defer_value){
  methods <- names(method_values$population_value)
  metrics <- c("population value",
               "value on non-deferred",
               "value loss vs ad-hoc",
               "worst-case loss in value",
               "true-CATE loss in value",
               "prop deferred",
               "prop wrong treatment among non-deferred",
               "effect bias",
               "effect mae",
               "effect mse")

  values <- c(
    method_values$population_value,
    method_values$nondeferred_value,
    method_values$population_value["Ad-hoc"] - method_values$population_value,
    method_values$worst_case_loss,
    method_values$true_cate_loss,
    method_values$prop_deferred,
    method_values$prop_wrong,
    method_values$effect_bias,
    method_values$effect_mae,
    method_values$effect_mse
  )

  data.frame(
    value=as.numeric(values),
    d=d,
    metric=rep(metrics, each=length(methods)),
    method=rep(methods, times=length(metrics)),
    index=i,
    n=n,
    eta=eta,
    lambda=0,
    xi=xi,
    delta=delta,
    setting1=setting1,
    setting2=setting2,
    scenario="identifiable_no_u_to_a",
    p_covariates=p_covariates,
    za_coef=za_coef,
    x1_to_a_coef=x1_coef,
    x2_to_a_coef=x2_coef,
    first_stage_at_x0=plogis(za_coef) - plogis(-za_coef),
    first_stage_mean=first_stage_mean,
    first_stage_min=unname(first_stage_bounds(za_coef, x1_coef, x2_coef)["min"]),
    first_stage_max=unname(first_stage_bounds(za_coef, x1_coef, x2_coef)["max"]),
    u_a_multiplier=0,
    defer_value=defer_value,
    optimization_status=rep(method_values$optimization_status,
                            times=length(metrics)),
    optimization_error=rep(method_values$optimization_error,
                           times=length(metrics)),
    stringsAsFactors=FALSE
  )
}

sim <- function(n=500, eta=0.5, setting1=1, setting2=1, d=0.5,
                p_covariates=as.integer(Sys.getenv("IVITR_IDENT_P", unset="5")),
                za_coef=as.numeric(Sys.getenv("IVITR_IDENT_ZA", unset="8")),
                x1_coef=as.numeric(Sys.getenv("IVITR_IDENT_X1", unset="1")),
                x2_coef=as.numeric(Sys.getenv("IVITR_IDENT_X2", unset="7")),
                defer_value=as.numeric(Sys.getenv("IVITR_IDENT_DEFER_VALUE",
                                                  unset="0")),
                n_rep=as.integer(Sys.getenv("IVITR_N_REP", unset="200")),
                n_eval=as.integer(Sys.getenv("IVITR_N_EVAL", unset="1000")),
                result_dir=Sys.getenv("IVITR_IDENT_DIR",
                                      unset="simulation_R2_ident")){
  configure_hpc_threads()
  library(doParallel)
  n_cores <- sim_cores()
  cl <- NULL
  if (n_cores > 1){
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    on.exit(stopCluster(cl), add=TRUE)
  } else {
    foreach::registerDoSEQ()
  }

  xi <- delta <- eta
  helper_exports <- c("first_stage_shift", "treatment_prob",
                      "first_stage_delta", "first_stage_bounds",
                      "draw_treatment", "outcome_linear_predictor",
                      "draw_outcome", "mean_outcome",
                      "cond_mean_outcome", "decision_from_effect",
                      "value_nonrejected", "value_population",
                      "wrong_treatment_prop", "effect_metric",
                      "true_cond_prob", "loss_components", "safe_proposed_decision",
                      "metric_rows")

  print(system.time(res <- foreach(i=1:n_rep, .packages="IVITR",
                                   .export=helper_exports,
                                   .errorhandling="stop",
                                   .combine=rbind)%dopar%{
    set.seed(i)

    x <- array(runif(n * p_covariates, min=-1, max=1),
               c(n, p_covariates))
    u <- runif(n, min=-1, max=1)
    z_tfd <- sample(c(-1L, 1L), n, replace=TRUE)
    trt <- draw_treatment(x, z_tfd, za_coef, x1_coef, x2_coef)
    y <- draw_outcome(x, u, trt, xi, delta, setting1, setting2)

    trt_tfd <- 2 * (trt - 0.5)
    y_tfd <- 2 * (y - 0.5)

    newx <- array(runif(n_eval * p_covariates, min=-1, max=1),
                  c(n_eval, p_covariates))
    newu <- runif(n_eval, min=-1, max=1)
    y1 <- mean_outcome(newx, newu, rep(1, n_eval), xi, delta, setting1, setting2)
    y0 <- mean_outcome(newx, newu, rep(0, n_eval), xi, delta, setting1, setting2)
    first_stage_x <- first_stage_delta(newx, za_coef, x1_coef, x2_coef)
    true_effect_x <- cond_mean_outcome(newx, 1, xi, delta, setting1, setting2) -
      cond_mean_outcome(newx, 0, xi, delta, setting1, setting2)

    # Use the population probability law for loss evaluation.  Estimating these
    # probabilities from the finite evaluation sample could make a nominal bound
    # fail to contain the known CATE merely through smoothing error.
    true_prob_bdd <- partial_bdd(
      true_cond_prob(newx, xi, delta, setting1, setting2,
                     za_coef, x1_coef, x2_coef),
      bdd_type="B-P"
    )
    lower_bdd <- apply(true_prob_bdd$p, 1, max)
    upper_bdd <- apply(true_prob_bdd$q, 1, min)
    if (any(true_effect_x < lower_bdd - 1e-10 |
            true_effect_x > upper_bdd + 1e-10)){
      stop("The population B-P bounds do not contain the true CATE.")
    }

    fit_result <- safe_proposed_decision(
      cf.iv_itr, y_tfd, trt_tfd, z_tfd, x, newx, d, debias=TRUE
    )
    fit0_result <- safe_proposed_decision(
      cf.iv_itr, y_tfd, trt_tfd, z_tfd, x, newx, d, debias=FALSE
    )
    fit_decision <- fit_result$decision
    fit0_decision <- fit0_result$decision

    fit_adhoc <- iv_itr.qiu(y_tfd, trt_tfd, z_tfd, x, newx=newx)
    adhoc_decision <- decision_from_effect(fit_adhoc$delta_eff)

    methods <- c("Augmented ERM", "ERM", "Ad-hoc")
    fit_loss <- loss_components(fit_decision, lower_bdd, upper_bdd,
                                true_effect_x, d)
    fit0_loss <- loss_components(fit0_decision, lower_bdd, upper_bdd,
                                 true_effect_x, d)
    adhoc_loss <- loss_components(adhoc_decision, lower_bdd, upper_bdd,
                                  true_effect_x, d)
    method_values <- list(
      population_value=setNames(c(
        value_population(y1, y0, fit_decision, defer_value),
        value_population(y1, y0, fit0_decision, defer_value),
        value_population(y1, y0, adhoc_decision, defer_value)
      ), methods),
      nondeferred_value=setNames(c(
        value_nonrejected(y1, y0, fit_decision),
        value_nonrejected(y1, y0, fit0_decision),
        value_nonrejected(y1, y0, adhoc_decision)
      ), methods),
      worst_case_loss=setNames(c(fit_loss["worst_case"],
                                 fit0_loss["worst_case"],
                                 adhoc_loss["worst_case"]), methods),
      true_cate_loss=setNames(c(fit_loss["true_cate"],
                                fit0_loss["true_cate"],
                                adhoc_loss["true_cate"]), methods),
      optimization_status=setNames(c(fit_result$status, fit0_result$status,
                                     "not_applicable"), methods),
      optimization_error=setNames(c(fit_result$error, fit0_result$error,
                                    NA_character_), methods),
      prop_deferred=setNames(c(
        mean(fit_decision == 0, na.rm=TRUE),
        mean(fit0_decision == 0, na.rm=TRUE),
        mean(adhoc_decision == 0, na.rm=TRUE)
      ), methods),
      prop_wrong=setNames(c(
        wrong_treatment_prop(fit_decision, true_effect_x),
        wrong_treatment_prop(fit0_decision, true_effect_x),
        wrong_treatment_prop(adhoc_decision, true_effect_x)
      ), methods),
      effect_bias=setNames(c(
        NA_real_,
        NA_real_,
        effect_metric(fit_adhoc$delta_eff, true_effect_x, "bias")
      ), methods),
      effect_mae=setNames(c(
        NA_real_,
        NA_real_,
        effect_metric(fit_adhoc$delta_eff, true_effect_x, "mae")
      ), methods),
      effect_mse=setNames(c(
        NA_real_,
        NA_real_,
        effect_metric(fit_adhoc$delta_eff, true_effect_x, "mse")
      ), methods)
    )

    metric_rows(method_values, d, i, n, eta, xi, delta, setting1, setting2,
                p_covariates, za_coef, x1_coef, x2_coef,
                mean(first_stage_x), defer_value)
  }))

  expected_rows <- n_rep * 30
  if (is.null(res) || nrow(res) != expected_rows){
    stop("Expected ", expected_rows, " result rows, but got ",
         if (is.null(res)) 0 else nrow(res), ".")
  }

  result_path <- resolve_result_dir(result_dir)
  dir.create(result_path, showWarnings=FALSE, recursive=TRUE)
  write.csv(
    res,
    file=file.path(result_path,
                   ident_result_name(n, eta, setting1, setting2, p_covariates,
                                     za_coef, x1_coef, x2_coef, d)),
    row.names=FALSE
  )
}

time_sim <- function(n=500, eta=0.5, setting1=1, setting2=1, d=0.5,
                     p_covariates=as.integer(Sys.getenv("IVITR_IDENT_P",
                                                        unset="5")),
                     za_coef=as.numeric(Sys.getenv("IVITR_IDENT_ZA",
                                                   unset="8")),
                     x1_coef=as.numeric(Sys.getenv("IVITR_IDENT_X1",
                                                   unset="1")),
                     x2_coef=as.numeric(Sys.getenv("IVITR_IDENT_X2",
                                                   unset="7")),
                     defer_value=as.numeric(Sys.getenv("IVITR_IDENT_DEFER_VALUE",
                                                       unset="0")),
                     result_dir=Sys.getenv("IVITR_IDENT_DIR",
                                           unset="simulation_R2_ident"),
                     timing_file=file.path(resolve_result_dir(result_dir),
                                           "runtime_log_ident.csv")){
  dir.create(dirname(timing_file), showWarnings=FALSE, recursive=TRUE)
  start_time <- Sys.time()
  setting_label <- paste0("identifiable n=", n, ", eta=", eta,
                          ", setting=(", setting1, ",", setting2, "), p=", p_covariates,
                          ", z_coef=", za_coef,
                          ", x1_coef=", x1_coef,
                          ", x2_coef=", x2_coef,
                          ", first_stage_range=(",
                          paste(round(first_stage_bounds(za_coef, x1_coef,
                                                         x2_coef), 3),
                                collapse=","),
                          ")",
                          ", u_a=0, d=", d)
  cat(sprintf("[%s] Starting identifiable sim: %s\n",
              format(start_time, "%Y-%m-%d %H:%M:%S %Z"), setting_label))
  flush.console()

  error_message <- NA_character_
  timing <- system.time(
    status <- tryCatch(
      {
        sim(n=n, eta=eta, setting1=setting1, setting2=setting2, d=d,
            p_covariates=p_covariates,
            za_coef=za_coef, x1_coef=x1_coef, x2_coef=x2_coef,
            defer_value=defer_value,
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

  cat(sprintf("[%s] Finished identifiable sim: %s | status=%s | elapsed=%.2f sec\n",
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
    n=n,
    eta=eta,
    lambda=0,
    xi=eta,
    delta=eta,
    setting1=setting1,
    setting2=setting2,
    d=d,
    p_covariates=p_covariates,
    za_coef=za_coef,
    x1_to_a_coef=x1_coef,
    x2_to_a_coef=x2_coef,
    first_stage_at_x0=plogis(za_coef) - plogis(-za_coef),
    first_stage_min=unname(first_stage_bounds(za_coef, x1_coef, x2_coef)["min"]),
    first_stage_max=unname(first_stage_bounds(za_coef, x1_coef, x2_coef)["max"]),
    u_a_multiplier=0,
    defer_value=defer_value
  )

  write.table(runtime, file=timing_file, sep=",", row.names=FALSE,
              col.names=!file.exists(timing_file), append=file.exists(timing_file))
  if (status != "success"){
    stop(error_message, call.=FALSE)
  }
  invisible(runtime)
}

run_driver <- function(args=commandArgs(trailingOnly=TRUE)){
  grid <- simulation_grid()
  task_id <- selected_task_id(args)

  cat("Identifiable simulation grid has", nrow(grid), "tasks.\n")
  if (!is.na(task_id)){
    if (task_id < 1 || task_id > nrow(grid)){
      stop("task-id must be between 1 and ", nrow(grid), ".")
    }
    grid <- grid[task_id,, drop=FALSE]
    cat("Running identifiable task", task_id, "only.\n")
  }
  cat("IVITR_NCORES=", sim_cores(), "\n", sep="")
  cat("IVITR_N_REP=", Sys.getenv("IVITR_N_REP", unset="200"), "\n", sep="")
  cat("IVITR_N_EVAL=", Sys.getenv("IVITR_N_EVAL", unset="1000"), "\n", sep="")
  cat("IVITR_IDENT_DIR=", Sys.getenv("IVITR_IDENT_DIR",
                                     unset="simulation_R2_ident"), "\n", sep="")
  cat("Identifiable design: manuscript outcome settings, no U-to-A term,",
      "first stage P(A=1|Z*,X)=expit(gamma_Z Z* + gamma_1 X1 - gamma_2 X2).\n")

  for (row_idx in seq_len(nrow(grid))){
    time_sim(
      n=grid$n[row_idx],
      eta=grid$eta[row_idx],
      setting1=grid$setting1[row_idx],
      setting2=grid$setting2[row_idx],
      d=grid$d[row_idx]
    )
  }
}

if (!interactive() && !identical(Sys.getenv("IVITR_SKIP_DRIVER"), "true")){
  run_driver()
}
