defer_decision_from_link <- function(link, d = 0.5){
  if (!is.numeric(link)){
    stop("link must be numeric.")
  }
  if (!is.numeric(d) || length(d) != 1 || is.na(d) || d <= 0 || d > 0.5){
    stop("d must be a single number in (0, 0.5] for link-threshold decisions.")
  }

  cutoff <- 0.5 - d
  as.integer((link > cutoff) - (link < -cutoff))
}

oracle_score_decision <- function(score, d = 0.5){
  defer_decision_from_link(score, d = d)
}

oracle_bdd_decision <- function(cond_prob_array = NULL, p = NULL, q = NULL,
                                d = 0.5, bdd_type = "B-P"){
  if (!is.numeric(d) || length(d) != 1 || is.na(d) || d <= 0 || d > 1){
    stop("d must be a single number in (0, 1] for bound-based oracle decisions.")
  }
  if (is.null(p) || is.null(q)){
    if (is.null(cond_prob_array)){
      stop("Provide either cond_prob_array or both p and q.")
    }
    pq <- partial_bdd(cond_prob_array = cond_prob_array, bdd_type = bdd_type)
    p <- pq$p
    q <- pq$q
  }

  p_max <- if (is.null(dim(p))) as.numeric(p) else apply(p, 1, max)
  q_min <- if (is.null(dim(q))) as.numeric(q) else apply(q, 1, min)
  if (length(p_max) != length(q_min)){
    stop("p and q must have the same number of observations.")
  }

  l_minus <- abs(p_max) * (p_max <= 0)
  u_plus <- abs(q_min) * (q_min >= 0)

  decision <- rep(NA_integer_, length(p_max))
  observed <- !(is.na(l_minus) | is.na(u_plus))
  decision[observed] <- 0L
  decision[observed & l_minus < d * u_plus] <- 1L
  decision[observed & u_plus < d * l_minus] <- -1L
  decision
}

normalize_defer_decision <- function(decision){
  if (!is.numeric(decision)){
    stop("decision must be numeric with values -1, 0, or 1.")
  }
  valid <- is.na(decision) | decision %in% c(-1, 0, 1)
  if (!all(valid)){
    stop("decision must only contain -1, 0, 1, or NA.")
  }
  as.integer(decision)
}

bdd_bounds <- function(cond_prob_array = NULL, p = NULL, q = NULL,
                       bdd_type = "B-P"){
  if (is.null(p) || is.null(q)){
    if (is.null(cond_prob_array)){
      stop("Provide either cond_prob_array or both p and q.")
    }
    pq <- partial_bdd(cond_prob_array = cond_prob_array, bdd_type = bdd_type)
    p <- pq$p
    q <- pq$q
  }

  lower <- if (is.null(dim(p))) as.numeric(p) else apply(p, 1, max)
  upper <- if (is.null(dim(q))) as.numeric(q) else apply(q, 1, min)
  if (length(lower) != length(upper)){
    stop("p and q must have the same number of observations.")
  }

  list(lower=lower, upper=upper)
}

defer_review_proportions <- function(estimated, y1 = NULL, y0 = NULL,
                                     true_effect = NULL,
                                     cond_prob_array = NULL,
                                     p = NULL, q = NULL,
                                     d = 0.5,
                                     estimated_type = c("link", "decision"),
                                     bdd_type = "B-P"){
  estimated_type <- match.arg(estimated_type)
  estimated_decision <- if (estimated_type == "link"){
    defer_decision_from_link(estimated, d = d)
  } else {
    normalize_defer_decision(estimated)
  }

  if (is.null(true_effect)){
    stop("Provide true_effect as E[Y(1)|X] - E[Y(0)|X].")
  }

  bounds <- bdd_bounds(cond_prob_array = cond_prob_array, p = p, q = q,
                       bdd_type = bdd_type)

  if (length(estimated_decision) != length(true_effect) ||
      length(estimated_decision) != length(bounds$lower)){
    stop("estimated, true_effect, and bounds must have the same length.")
  }

  observed <- !(is.na(estimated_decision) | is.na(true_effect) |
                  is.na(bounds$lower) | is.na(bounds$upper))
  estimated_decision <- estimated_decision[observed]
  true_effect <- true_effect[observed]
  lower <- bounds$lower[observed]
  upper <- bounds$upper[observed]

  deferred <- estimated_decision == 0
  nondeferred <- estimated_decision != 0
  bounds_cross_zero <- lower <= 0 & upper >= 0
  wrong_treatment <- estimated_decision * true_effect < 0

  safe_prop <- function(event, group){
    if (!any(group)){
      return(NA_real_)
    }
    mean(event[group])
  }

  c(
    prop_bounds_cross_zero_among_deferred = safe_prop(bounds_cross_zero, deferred),
    prop_wrong_treatment_among_nondeferred = safe_prop(wrong_treatment, nondeferred),
    n_deferred = sum(deferred),
    n_nondeferred = sum(nondeferred),
    n_complete = length(estimated_decision)
  )
}

defer_assignment_proportions <- function(estimated, oracle = NULL,
                                         y1 = NULL, y0 = NULL,
                                         cond_prob_array = NULL,
                                         p = NULL, q = NULL,
                                         d = 0.5,
                                         estimated_type = c("link", "decision"),
                                         oracle_type = c("bdd", "score", "decision"),
                                         bdd_type = "B-P"){
  estimated_type <- match.arg(estimated_type)
  oracle_type <- match.arg(oracle_type)

  estimated_decision <- if (estimated_type == "link"){
    defer_decision_from_link(estimated, d = d)
  } else {
    normalize_defer_decision(estimated)
  }

  oracle_decision <- if (oracle_type == "bdd"){
    oracle_bdd_decision(cond_prob_array = cond_prob_array, p = p, q = q,
                        d = d, bdd_type = bdd_type)
  } else if (oracle_type == "score"){
    if (is.null(oracle)){
      if (is.null(y1) || is.null(y0)){
        stop("Provide oracle score or both y1 and y0 for score-based oracle decisions.")
      }
      oracle <- y1 - y0
    }
    oracle_score_decision(oracle, d = d)
  } else {
    if (is.null(oracle)){
      stop("Provide oracle decisions when oracle_type = 'decision'.")
    }
    normalize_defer_decision(oracle)
  }

  if (length(estimated_decision) != length(oracle_decision)){
    stop("estimated and oracle decisions must have the same length.")
  }

  observed <- !(is.na(estimated_decision) | is.na(oracle_decision))
  estimated_decision <- estimated_decision[observed]
  oracle_decision <- oracle_decision[observed]

  oracle_defer <- oracle_decision == 0
  oracle_active <- oracle_decision != 0

  safe_prop <- function(event, group){
    if (!any(group)){
      return(NA_real_)
    }
    mean(event[group])
  }

  c(
    prop_active_given_oracle_defer = safe_prop(estimated_decision != 0, oracle_defer),
    prop_defer_given_oracle_active = safe_prop(estimated_decision == 0, oracle_active),
    n_oracle_defer = sum(oracle_defer),
    n_oracle_active = sum(oracle_active),
    n_complete = length(oracle_decision)
  )
}
