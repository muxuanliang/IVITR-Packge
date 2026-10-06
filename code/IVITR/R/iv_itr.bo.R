iv_itr.bo <- function(y, trt, z, x, split_index = NULL, lambda=0.01, bandwidth=0.1, bdd_type=c("B-P", "Siddique")){
  sample_size <- length(y)
  if (is.null(split_index)){
    split_index <- rnorm(sample_size)>0
  }
  data_ns <- list(y=y[split_index], trt=trt[split_index], z=z[split_index], x=x[split_index,])
  data_solve <- list(y=y[!split_index], trt=trt[!split_index], z=z[!split_index], x=x[!split_index,])

  cond_prob <- get_cond_prob(y=data_ns$y, trt=data_ns$trt, z=data_ns$z, x=data_ns$x, newx=data_solve$x)
  pq <- partial_bdd(cond_prob_array = cond_prob, bdd_type = bdd_type)
  k_matrix <- rbf_gauss(data_solve$x , bandwidth = bandwidth)
  p_max <- apply(pq$p,1,max)
  q_min <- apply(pq$q,1,min)

  # pseudo treatment and weights
  pseudo_x <- data_solve$x
  pseudo_trt <- -((q_min < 0) - (p_max > 0) - sign(abs(q_min)-abs(p_max)) * (q_min >= 0) * (p_max <= 0))
  pseudo_weight <- abs(q_min) * (p_max > 0) + abs(p_max) * (q_min < 0) + abs(abs(q_min)-abs(p_max)) * (q_min >= 0) * (p_max <= 0)

  # warm-up
  fit_init <- iv_itr_bo(weight=pseudo_weight, trt=pseudo_trt, x=pseudo_x, lambda=lambda, bandwidth=bandwidth)
  alpha_updated <- fit_init$alpha_est

  list(alpha_est=alpha_updated * (abs(alpha_updated)>1e-5), basis_x=data_solve$x,
       status="optimal", alpha_warm=fit_init$alpha_est, bandwidth=bandwidth,
       bdd_type=bdd_type)
}
