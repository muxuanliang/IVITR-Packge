value_bdd <- function(fit,y, trt, z, x, split_index = NULL, d=0.5){
  sample_size <- length(y)
  if (is.null(split_index)){
    split_index <- TRUE
  }
  data_evaluate <- list(y=y[split_index], trt=trt[split_index], z=z[split_index], x=x[split_index,])
  k_matrix <- rbf_gauss(x=fit$basis_x, newx=data_evaluate$x, bandwidth = fit$bandwidth)
  predict.iv_itr <- k_matrix %*% fit$alpha_est
  predict.decision <- (predict.iv_itr>0.5)-(predict.iv_itr< -0.5)
  
  # probability bound
  cond_prob <- get_cond_prob(y=data_evaluate$y, trt=data_evaluate$trt, 
                             z=data_evaluate$z, x=data_evaluate$x, newx=data_evaluate$x)
  pq <- partial_bdd(cond_prob_array = cond_prob, bdd_type = fit$bdd_type)
  p_max <- apply(pq$p,1,max)
  q_min <- apply(pq$q,1,min)
  
  # objective
  obj <- mean(d*(q_min-p_max)*(q_min>0 & p_max<0)*(predict.decision==0)+
    (predict.decision==1)*(p_max<0) * abs(p_max)+
    (predict.decision==-1)*(q_min>0) * abs(q_min))
  
  obj
}