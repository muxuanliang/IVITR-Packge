value_bdd.qiu <- function(delta_est,cutoff,y, trt, z, x, d=0.5, bdd_type="B-P"){
  sample_size <- length(y)
  data_evaluate <- list(y=y, trt=trt, z=z, x=x)
  predict.decision <- (delta_est>cutoff)-(delta_est< -cutoff)

  # probability bound
  cond_prob <- get_cond_prob(y=data_evaluate$y, trt=data_evaluate$trt,
                             z=data_evaluate$z, x=data_evaluate$x, newx=data_evaluate$x)
  pq <- partial_bdd(cond_prob_array = cond_prob, bdd_type = bdd_type)
  p_max <- apply(pq$p,1,max)
  q_min <- apply(pq$q,1,min)

  # objective
  obj <- mean(d*(q_min-p_max)*(q_min>0 & p_max<0)*(predict.decision==0)+
                (predict.decision==1)*(p_max<0) * abs(p_max)+
                (predict.decision==-1)*(q_min>0) * abs(q_min))

  obj
}
