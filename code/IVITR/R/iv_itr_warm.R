iv_itr_warm <- function(y, trt, z, x, split_index = NULL, lambda=0.01, d=0.5, bandwidth=0.1, 
                        bdd_type=c("B-P", "Siddique"), solver="CVXR", max_iter=1000){
  sample_size <- length(y)
  if (is.null(split_index)){
    split_index <- rnorm(sample_size)>0
  }
  data_ns <- list(y=y[split_index], trt=trt[split_index], z=z[split_index], x=x[split_index,])
  data_solve <- list(y=y[!split_index], trt=trt[!split_index], z=z[!split_index], x=x[!split_index,])

  cond_prob <- get_cond_prob(y=data_ns$y, trt=data_ns$trt, z=data_ns$z, x=data_ns$x, newx=data_solve$x)
  pq <- partial_bdd(cond_prob_array = cond_prob, bdd_type = bdd_type)
  k_matrix <- rbf_gauss(data_solve$x, bandwidth = bandwidth)
  p_max <- apply(pq$p,1,max)
  q_min <- apply(pq$q,1,min)
  positive_weights <- abs(q_min)*(q_min>=0)
  negative_weights <- abs(p_max)*(p_max<=0)
  
  if (solver=="CVXR"){
    # square-root of K
    svd_k <- svd(k_matrix)
    chol_K <- diag(sqrt(svd_k$d)+0.01) %*% t(svd_k$u)
    
    # use convex optimization
    alpha <- CVXR::Variable(sum(!split_index))
    slack_positive <- CVXR::Variable(sum(!split_index))
    slack_negative <- CVXR::Variable(sum(!split_index))
    objective <- CVXR::Minimize(mean(positive_weights*slack_positive)+mean(negative_weights*slack_negative)+lambda* sum((chol_K %*% alpha)^2))
    constraint1 <- slack_positive >= 0
    constraint2 <- slack_positive >= 1-k_matrix %*% alpha
    constraint3 <- slack_positive >= 1-(1-d)/d * k_matrix %*% alpha
    constraint4 <- slack_negative >= 0
    constraint5 <- slack_negative >= 1+k_matrix %*% alpha
    constraint6 <- slack_negative >= 1+(1-d)/d * k_matrix %*% alpha
    problem <- CVXR::Problem(objective, constraints = list(constraint1, constraint2,constraint3, constraint4,
                                                           constraint5, constraint6))
    result <- CVXR::solve(problem, solver = "OSQP")
    alpha_value <- cvxr_value(result, alpha)
    list(alpha_est=alpha_value * (abs(alpha_value)>1e-5), basis_x=data_solve$x)
  } else if (solver=="HD"){
    # square-root of K
      alpha_positive_prev1 <- alpha_negative_prev1 <- rep(0, sum(!split_index))
      alpha_positive_prev2 <- alpha_negative_prev2 <- rep(0, sum(!split_index))
      alpha_prev <- alpha_positive_prev1+(1-d)/d *alpha_positive_prev2 -alpha_negative_prev1-(1-d)/d *alpha_negative_prev2
      iter <- 0
      diff_value <- 1
      batch <- floor(sum(!split_index)/1000)+1
      sample_index <- sample(1:batch, sum(!split_index), replace = TRUE)
      while(iter <= max_iter & abs(diff_value)>1e-5){
        update_index <- iter%%batch+1
        svd_k <- svd(k_matrix[sample_index==update_index, sample_index==update_index])
        chol_K <- diag(sqrt(svd_k$d)+0.01) %*% t(svd_k$u)
        # use convex optimization
        alpha_positive1 <- CVXR::Variable(sum(sample_index==update_index))
        alpha_negative1 <- CVXR::Variable(sum(sample_index==update_index))
        alpha_positive2 <- CVXR::Variable(sum(sample_index==update_index))
        alpha_negative2 <- CVXR::Variable(sum(sample_index==update_index))
        alpha <- CVXR::Variable(sum(sample_index==update_index))
        objective <- CVXR::Minimize(0.5* sum((chol_K %*% alpha)^2)
                                      +t(alpha_positive_prev1[sample_index!=update_index]+(1-d)/d *alpha_positive_prev2[sample_index!=update_index] -alpha_negative_prev1[sample_index!=update_index]-(1-d)/d *alpha_negative_prev2[sample_index!=update_index]) %*% 
                                      k_matrix[sample_index!=update_index, sample_index==update_index]%*% alpha-
                                      sum(alpha_positive1+alpha_negative1+alpha_positive2+alpha_negative2))
        constraint0 <- alpha==alpha_positive1+(1-d)/d*alpha_positive2-alpha_negative1-(1-d)/d *alpha_negative2
        constraint1 <- alpha_positive1 >= 0
        constraint2 <- alpha_positive2 >= 0
        constraint3 <- (alpha_positive1+alpha_positive2) <= positive_weights[sample_index==update_index]/(2*lambda*sum(!split_index))
        constraint4 <- alpha_negative1 >= 0
        constraint5 <- alpha_negative2 >= 0
        constraint6 <- (alpha_negative1+alpha_negative2) <= negative_weights[sample_index==update_index]/(2*lambda*sum(!split_index))
        problem <- CVXR::Problem(objective, constraints = list(constraint0,constraint1, constraint2,constraint3, constraint4,constraint5, constraint6))
        result <- CVXR::solve(problem, solver = "OSQP")
        alpha_positive_updated1 <- alpha_positive_prev1
        alpha_negative_updated1 <- alpha_negative_prev1
        alpha_positive_updated2 <- alpha_positive_prev2
        alpha_negative_updated2 <- alpha_negative_prev2
        alpha_positive_updated1[sample_index==update_index] <- cvxr_value(result, alpha_positive1)
        alpha_negative_updated1[sample_index==update_index] <- cvxr_value(result, alpha_negative1)
        alpha_positive_updated2[sample_index==update_index] <- cvxr_value(result, alpha_positive2)
        alpha_negative_updated2[sample_index==update_index] <- cvxr_value(result, alpha_negative2)
        alpha_updated <- alpha_positive_updated1+(1-d)/d *alpha_positive_updated2 -alpha_negative_updated1-(1-d)/d *alpha_negative_updated2
        diff_value <- mean(positive_weights*pmax(0,(1-k_matrix %*% alpha_updated ),(1-(1-d)/d * k_matrix %*% alpha_updated)))+
          mean(negative_weights*pmax(0,(1+k_matrix %*% alpha_updated),(1+(1-d)/d * k_matrix %*% alpha_updated)))+
          lambda* t(alpha_updated) %*% k_matrix %*% alpha_updated-mean(positive_weights*pmax(0,(1-k_matrix %*% alpha_prev),(1-(1-d)/d * k_matrix %*% alpha_prev)))-
          mean(negative_weights*pmax(0,(1+k_matrix %*% alpha_prev),(1+(1-d)/d * k_matrix %*% alpha_prev)))-
          lambda* t(alpha_prev) %*% k_matrix %*% alpha_prev
        alpha_prev <- alpha_updated
        alpha_positive_prev1 <- alpha_positive_updated1
        alpha_positive_prev2 <- alpha_positive_updated2
        alpha_negative_prev1 <- alpha_negative_updated1
        alpha_negative_prev2 <- alpha_negative_updated2
        iter <- iter+1
      }
    list(alpha_est=alpha_updated * (abs(alpha_updated)>1e-5), basis_x=data_solve$x, alpha_positive_prev1=alpha_positive_prev1,
         alpha_positive_prev2=alpha_positive_prev2,alpha_negative_prev1=alpha_negative_prev1, alpha_negative_prev2=alpha_negative_prev2)
  }
}
