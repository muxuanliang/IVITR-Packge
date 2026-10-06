iv_itr_bo <- function(weights, trt, x, lambda=0.01, bandwidth=0.1){
  k_matrix <- rbf_gauss(x, bandwidth = bandwidth)
  chol_K <- chol(k_matrix)

  # use convex optimization
  alpha <- CVXR::Variable(length(trt))
  slack <- CVXR::Variable(length(trt))
  objective <- CVXR::Minimize(mean(weights*slack)+lambda* sum((chol_K %*% alpha)^2))
  constraint1 <- slack >= 0
  constraint2 <- slack >= 1-k_matrix %*% alpha
  problem <- CVXR::Problem(objective, constraints = list(constraint1, constraint2))
  result <- CVXR::solve(problem, solver = "OSQP")

  alpha_value <- cvxr_value(result, alpha)
  list(alpha_est=alpha_value * (abs(alpha_value)>1e-5), basis_x=x)
}
