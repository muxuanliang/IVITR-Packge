iv_itr.qiu <- function(y, trt, z, x, newx){
  # estimate delta_y and delta_trt
  propensity <- ks(x, (z+1)/2, x)
  delta_y <- ks(x, z/propensity * y, newx)
  delta_trt <- ks(x, z/propensity * trt, newx)

  # estimate effect
  delta_eff <- delta_y/delta_trt

  list(delta_eff=delta_eff)
}
