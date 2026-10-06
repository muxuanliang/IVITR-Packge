cvxr_value <- function(result, variable){
  status <- tryCatch(result$status, error = function(e) NA_character_)
  if (length(status) > 0 && !is.na(status) &&
      !status %in% c("optimal", "optimal_inaccurate")){
    stop("CVXR solver returned status ", status, ".")
  }

  if (!is.null(result$getValue) && is.function(result$getValue)){
    return(result$getValue(variable))
  }

  if ("value" %in% getNamespaceExports("CVXR")){
    return(CVXR::value(variable))
  }

  stop("Could not extract CVXR solution value. ",
       "The CVXR result has no callable getValue() method and ",
       "CVXR::value() is not available.")
}
