ignore_me <- setMethod("show", "Rcpp_RankedRSMatrix", function(object) {
  taille <- object$nrow() * object$ncol()
  if (taille == 0) return (cat("Null Ranked Row Sparse Matrix."))
  tx <- (taille - object$nvalid() )/taille
  
  cat("Ranked Row Sparse Matrix : ", object$nrow(), "rows x ", object$ncol(), "columns.\n")
  cat("Ratio of missing values : ", tx, "\n")
  
  n_elem <- min(10, object$nvalid())
  
  cat("common p : ", object$p[1:min(10, length(object$p))], "\n")
  cat("ranked j : ", object$jr[1:n_elem], "\n")
  cat("ranked x : ", object$xr[1:n_elem], "\n")
})