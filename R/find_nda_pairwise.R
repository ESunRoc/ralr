#' Pairwise reference identitication
#' 
#' A helper function to identify a set of reference taxa for use in the robust additive log-ratio (ralr) transformation.
#' 
#' @param M Matrix/data frame; relative abundance matrix of dimension `(n1+n2)`-by-`p` with subjects down the rows and features (taxa) across the columns.
#' @param n1 Integer; the number of observations at the first visit being evaluated.
#' @param n2 Integer; the number of observations at the second visit being evaluated.
#' @param alpha Numeric; the confidence bandwidth parameter taking values in the interval \eqn{[0,1]}. Defaults to 0.05.
#' @param n_theta Integer; the number of angles in the candidate grid. Defaults to 180, but can be adjusted depending on the number of taxa.
#' @param theta Numeric vector; the grid of angles parameterizing the confidence band. Defaults to `seq(0, atan(10), length = n_theta)`.
#' 
#' @return A list containing: `theta_selected`, the selected confidence band angle; and `id.nda`, the taxa identified for the reference set.
#' 
#' @export


find_nda_pairwise <- function(M, n1, n2, alpha = 0.05, n_theta = 180, theta = seq(0, atan(10), length = n_theta)){

  z <- qnorm(1-alpha/2)
  sin_theta <- sin(theta); cos_theta <- cos(theta)
  
  Xbar <- colMeans(M[1:n1,]); Ybar <- colMeans(M[(n1+1):(n1+n2),])
  
  which_comps_within <- matrix(ncol = ncol(M), nrow = n_theta)
  
  for(k in 1:n_theta){
    D_k <- atan(Ybar/Xbar) - theta[k]
    se_Dk_num <- sqrt(n1*Xbar * n2^2*Ybar^2 * (1-Xbar) + n1^2*Xbar^2 * n2*Ybar * (1-Ybar))
    se_Dk_denom <- (n1^2*Xbar^2 + n2^2*Ybar^2)
    
    which_comps_within[k,] <- (abs(D_k) <= z*(se_Dk_num/se_Dk_denom))
  }
  
  theta_idmax <- which.max(rowSums(which_comps_within))
  theta_max <- theta[theta_idmax]
  id_for_ref <- which(which_comps_within[theta_idmax,])
  
  return(list("theta_selected" = theta_max, "id.nda" = id_for_ref))
}
