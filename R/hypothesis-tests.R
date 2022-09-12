#' Wald test for inputs
#'
#'
#' @param X Data
#' @param y Response
#' @param W Weight vector
#' @param q Number of hidden nodes
#' @return Wald hypothesis test for each input
#' @export
wald_test <- function(X, y, W, q) {

  p <- ncol(X)
  n <- nrow(X)

  nn <- nnet::nnet(X, y, size = q, linout = T, Hess = T, maxit = 0, trace = F,
             Wts = W)

  sigma2 <- nn$value / n  # estimate \sigma^2

  Sigma_inv <- nn$Hessian / (2*sigma2) # $\Sigma^-1 = I(\theta)$

  Sigma <- solve(Sigma_inv)

  p_values <- rep(NA, p)
  chisq <- rep(NA, p)

  for (i in 1:p) {

    # stores which weights correspond input unit i
    ind_vec <- sapply(X = 1:q,
                      FUN = function(x) (x - 1)*(p + 1) + 1 + i)

    theta_x <- W[ind_vec]
    Sigma_inv_x <- solve(Sigma[ind_vec, ind_vec])

    chisq[i] <- t(theta_x) %*% Sigma_inv_x %*% theta_x

    p_values[i] <- 1 - pchisq(chisq[i], df = q)
  }

  return(list("chisq" = chisq, "p_value" = p_values))

}
