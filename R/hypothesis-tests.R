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


#' Likelihood ratio test for inputs
#'
#'
#' @param X Data
#' @param y Response
#' @param W Weight vector
#' @param q Number of hidden nodes
#' @param n_init Number of initialisations
#' @param unif Value for generating random weights
#' @param maxit Maximum number of iterations for nnet
#' @return Wald hypothesis test for each input
#' @export
lr_test <- function(X, y, W, q, n_init = 1, unif = 3, maxit = 1000, ...) {

  n <- nrow(X)
  p <- ncol(X)

  y_hat_full <- nn_pred(X, W, q)

  k_full <- (p + 2) * q + 1

  RSS_full <- sum((y - y_hat_full) ^ 2)

  sigma2 <- RSS_full / n

  log_like_full = (- n / 2) * log(2 * pi * sigma2) - (RSS / (2 * sigma2))

  k_rem <- (p + 1) * q + 1

  deg_freedom <- k_full - k_rem

  p_values <- rep(NA, p)
  chisq <- rep(NA, p)

  for (i in 1:p) {

    X_temp <- X[, -i]

    weight_matrix_init <- matrix(
      runif(k_rem * n_init, min = - unif, max = unif),
      nrow = n_init,
      byrow = T)

    log_like <- rep(NA, n_init)

    for (j in 1:n_init) {

      nn_model <-  nnet::nnet(X_temp, y, size = q, trace = FALSE, linout = TRUE,
                              Wts = weight_matrix_init[j, ], maxit = maxit,
                              ...)

      log_like[j] <- nn_loglike(nn_model)

    }

    log_like_rem <- max(log_like)

    chisq[i] <- - 2 * (log_like_rem - log_like_full)

    p_values[i] = pchisq(chisq[i], df = deg_freedom, lower.tail = F)
  }
  return(list("chisq" = chisq, "p_value" = p_values))
}