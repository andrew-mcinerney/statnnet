#' Network Output Function
#'
#'
#' @param X Data
#' @param W Weight vector
#' @param q Number of hidden units
#' @return Prediction for given inputs
#' @export
nn_pred <- function(X, W, q) {
  n <- nrow(X)
  p <- ncol(X)

  if(length(W) == ((p + 2)*q + 1)){
    X <- cbind(rep(1, n), X)

    h_input <- X %*% t(matrix(W[1:((p + 1)*q)], nrow = q, byrow = TRUE))

    h_act <- cbind(rep(1, n), sigmoid(h_input))

    y_hat <- h_act %*% matrix(W[c((length(W) - q):length(W))], ncol = 1)

    return(y_hat)

  }else{
    return(print('Error: Incorrect number of weights for NN structure'))
  }
}

#' Sigmoid activation function
#'
#'
#' @param x Input
#' @return Sigmoid function
#' @export
sigmoid <- function(x) 1 / (1 + exp(- x))

#' Neural Network Normal Log-likelihood Value
#'
#'
#' @param X Data
#' @param W Weight vector
#' @param q Number of hidden units
#' @return Prediction for given inputs
#' @export
nn_loglike <- function(nn_model){

  n <- nrow(nn_model$residuals)

  RSS <- nn_model$value

  sigma2 <- RSS / n

  log_like <- (- n / 2) * log(2 * pi * sigma2) - RSS / (2 * sigma2)

  return(log_like)

}

#' Difference in average prediction for values above and below median
#'
#'
#' @param X Data
#' @param W Weight vector
#' @param q Number of hidden units
#' @return Effect for each input
#' @export
covariate_eff <- function(X, W, q) {
  eff <- rep(NA, ncol(X))
  for (col in 1:ncol(X)) {
    low <- X[X[, col] <= median(X[, col]), ]
    high <- X[X[, col] > median(X[, col]), ]

    eff[col] <- mean(nn_pred(high, W, q)) - mean(nn_pred(low, W, q))
  }
  names(eff) <- colnames(X)
  return(eff)
}

#' Partial Dependence Plot for one std. dev. increase
#'
#'
#' @param W Weight vector
#' @param X Data
#' @param q Number of hidden units
#' @param ind index of column to plot
#' @return Effect for each input
#' @export
pdp_effect <- function(W, X, q, ind, x_r = c(-3, 3), len = 301){
  sd_m <- matrix(0, ncol = ncol(X), nrow = nrow(X))
  sd_m[, ind] <- sd(X[, ind])

  x <- seq(from = x_r[1], to = x_r[2], length.out = len)

  eff <- rep(NA, len)

  for (i in 1:len){
    X[, ind] <- x[i]
    eff[i] <- mean(nn_pred(X + sd_m, W, q) - nn_pred(X, W, q))
  }
  return(eff)
}
