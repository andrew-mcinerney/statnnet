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
#' @param x_r x-axis range
#' @param len number of breaks for x-axis
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

#' Perform mle simulation for a function FUN to calculate associated uncertainty
#'
#'
#' @param W Weight vector
#' @param X Data
#' @param y Response
#' @param q Number of hidden units
#' @param ind index of column to plot
#' @param B number of replicates
#' @param alpha significance level
#' @param x_r x-axis range
#' @param len number of breaks for x-axis
#' @return Effect for each input
#' @export
mlesim <- function (W, X, y, q, ind, FUN, B = 1000, alpha = 0.05, x_r = c(-3, 3),
                    len = 301) {

  nn <- nnet::nnet(y~., data = data.frame(X, y), size = q, Wts = W,
                   linout = TRUE, trace = FALSE, maxit = 0, Hess = TRUE)

  sigma2 <- nn$value/n # nn$value = RSS

  Sigma_inv <- nn$Hessian/(2*sigma2)

  Sigma_hat <- solve(Sigma_inv)

  sim <- MASS::mvrnorm(n = B, mu = W, Sigma = Sigma_hat)

  pred <- apply(sim, 1, function(x) FUN(x, X = X, ind = ind, q = q,
                                        x_r = x_r, len = len))


  lower <- apply(pred, 1, quantile, probs = alpha / 2)
  upper <- apply(pred, 1, quantile, probs = 1 - alpha / 2)

  return(list('upper' = upper, 'lower' = lower))
}

#' Perform delta method for a function FUN to calculate associated uncertainty
#'
#'
#' @param W Weight vector
#' @param X Data
#' @param y Response
#' @param q Number of hidden units
#' @param ind index of column to plot
#' @param alpha significance level
#' @param x_r x-axis range
#' @param len number of breaks for x-axis
#' @return Effect for each input
#' @export
delta_method <- function(W, X, y, q, ind, FUN, alpha = 0.05, x_r = c(-3, 3),
                         len = 301, ...){

  nn <- nnet::nnet(X, y, size = q, Wts = W, linout = TRUE, Hess = TRUE,
                   maxit = 0, trace = FALSE)

  sigma2 <- nn$value / nrow(X)  # estimate \sigma^2

  Sigma_inv <- nn$Hessian/(2*sigma2)

  Sigma_hat <- solve(Sigma_inv)

  gradient <- numDeriv::jacobian(func = FUN,
                                 x = W,
                                 X = X,
                                 ind = ind,
                                 q = q,
                                 x_r = x_r,
                                 len = len,
                                 ...)

  var_est <- as.matrix(gradient) %*% Sigma_hat %*% t(as.matrix(gradient))

  pred <- FUN(W = W, X = X, q = q, ind = ind, x_r = x_r, len = len, ...)

  upper <- pred + qnorm(1 - alpha / 2) * sqrt(diag(var_est))
  lower <- pred + qnorm(alpha / 2) * sqrt(diag(var_est))

  return(list('upper' = upper, 'lower' = lower))
}
