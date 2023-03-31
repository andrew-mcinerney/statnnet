#' Neural network prediction
#'
#'
#' @param X Data
#' @param W Weight vector
#' @param q Number of hidden nodes
#' @param response Response type: `"continuous"` (default) or
#'  `"binary"`
#' @return Prediction for given inputs
#' @export
nn_pred <- function(X, W, q, response = "continuous") {
  n <- nrow(X)
  p <- ncol(X)

  k <- (p + 2) * q + 1

  if (length(W) == k) {
    X <- cbind(rep(1, n), X)

    h_input <- as.matrix(X) %*% t(matrix(W[1:((p + 1) * q)], nrow = q, byrow = TRUE))

    h_act <- cbind(rep(1, n), sigmoid(h_input))

    if (response == "continuous") {
      y_hat <- h_act %*% matrix(W[c((length(W) - q):length(W))], ncol = 1)
    } else if (response == "binary") {
      y_hat <- sigmoid(
        h_act %*% matrix(W[c((length(W) - q):length(W))], ncol = 1)
      )
    } else {
      stop(
        sprintf(
          "Error: %s not recognised as available response type.",
          output
        )
      )
    }

    return(y_hat)
  } else {
    stop(sprintf(
      "Error: Incorrect number of weights for NN structure. W should have
      %s weights (%s weights supplied).", k, length(W)
    ))
  }
}

#' Sigmoid activation function
#'
#'
#' @param x Input
#' @return Sigmoid function
#' @export
sigmoid <- function(x) 1 / (1 + exp(-x))

#' Neural Network Normal Log-likelihood Value
#'
#'
#' @param nn_model nnet object
#' @return Log-Likelihhod value
#' @export
nn_loglike <- function(nn_model) {
  n <- nrow(nn_model$residuals)

  RSS <- nn_model$value

  sigma2 <- RSS / n

  log_like <- (-n / 2) * log(2 * pi * sigma2) - RSS / (2 * sigma2)

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
    low <- X[X[, col] <= stats::median(X[, col]), ]
    high <- X[X[, col] > stats::median(X[, col]), ]

    eff[col] <- mean(nn_pred(high, W, q)) - mean(nn_pred(low, W, q))
  }
  names(eff) <- colnames(X)
  return(eff)
}



#' Perform m.l.e. simulation for a function FUN to calculate associated uncertainty
#'
#'
#' @param W Weight vector
#' @param X Data
#' @param y Response
#' @param q Number of hidden units
#' @param ind index of column to plot
#' @param FUN function for m.l.e. simulation
#' @param B number of replicates
#' @param alpha significance level
#' @param x_r x-axis range
#' @param len number of breaks for x-axis
#' @param lambda Ridge penalty. Default is 0.
#' @param response Response type: `"continuous"` (default) or
#'  `"binary"`
#' @param ... additional arguments to FUN
#' @return Effect for each input
#' @export
mlesim <- function(W, X, y, q, ind, FUN, B = 1000, alpha = 0.05, x_r = c(-3, 3),
                   len = 301, lambda = 0, response = "continuous", ...) {
  vc <- VC(W, X, y, q, lambda = lambda, response = response)

  sim <- MASS::mvrnorm(n = B, mu = W, Sigma = vc)

  pred <- apply(sim, 1, function(x) {
    FUN(x,
      X = X, q = q,
      x_r = x_r, len = len, ...
    )
  })


  lower <- apply(pred, 1, stats::quantile, probs = alpha / 2)
  upper <- apply(pred, 1, stats::quantile, probs = 1 - alpha / 2)

  return(list("upper" = upper, "lower" = lower))
}

#' Perform delta method for a function FUN to calculate associated uncertainty
#'
#'
#' @param W Weight vector
#' @param X Data
#' @param y Response
#' @param q Number of hidden units
#' @param ind index of column to plot
#' @param FUN function for delta method
#' @param alpha significance level
#' @param x_r x-axis range
#' @param len number of breaks for x-axis
#' @param lambda Ridge penalty. Default is 0.
#' @param response Response type: `"continuous"` (default) or
#'  `"binary"`
#' @param ... additional arguments to FUN
#' @return Effect for each input
#' @export
delta_method <- function(W, X, y, q, ind, FUN, alpha = 0.05, x_r = c(-3, 3),
                         len = 301, lambda = 0, response = "continuous", ...) {

  vc <- VC(W, X, y, q, lambda = lambda, response = response)

  gradient <- numDeriv::jacobian(
    func = FUN,
    x = W,
    X = X,
    q = q,
    x_r = x_r,
    len = len,
    ...
  )

  var_est <- as.matrix(gradient) %*% vc %*% t(as.matrix(gradient))

  pred <- FUN(W = W, X = X, q = q, x_r = x_r, len = len, ...)

  upper <- pred + stats::qnorm(1 - alpha / 2) * sqrt(diag(var_est))
  lower <- pred + stats::qnorm(alpha / 2) * sqrt(diag(var_est))

  return(list("upper" = upper, "lower" = lower))
}


#' Calculate variance-covariance matrix for nnet object
#'
#'
#' @param W Weight vector
#' @param X Data
#' @param y Response
#' @param q Number of hidden units
#' @param lambda Ridge penalty. Default is 0.
#' @param response Response type: `"continuous"` (default) or
#'  `"binary"`
#' @param ... additional arguments to nnet
#' @return Effect for each input
#' @export
VC <- function(W, X, y, q, lambda = 0, response = "continuous") {

  if (response == "continuous") {
    linout <- TRUE
    entropy <- FALSE
  } else if (task == "binary") {
    linout <- FALSE
    entropy <- TRUE
  } else {
    stop(sprintf(
      "Error: %s not recognised as response. Please choose continuous or binary",
      response
    ))
  }

  nn_0 <- nnet::nnet(x = X, y = y, size = q, Wts = W, decay = lambda, linout = linout,
                     entropy = entropy, maxit = 0, Hess = TRUE, trace = FALSE)

  if (response == "continuous") {
    sigma2 <- nn_0$value / length(y)
    I_0 <- nn_0$Hessian / (2 * sigma2)
  } else if (task == "binary") {
    I_0 <-  nn_0$Hessian
  }

  P <- diag(length(W)) * 2 * lambda

  I_p <- I_0 + P

  vc <- solve(I_p) %*% I_0 %*% solve(I_p)

  return(vc)

}
