library(nnet)

# prep data ---------------------------------------------------------------

set.seed(1)
n <- 500
p <- 4
q <- 2
K <- (p + 2) * q + 1

X <- matrix(rnorm(p * n), ncol = p)

W <- nnic::my_runif(K, 3, 1)

y <- nnic::nn_pred(X, W, q) + rnorm(n)

# nnet --------------------------------------------------------------------

nn <- nnet(X, y, size = q, linout = TRUE, trace = FALSE)


