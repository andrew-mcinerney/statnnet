make_gaussian_data <- function(n = 100L) {
  set.seed(4101)
  data <- data.frame(
    y = stats::rnorm(n),
    x1 = stats::rnorm(n),
    x2 = stats::runif(n, -1, 1),
    f = factor(rep(c("a", "b", "c", "a"), length.out = n))
  )
  data$y <- 1 + 1.2 * data$x1 - 0.7 * data$x2 +
    c(a = 0, b = 0.4, c = -0.3)[data$f] + stats::rnorm(n, sd = 0.5)
  data
}

fit_gaussian <- function(data, decay = 0.1, Hess = TRUE, size = 1L) {
  set.seed(4102)
  nnet::nnet(
    y ~ x1 + x2 + f,
    data = data,
    size = size,
    linout = TRUE,
    decay = decay,
    Hess = Hess,
    maxit = 1000,
    trace = FALSE
  )
}

make_binary_data <- function(n = 120L) {
  set.seed(4103)
  x <- stats::rnorm(n)
  f <- factor(rep(c("low", "middle", "high"), length.out = n))
  probability <- stats::plogis(-0.2 + x + c(low = -0.4, middle = 0, high = 0.5)[f])
  data.frame(
    y = factor(ifelse(stats::runif(n) < probability, "yes", "no"),
               levels = c("no", "yes")),
    x = x,
    f = f
  )
}

fit_binary <- function(data, decay = 0.1, Hess = TRUE, size = 1L) {
  set.seed(4104)
  nnet::nnet(
    y ~ x + f,
    data = data,
    size = size,
    decay = decay,
    Hess = Hess,
    maxit = 1000,
    trace = FALSE
  )
}
