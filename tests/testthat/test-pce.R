test_that("analytic prediction gradients agree with finite differences", {
  data <- make_binary_data()
  fit <- fit_binary(data)
  model <- suppressWarnings(interpretnn(fit, data = data))
  x <- model$x[1:5, , drop = FALSE]
  analytic <- statnnet:::.nn_gradient_matrix(
    x,
    model$weights,
    hidden = model$architecture$hidden,
    response = model$response
  )
  step <- 1e-6
  numerical <- sapply(seq_along(model$weights), function(index) {
    upper <- lower <- model$weights
    upper[index] <- upper[index] + step
    lower[index] <- lower[index] - step
    (
      statnnet:::.nn_predict_matrix(x, upper, 1, "binary") -
        statnnet:::.nn_predict_matrix(x, lower, 1, "binary")
    ) / (2 * step)
  })
  expect_equal(unname(analytic), unname(numerical), tolerance = 1e-6)
})

test_that("binary PCE is the direct zero-to-one partial-dependence contrast", {
  set.seed(4107)
  data <- data.frame(
    y = stats::rnorm(100),
    x = stats::rnorm(100),
    z = rep(0:1, 50)
  )
  data$y <- 0.4 + data$x + 0.8 * data$z + stats::rnorm(100, sd = 0.4)
  fit <- nnet::nnet(
    y ~ x + z, data = data, size = 1, linout = TRUE,
    decay = 0.1, Hess = TRUE, maxit = 1000, trace = FALSE
  )
  model <- suppressWarnings(interpretnn(fit, data = data))
  result <- pce(model, "z", uncertainty = "none")
  low <- high <- data
  low$z <- 0
  high$z <- 1
  expected <- mean(statnnet:::.nn_predict_matrix(
    statnnet:::.build_model_matrix(model, high), model$weights, 1, "continuous"
  )) - mean(statnnet:::.nn_predict_matrix(
    statnnet:::.build_model_matrix(model, low), model$weights, 1, "continuous"
  ))

  expect_equal(result$estimate, expected)
  expect_equal(result$contrast, "1 - 0")
})

test_that("delta-method PCE standard errors match direct matrix calculations", {
  data <- make_gaussian_data()
  fit <- fit_gaussian(data)
  model <- suppressWarnings(interpretnn(fit, data = data))
  skip_if_not(model$diagnostics$covariance_available)
  result <- pce(model, "x1", values = 0, d = 0.5)
  scenarios <- statnnet:::.pce_scenarios(model, "x1", 0, 0.5, 2, FALSE)
  x_low <- statnnet:::.build_model_matrix(model, scenarios$items[[1]]$low)
  x_high <- statnnet:::.build_model_matrix(model, scenarios$items[[1]]$high)
  gradient <- colMeans(
    statnnet:::.nn_gradient_matrix(x_high, model$weights, 1, model$response) -
      statnnet:::.nn_gradient_matrix(x_low, model$weights, 1, model$response)
  )
  expected <- sqrt(drop(gradient %*% vcov(model) %*% gradient))

  expect_equal(result$std_error, expected)
})

test_that("simulation intervals are reproducible with a fixed seed", {
  data <- make_gaussian_data()
  fit <- fit_gaussian(data)
  model <- suppressWarnings(interpretnn(fit, data = data))
  skip_if_not(model$diagnostics$covariance_available)
  first <- pce(
    model, "x1", values = c(-0.5, 0.5), d = 0.25,
    uncertainty = "simulation", nsim = 50, seed = 99
  )
  second <- pce(
    model, "x1", values = c(-0.5, 0.5), d = 0.25,
    uncertainty = "simulation", nsim = 50, seed = 99
  )
  expect_equal(first, second)
})

test_that("factor PCEs preserve fitted levels and contrasts", {
  data <- make_gaussian_data()
  contrasts(data$f) <- contr.sum(3)
  fit <- fit_gaussian(data)
  model <- suppressWarnings(interpretnn(fit, data = data))
  result <- pce(model, "f", uncertainty = "none")

  expect_equal(nrow(result), 2L)
  expect_equal(result$contrast, c("b - a", "c - a"))
})
