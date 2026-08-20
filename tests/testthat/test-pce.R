test_that("analytic prediction gradients agree with finite differences", {
  data <- make_binary_data()
  fit <- fit_binary(data)
  model <- suppressWarnings(statnnet(fit, data = data))
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
  model <- suppressWarnings(statnnet(fit, data = data))
  result <- pce(model, "z", uncertainty = "none")
  summary_result <- summary(model)
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
  expect_equal(
    summary_result$pce$estimate[summary_result$pce$variable == "z"],
    expected
  )
})

test_that("delta-method PCE standard errors match direct matrix calculations", {
  data <- make_gaussian_data()
  fit <- fit_gaussian(data)
  model <- suppressWarnings(statnnet(fit, data = data))
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
  model <- suppressWarnings(statnnet(fit, data = data))
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
  model <- suppressWarnings(statnnet(fit, data = data))
  result <- pce(model, "f", uncertainty = "none")

  expect_equal(nrow(result), 2L)
  expect_equal(result$contrast, c("b - a", "c - a"))
})

test_that("continuous by predictors condition PCE curves at mean plus or minus one SD", {
  data <- make_gaussian_data()
  fit <- fit_gaussian(data)
  model <- suppressWarnings(statnnet(fit, data = data))
  evaluation_values <- c(-0.5, 0.5)
  result <- pce(
    model,
    "x1",
    values = evaluation_values,
    d = 0.25,
    uncertainty = "none",
    by = "x2"
  )
  expected_by_values <- mean(data$x2) + c(-1, 1) * stats::sd(data$x2)

  expect_equal(nrow(result), 4L)
  expect_equal(result$value, rep(evaluation_values, 2L))
  expect_equal(unique(result$by_value), expected_by_values)
  expect_equal(unique(result$by_variable), "x2")

  low <- high <- data
  low$x1 <- evaluation_values[1L]
  high$x1 <- evaluation_values[1L] + 0.25
  low$x2 <- expected_by_values[1L]
  high$x2 <- expected_by_values[1L]
  expected <- mean(statnnet:::.nn_predict_matrix(
    statnnet:::.build_model_matrix(model, high), model$weights, 1, "continuous"
  )) - mean(statnnet:::.nn_predict_matrix(
    statnnet:::.build_model_matrix(model, low), model$weights, 1, "continuous"
  ))
  expect_equal(result$estimate[1L], expected)

  if (isTRUE(model$diagnostics$covariance_available)) {
    delta_result <- pce(
      model,
      "x1",
      values = evaluation_values,
      d = 0.25,
      uncertainty = "delta",
      by = "x2"
    )
    expect_true(all(is.finite(delta_result$std_error)))
  }
})

test_that("binary by predictors condition PCE curves at zero and one", {
  set.seed(4108)
  data <- data.frame(
    y = stats::rnorm(100),
    x = stats::rnorm(100),
    z = rep(0:1, 50)
  )
  data$y <- 0.5 + data$x + data$x * data$z + stats::rnorm(100, sd = 0.4)
  fit <- nnet::nnet(
    y ~ x + z, data = data, size = 1, linout = TRUE,
    decay = 0.1, Hess = TRUE, maxit = 1000, trace = FALSE
  )
  model <- suppressWarnings(statnnet(fit, data = data))
  result <- pce(
    model, "x", length_out = 5L, uncertainty = "none", by = "z"
  )

  expect_equal(nrow(result), 10L)
  expect_equal(unique(result$by_value), c(0, 1))
  expect_equal(unique(result$by_label), c("z = 0", "z = 1"))
})

test_that("plot returns both conditioned curves", {
  data <- make_gaussian_data()
  fit <- fit_gaussian(data)
  model <- suppressWarnings(statnnet(fit, data = data))
  plot_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file)
  on.exit({
    if (grDevices::dev.cur() > 1L) grDevices::dev.off()
    unlink(plot_file)
  }, add = TRUE)

  result <- plot(
    model,
    variable = "x1",
    by = "x2",
    uncertainty = "none",
    length_out = 5L
  )
  grDevices::dev.off()

  expect_s3_class(result, "statnnet_pce")
  expect_equal(nrow(result), 10L)
  expect_equal(length(unique(result$by_value)), 2L)
})
