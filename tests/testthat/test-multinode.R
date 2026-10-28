test_that("multi-node Gaussian fits preserve predictions and term groups", {
  data <- make_gaussian_data(140)
  fit <- fit_gaussian(data, size = 3L)
  model <- suppressWarnings(statnnet(fit, data = data))

  expect_equal(model$architecture$hidden, 3L)
  expect_false(model$diagnostics$architecture_calibration_warning)
  expect_equal(model$diagnostics$hidden_to_input_ratio, 3 / 4)
  expect_equal(predict(model), predict(fit), tolerance = 1e-10)
  expect_equal(model$groups$x1$weight_indices, c(2L, 7L, 12L))
  expect_equal(model$groups$f$weight_indices, c(4L, 5L, 9L, 10L, 14L, 15L))
  expect_equal(anova(model)$n_weights[anova(model)$term == "f"], 6L)
})

test_that("multi-node binary gradients agree with central differences", {
  data <- make_binary_data(160)
  fit <- fit_binary(data, size = 2L)
  model <- suppressWarnings(statnnet(fit, data = data))
  x <- model$x[1:6, , drop = FALSE]
  analytic <- statnnet:::.nn_gradient_matrix(
    x,
    model$weights,
    hidden = model$architecture$hidden,
    response = model$response
  )
  step <- 1e-6
  numerical <- vapply(seq_along(model$weights), function(index) {
    upper <- lower <- model$weights
    upper[index] <- upper[index] + step
    lower[index] <- lower[index] - step
    (
      statnnet:::.nn_predict_matrix(x, upper, 2L, "binary") -
        statnnet:::.nn_predict_matrix(x, lower, 2L, "binary")
    ) / (2 * step)
  }, numeric(nrow(x)))

  expect_equal(model$architecture$hidden, 2L)
  expect_true(model$diagnostics$binary_calibration_unvalidated)
  expect_match(
    capture.output(print(model)),
    "binary-response Wald p-values",
    all = FALSE
  )
  expect_equal(predict(model), predict(fit), tolerance = 1e-10)
  expect_equal(unname(analytic), unname(numerical), tolerance = 1e-6)
})

test_that("multi-node PCE gradients and variances use the correct weights", {
  data <- make_gaussian_data(120)
  fit <- fit_gaussian(data, size = 2L)
  model <- suppressWarnings(statnnet(fit, data = data))
  skip_if_not(model$diagnostics$covariance_available)

  result <- pce(
    model, "x1", values = 0, d = 0.4, uncertainty = "delta"
  )
  scenarios <- statnnet:::.pce_scenarios(model, "x1", 0, 0.4, 2L, FALSE)
  x_low <- statnnet:::.build_model_matrix(model, scenarios$items[[1L]]$low)
  x_high <- statnnet:::.build_model_matrix(model, scenarios$items[[1L]]$high)
  gradient <- colMeans(
    statnnet:::.nn_gradient_matrix(x_high, model$weights, 2L, model$response) -
      statnnet:::.nn_gradient_matrix(x_low, model$weights, 2L, model$response)
  )
  expected_se <- sqrt(drop(gradient %*% vcov(model) %*% gradient))

  expect_equal(result$std_error, expected_se, tolerance = 1e-10)
})

test_that("ill-conditioned multi-node fits have no covariance summaries", {
  data <- make_gaussian_data(100)
  fit <- fit_gaussian(data, decay = 0.1, size = 2L)
  n_parameters <- length(fit$wts)
  fit$Hessian[,] <- 0
  diag(fit$Hessian) <- c(1, rep(1e-14, n_parameters - 1L))

  expect_warning(
    model <- statnnet(fit, data = data, covariance_tol = 1e-10),
    "Covariance-based summaries are unavailable"
  )
  expect_false(model$diagnostics$covariance_available)
  expect_true(all(is.na(anova(model)$p_value)))
})

test_that("a hidden layer larger than the input layer warns", {
  data <- make_gaussian_data(140)[c("y", "x1")]
  set.seed(4111)
  fit <- nnet::nnet(
    y ~ x1,
    data = data,
    size = 2,
    linout = TRUE,
    decay = 0.1,
    Hess = TRUE,
    maxit = 1000,
    trace = FALSE
  )

  expect_warning(
    model <- statnnet(fit, data = data),
    "more hidden nodes"
  )
  expect_true(model$diagnostics$architecture_calibration_warning)
  expect_equal(model$diagnostics$n_input_nodes, 1L)
  expect_equal(model$diagnostics$n_hidden_nodes, 2L)
  expect_match(
    capture.output(print(model)),
    "Calibration warning",
    all = FALSE
  )
})
