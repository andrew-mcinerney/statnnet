test_that("multi-node Gaussian fits preserve prediction and term mappings", {
  data <- make_gaussian_data(140)
  fit <- fit_gaussian(data, size = 3L)
  model <- suppressWarnings(statnnet(fit, data = data))

  expect_equal(model$architecture$hidden, 3L)
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
  expect_equal(predict(model), predict(fit), tolerance = 1e-10)
  expect_equal(unname(analytic), unname(numerical), tolerance = 1e-6)
})

test_that("multi-node PCE gradients and delta variances are mapped correctly", {
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

test_that("ill-conditioned multi-node curvature disables covariance summaries", {
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
