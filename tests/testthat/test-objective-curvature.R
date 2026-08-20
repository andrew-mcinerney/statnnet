test_that("Gaussian nnet objectives are reconstructed exactly", {
  data <- make_gaussian_data()
  for (decay in c(0, 0.2)) {
    fit <- fit_gaussian(data, decay = decay)
    model <- suppressWarnings(statnnet(fit, data = data, response = "continuous"))
    rss <- sum((data$y - as.numeric(fit$fitted.values))^2)

    expect_equal(fit$value, rss + decay * sum(fit$wts^2), tolerance = 1e-7)
    expect_equal(model$objective$unpenalised, rss, tolerance = 1e-7)
    expect_equal(model$objective$penalty, decay * sum(fit$wts^2), tolerance = 1e-12)
    expect_true(model$objective$agrees)
  }
})

test_that("binary entropy objectives are reconstructed exactly", {
  data <- make_binary_data()
  y <- as.numeric(data$y) - 1
  for (decay in c(0, 0.2)) {
    fit <- fit_binary(data, decay = decay)
    model <- suppressWarnings(statnnet(fit, data = data, response = "binary"))
    probability <- pmin(
      pmax(as.numeric(fit$fitted.values), .Machine$double.eps),
      1 - .Machine$double.eps
    )
    negative_log_likelihood <- -sum(
      y * log(probability) + (1 - y) * log1p(-probability)
    )

    expect_equal(
      fit$value,
      negative_log_likelihood + decay * sum(fit$wts^2),
      tolerance = 1e-7
    )
    expect_equal(model$objective$unpenalised, negative_log_likelihood)
  }
})

test_that("penalty removal and likelihood scaling follow the documented convention", {
  data <- make_gaussian_data()
  decay <- 0.2
  fit <- fit_gaussian(data, decay = decay)
  model <- suppressWarnings(statnnet(fit, data = data))
  expected_h0 <- fit$Hessian - 2 * decay * diag(length(fit$wts))
  sigma2 <- sum((data$y - as.numeric(fit$fitted.values))^2) / nrow(data)

  expect_equal(unname(model$hessian_unpenalised), expected_h0, tolerance = 1e-10)
  expect_equal(
    unname(model$information_unpenalised),
    expected_h0 / (2 * sigma2),
    tolerance = 1e-10
  )

  binary_data <- make_binary_data()
  binary_fit <- fit_binary(binary_data, decay = decay)
  binary_model <- suppressWarnings(statnnet(binary_fit, data = binary_data))
  expect_equal(
    unname(binary_model$information_unpenalised),
    binary_fit$Hessian - 2 * decay * diag(length(binary_fit$wts)),
    tolerance = 1e-10
  )
})

test_that("stored nnet Hessian agrees with a numerical objective Hessian", {
  set.seed(4105)
  data <- data.frame(y = stats::rnorm(30), x = stats::rnorm(30))
  fit <- nnet::nnet(
    y ~ x,
    data = data,
    size = 1,
    linout = TRUE,
    decay = 0.15,
    Hess = TRUE,
    maxit = 1000,
    trace = FALSE
  )
  model <- suppressWarnings(statnnet(fit, data = data))
  objective <- function(weights) {
    prediction <- statnnet:::.nn_predict_matrix(
      model$x, weights, hidden = 1, response = "continuous"
    )
    sum((model$y - prediction)^2) + fit$decay * sum(weights^2)
  }
  weights <- fit$wts
  step <- 1e-4
  numerical <- matrix(0, length(weights), length(weights))
  for (i in seq_along(weights)) {
    for (j in seq_along(weights)) {
      ei <- ej <- rep(0, length(weights))
      ei[i] <- step
      ej[j] <- step
      numerical[i, j] <- (
        objective(weights + ei + ej) - objective(weights + ei - ej) -
          objective(weights - ei + ej) + objective(weights - ei - ej)
      ) / (4 * step^2)
    }
  }
  expect_equal(unname(fit$Hessian), numerical, tolerance = 2e-4)
})

test_that("covariance is symmetric and zero decay reduces to inverse information", {
  data <- make_gaussian_data(160)
  fit <- fit_gaussian(data, decay = 0)
  model <- suppressWarnings(statnnet(fit, data = data))
  skip_if_not(model$diagnostics$covariance_available)

  expect_equal(vcov(model), t(vcov(model)), tolerance = 1e-10)
  expect_equal(
    unname(vcov(model)),
    chol2inv(chol(model$information_unpenalised)),
    tolerance = 1e-7
  )
})

test_that("singular curvature prevents covariance reporting", {
  data <- make_gaussian_data()
  fit <- fit_gaussian(data)
  fit$Hessian[,] <- 0
  expect_warning(
    model <- statnnet(fit, data = data),
    "Covariance-based summaries are unavailable"
  )
  expect_false(model$diagnostics$covariance_available)
  expect_error(vcov(model), "unavailable")
})

test_that("a missing stored Hessian is recomputed from the original data", {
  data <- make_gaussian_data()
  fit <- fit_gaussian(data, Hess = FALSE)
  model <- suppressWarnings(statnnet(fit, data = data))

  expect_identical(model$diagnostics$hessian_source, "recomputed")
  expect_equal(
    unname(model$hessian_penalised),
    nnet::nnetHess(fit, model$x, model$y),
    tolerance = 1e-10
  )
})
