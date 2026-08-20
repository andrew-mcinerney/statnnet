test_that("the original fit is retained and prediction delegates to nnet", {
  data <- make_gaussian_data()
  fit <- fit_gaussian(data)
  model <- suppressWarnings(statnnet(fit, data = data))
  newdata <- data[1:8, ]
  newdata$f <- factor(newdata$f, levels = rev(levels(data$f)))

  expect_identical(model$fit, fit)
  expect_equal(predict(model), predict(fit))
  expect_equal(predict(model, newdata = newdata), predict(fit, newdata = newdata))
})

test_that("fitted formula transformations are retained for new data", {
  set.seed(4110)
  data <- data.frame(y = stats::rnorm(80), x = stats::rnorm(80))
  data$y <- 0.5 + data$x - 0.2 * data$x^2 + stats::rnorm(80, sd = 0.3)
  fit <- nnet::nnet(
    y ~ poly(x, 2),
    data = data,
    size = 1,
    linout = TRUE,
    decay = 0.1,
    Hess = TRUE,
    maxit = 1000,
    trace = FALSE
  )
  model <- suppressWarnings(statnnet(fit, data = data))
  newdata <- data.frame(x = c(-1, 0, 1))

  expect_equal(predict(model, newdata), predict(fit, newdata))
  expect_equal(
    colnames(statnnet:::.build_model_matrix(model, newdata)),
    fit$coefnames
  )
})

test_that("factor columns are grouped from model-matrix assignments", {
  data <- make_gaussian_data()
  contrasts(data$f) <- contr.sum(3)
  fit <- fit_gaussian(data)
  model <- suppressWarnings(statnnet(fit, data = data))
  grouped <- anova(model)
  factor_row <- grouped[grouped$term == "f", ]

  expect_equal(length(model$groups$f$column_indices), 2L)
  expect_equal(length(model$groups$f$weight_indices), 2L)
  expect_equal(factor_row$n_weights, 2L)
  expect_match(factor_row$model_columns, "f")

  absent_level_data <- droplevels(subset(data, f != "c"))[1:5, ]
  expect_equal(
    predict(model, newdata = absent_level_data),
    predict(fit, newdata = absent_level_data)
  )
})

test_that("grouped Wald calculations use effective degrees of freedom", {
  data <- make_gaussian_data()
  fit <- fit_gaussian(data)
  model <- suppressWarnings(statnnet(fit, data = data))
  skip_if_not(model$diagnostics$covariance_available)
  result <- anova(model)
  group <- model$groups$x1
  indices <- group$weight_indices
  group_covariance <- vcov(model)[indices, indices, drop = FALSE]
  theta <- coef(model)[indices]
  expected_statistic <- drop(crossprod(theta, solve(group_covariance, theta)))
  expected_df <- sum(diag(model$shrinkage[indices, indices, drop = FALSE]))

  expect_equal(result$wald_chisq[result$term == "x1"], expected_statistic)
  expect_equal(result$effective_df[result$term == "x1"], expected_df)
})

test_that("unsupported nnet structures and backends fail early", {
  data <- make_gaussian_data(40)
  x <- model.matrix(y ~ x1 + x2, data)[, -1, drop = FALSE]
  set.seed(4106)
  skip_fit <- nnet::nnet(
    x, data$y, size = 1, linout = TRUE, skip = TRUE, trace = FALSE
  )
  expect_error(
    statnnet(skip_fit, formula = y ~ x1 + x2, data = data),
    "Skip-layer"
  )
  expect_error(statnnet(list()), "must inherit")

  fit <- fit_gaussian(data)
  fit$call$mask <- quote(rep(TRUE, length(wts)))
  expect_error(statnnet(fit, data = data), "mask")

  class_data <- transform(data, y = factor(rep(c("a", "b", "c", "a"), 10)))
  set.seed(4108)
  softmax_fit <- nnet::nnet(y ~ x1, data = class_data, size = 1, trace = FALSE)
  expect_error(statnnet(softmax_fit, data = class_data), "Softmax")
})

test_that("non-unit case weights are rejected", {
  data <- make_gaussian_data()
  data$case_weight <- rep(c(1, 2), length.out = nrow(data))
  set.seed(4109)
  fit <- nnet::nnet(
    y ~ x1 + x2,
    data = data,
    weights = case_weight,
    size = 1,
    linout = TRUE,
    Hess = TRUE,
    trace = FALSE
  )
  expect_error(
    statnnet(fit, data = data),
    "Case weights other than one"
  )
})

test_that("non-convergence disables covariance summaries", {
  data <- make_gaussian_data()
  fit <- fit_gaussian(data)
  fit$convergence <- 1L
  expect_warning(
    model <- statnnet(fit, data = data),
    "did not converge"
  )
  expect_false(model$diagnostics$covariance_available)
  expect_true(all(is.na(anova(model)$p_value)))
})

test_that("summary and anova return structured objects", {
  data <- make_gaussian_data()
  fit <- fit_gaussian(data)
  model <- suppressWarnings(statnnet(fit, data = data))
  result <- summary(model, weights = TRUE)

  expect_s3_class(result, "summary.statnnet")
  expect_s3_class(result$grouped_wald, "anova.statnnet")
  expect_true(is.data.frame(result$pce))
  expect_true(is.data.frame(result$parameter_wald))
})

test_that("summary uses fast observed-row effects by default", {
  data <- make_gaussian_data()
  fit <- fit_gaussian(data)
  model <- suppressWarnings(statnnet(fit, data = data))
  result <- summary(model)

  d <- stats::sd(data$x1)
  high <- low <- data
  high$x1 <- high$x1 + d
  expected <- mean(
    statnnet:::.nn_predict_matrix(
      statnnet:::.build_model_matrix(model, high),
      model$weights, 1, model$response
    ) - statnnet:::.nn_predict_matrix(
      statnnet:::.build_model_matrix(model, low),
      model$weights, 1, model$response
    )
  )

  expect_identical(result$effects, "rowwise")
  expect_equal(result$pce$estimate[result$pce$variable == "x1"], expected)
})

test_that("summary effect modes preserve exact PCEs or omit effects", {
  data <- make_gaussian_data(40)
  fit <- fit_gaussian(data)
  model <- suppressWarnings(statnnet(fit, data = data))

  partial <- summary(model, effects = "partial")
  expected <- do.call(rbind, lapply(model$variables, function(variable) {
    pce(model, variable, type = "average", uncertainty = "delta")
  }))
  rownames(expected) <- NULL
  omitted <- summary(model, effects = "none")

  expect_equal(partial$pce, expected)
  expect_identical(partial$effects, "partial")
  expect_identical(omitted$effects, "none")
  expect_equal(nrow(omitted$pce), 0L)
  expect_error(summary(model, effects = "unknown"), "one of")
})
