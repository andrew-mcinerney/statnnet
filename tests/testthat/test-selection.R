make_selection_data <- function(n = 80L) {
  set.seed(5101)
  data <- data.frame(
    x1 = stats::rnorm(n),
    x2 = stats::rnorm(n),
    f = factor(rep(c("a", "b", "c", "a"), length.out = n))
  )
  data$y <- 1 + 1.8 * data$x1 + ifelse(data$f == "b", 0.8, 0) +
    stats::rnorm(n, sd = 0.45)
  data
}

test_that("selection stores the fitted model and search history", {
  data <- make_selection_data()
  selection <- select_nnet(
    y ~ x1 + x2 + f,
    data = data,
    max_size = 2,
    n_starts = 2,
    seed = 5102,
    progress = FALSE,
    maxit = 500
  )

  expect_s3_class(selection, "statnnet_selection")
  expect_s3_class(selection$fit, "nnet")
  expect_identical(selection$response, "continuous")
  expect_equal(selection$fit$decay, 0)
  expect_true(selection$fit$convergence == 0L)
  expect_true(selection$size %in% 1:3)
  expect_true(all(selection$selected_terms %in% c("x1", "x2", "f")))
  expect_setequal(
    c(selection$selected_terms, selection$excluded_terms),
    c("x1", "x2", "f")
  )

  path <- selection_path(selection)
  expect_true(is.data.frame(path))
  expect_true(all(c(
    "phase", "candidate", "included_terms", "size", "n_parameters",
    "criterion_value", "converged_starts", "selected"
  ) %in% names(path)))
  expect_equal(sum(path$selected), 1L)
  expect_equal(path$criterion_value[path$selected], selection$criterion_value)
  expect_true(all(selection$start_diagnostics$evaluation %in% path$evaluation))
})

test_that("fixed seeds reproduce selection", {
  data <- make_selection_data(60)
  first <- select_nnet(
    y ~ x1 + x2,
    data = data,
    max_size = 1,
    n_starts = 1,
    seed = 5103,
    progress = FALSE,
    maxit = 300
  )
  second <- select_nnet(
    y ~ x1 + x2,
    data = data,
    max_size = 1,
    n_starts = 1,
    seed = 5103,
    progress = FALSE,
    maxit = 300
  )

  expect_equal(first$criterion_value, second$criterion_value)
  expect_equal(first$selected_terms, second$selected_terms)
  expect_equal(first$fit$wts, second$fit$wts)
  progress_messages <- capture.output(
    invisible(
      select_nnet(
        y ~ x1,
        data = data,
        max_size = 1,
        n_starts = 1,
        seed = 5104,
        progress = TRUE,
        maxit = 300
      )
    ),
    type = "message"
  )
  expect_true(any(grepl("candidate", progress_messages, fixed = TRUE)))
})

test_that("factor columns are selected together", {
  data <- make_selection_data()
  prepared <- statnnet:::.prepare_selection_formula(
    y ~ x1 + f,
    data,
    environment()
  )

  expect_equal(names(prepared$groups), c("x1", "f"))
  expect_equal(length(prepared$groups$f), nlevels(data$f) - 1L)

  selection <- select_nnet(
    y ~ x1 + f,
    data = data,
    max_size = 1,
    n_starts = 1,
    seed = 5105,
    progress = FALSE,
    maxit = 300
  )
  expect_false(any(grepl("fb|fc", selection_path(selection)$included_terms)))
})

test_that("factor predictions retain levels and contrasts", {
  set.seed(5109)
  data <- data.frame(
    f = factor(rep(c("a", "b", "c"), 30L))
  )
  data$y <- c(a = 0, b = 2, c = -1)[data$f] +
    stats::rnorm(nrow(data), sd = 0.2)
  selection <- select_nnet(
    y ~ f,
    data = data,
    max_size = 1,
    n_starts = 2,
    seed = 5110,
    progress = FALSE,
    maxit = 500
  )
  model <- suppressWarnings(statnnet(selection))
  newdata <- data.frame(f = factor(c("c", "a"), levels = levels(data$f)))

  expect_equal(
    predict(selection, newdata = newdata),
    predict(model, newdata = newdata)
  )
})

test_that("transformed terms retain their prediction basis", {
  set.seed(5111)
  data <- data.frame(x = seq(-2, 2, length.out = 80L))
  data$y <- 0.5 + data$x - 0.8 * data$x^2 + stats::rnorm(80L, sd = 0.2)
  selection <- select_nnet(
    y ~ poly(x, 2),
    data = data,
    max_size = 1,
    n_starts = 2,
    seed = 5112,
    progress = FALSE,
    maxit = 500
  )
  model <- suppressWarnings(statnnet(selection))
  newdata <- data.frame(x = c(-1.25, 0.25, 1.5))

  expect_equal(
    predict(selection, newdata = newdata),
    predict(model, newdata = newdata)
  )
})

test_that("selected models convert to statnnet objects", {
  data <- make_selection_data(90)
  selection <- select_nnet(
    y ~ x1 + x2,
    data = data,
    max_size = 1,
    n_starts = 2,
    seed = 5106,
    progress = FALSE,
    maxit = 500
  )
  model <- suppressWarnings(statnnet(selection))

  expect_s3_class(model, "statnnet")
  expect_identical(model$selection$criterion, selection$criterion)
  expect_equal(model$selection$criterion_value, selection$criterion_value)
  expect_equal(predict(model), predict(selection))
  expect_equal(
    predict(model, newdata = data[1:5, , drop = FALSE]),
    predict(selection, newdata = data[1:5, , drop = FALSE])
  )
  expect_true(any(grepl(
    "conditional on the selected model",
    capture.output(print(model)),
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "conditional on the selected model",
    capture.output(print(summary(model, effects = "none"))),
    fixed = TRUE
  )))
})

test_that("selection methods use the retained fit", {
  data <- make_selection_data(70)
  selection <- select_nnet(
    y ~ x1 + x2,
    data = data,
    max_size = 1,
    n_starts = 1,
    seed = 5107,
    progress = FALSE,
    maxit = 300
  )

  expect_equal(unname(coef(selection)), unname(coef(selection$fit)))
  expect_equal(
    predict(selection, newdata = data[1:5, , drop = FALSE]),
    predict(
      selection$fit,
      statnnet:::.selection_newdata_matrix(
        selection,
        data[1:5, , drop = FALSE]
      )
    )
  )
  expect_s3_class(summary(selection), "summary.statnnet_selection")
  expect_true(any(grepl(
    "Selected model",
    capture.output(print(selection)),
    fixed = TRUE
  )))

  plot_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file)
  on.exit({
    if (grDevices::dev.cur() > 1L) grDevices::dev.off()
    unlink(plot_file)
  }, add = TRUE)
  plotted <- plot(selection)
  grDevices::dev.off()
  expect_equal(plotted, selection_path(selection))
})

test_that("matrix inputs are selected by column", {
  data <- make_selection_data(70)
  x <- as.matrix(data[c("x1", "x2")])
  selection <- select_nnet(
    x,
    data$y,
    max_size = 1,
    n_starts = 1,
    seed = 5108,
    progress = FALSE,
    maxit = 300
  )

  expect_identical(selection$interface, "matrix")
  expect_true(all(selection$selected_terms %in% colnames(x)))
  expect_equal(
    length(predict(selection, newdata = x[1:4, , drop = FALSE])),
    4L
  )
  expect_s3_class(suppressWarnings(statnnet(selection)), "statnnet")
})

test_that("matrix input omits incomplete rows consistently", {
  data <- make_selection_data(50)
  x <- as.matrix(data[c("x1", "x2")])
  x[7, 2] <- NA_real_
  selection <- select_nnet(
    x,
    data$y,
    max_size = 1,
    n_starts = 1,
    seed = 5113,
    progress = FALSE,
    maxit = 300
  )

  expect_equal(nrow(selection$training_data), 49L)
  expect_s3_class(suppressWarnings(statnnet(selection)), "statnnet")
})

test_that("invalid selection requests are rejected", {
  data <- make_selection_data(40)
  class_data <- transform(data, y = factor(ifelse(y > median(y), "a", "b")))

  expect_error(
    select_nnet(y ~ x1, class_data, max_size = 1),
    "numeric Gaussian"
  )
  expect_error(
    select_nnet(y ~ x1 * x2, data, max_size = 1),
    "Interaction terms"
  )
  expect_error(
    select_nnet(y ~ x1 - 1, data, max_size = 1),
    "include an intercept"
  )
  expect_error(
    select_nnet(y ~ x1, data, max_size = 1, decay = 0.1),
    "controlled"
  )
  expect_error(
    select_nnet(y ~ x1, data, max_size = 1, weights = rep(1, nrow(data))),
    "controlled"
  )
  expect_error(
    select_nnet(y ~ x1, data, max_size = 0),
    "positive integer"
  )
  expect_error(selection_path(list()), "statnnet_selection")
})
