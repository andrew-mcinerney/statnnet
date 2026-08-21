test_that("plotnn maps every multi-node edge to the fitted nnet weights", {
  data <- make_gaussian_data(120)
  fit <- fit_gaussian(data, size = 2L)
  model <- suppressWarnings(statnnet(fit, data = data))
  plot_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file)
  on.exit({
    if (grDevices::dev.cur() > 1L) grDevices::dev.off()
    unlink(plot_file)
  }, add = TRUE)

  edges <- plotnn(model, annotate = "none")
  grDevices::dev.off()

  expect_equal(nrow(edges), 10L)
  expect_equal(
    edges$weight_index,
    c(2L:5L, 7L:10L, 12L:13L)
  )
  expect_equal(edges$estimate, unname(model$weights[edges$weight_index]))
  expect_equal(edges$from[1:4], colnames(model$x))
  expect_equal(edges$to[1:4], rep("h1", 4L))
})

test_that("plotnn supports all documented annotation modes", {
  data <- make_binary_data(140)
  fit <- fit_binary(data, size = 2L)
  model <- suppressWarnings(statnnet(fit, data = data))
  plot_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file)
  on.exit({
    if (grDevices::dev.cur() > 1L) grDevices::dev.off()
    unlink(plot_file)
  }, add = TRUE)

  for (mode in c("estimate", "p_value", "none")) {
    edges <- plotnn(model, annotate = mode, alpha = 0.1)
    expect_true(is.data.frame(edges))
    expect_equal(nrow(edges), model$architecture$hidden * (ncol(model$x) + 1L))
  }
  grDevices::dev.off()

  expect_error(plotnn(model, annotate = "unknown"), "one of")
  expect_error(plotnn(model, alpha = 0), "between zero and one")
  expect_error(plotnn(list()), "must inherit")
})
