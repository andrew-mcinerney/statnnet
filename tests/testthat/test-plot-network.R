test_that("network edges map to the fitted weights", {
  data <- make_gaussian_data(120)
  fit <- fit_gaussian(data, size = 2L)
  model <- suppressWarnings(statnnet(fit, data = data))
  plot_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file)
  on.exit({
    if (grDevices::dev.cur() > 1L) grDevices::dev.off()
    unlink(plot_file)
  }, add = TRUE)

  edges <- plot_network(model, annotate = "none")
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

test_that("plot_network supports all documented annotation modes", {
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
    edges <- plot_network(model, annotate = mode, alpha = 0.1)
    expect_true(is.data.frame(edges))
    expect_equal(nrow(edges), model$architecture$hidden * (ncol(model$x) + 1L))
  }
  grDevices::dev.off()

  expect_error(plot_network(model, annotate = "unknown"), "one of")
  expect_error(plot_network(model, alpha = 0), "between zero and one")
  expect_error(plot_network(list()), "must inherit")
})
