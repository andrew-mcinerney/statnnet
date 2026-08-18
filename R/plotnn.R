#' Plot a statnnet network with statistical annotations
#'
#' Draws the supported single-hidden-layer network. Input-to-hidden edges can
#' be annotated with individual-weight estimates or Wald p-values. Individual
#' weight results are diagnostics; grouped results from [anova.statnnet()] are
#' the intended covariate-level summaries.
#'
#' @param x A `statnnet` object.
#' @param annotate One of `"estimate"`, `"p_value"`, or `"none"`.
#' @param alpha Significance threshold used to colour p-value annotations.
#' @param digits Number of annotation digits.
#' @param ... Additional graphical parameters passed to [graphics::plot()].
#'
#' @return Edge information, invisibly.
#' @export
plotnn <- function(x, annotate = c("estimate", "p_value", "none"),
                   alpha = 0.05, digits = 3L, ...) {
  if (!inherits(x, "statnnet")) {
    stop("`x` must inherit from \"statnnet\".", call. = FALSE)
  }
  annotate <- match.arg(annotate)
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be between zero and one.", call. = FALSE)
  }

  n_inputs <- x$architecture$inputs
  n_hidden <- x$architecture$hidden
  input_y <- seq(0.9, 0.1, length.out = n_inputs)
  hidden_y <- seq(0.8, 0.2, length.out = n_hidden)
  output_y <- 0.5
  weight_table <- .parameter_wald(x)

  graphics::plot(
    c(0, 1),
    c(0, 1),
    type = "n",
    axes = FALSE,
    xlab = "",
    ylab = "",
    asp = 1,
    ...
  )
  graphics::text(0.05, input_y, labels = colnames(x$x), pos = 2, xpd = TRUE)
  graphics::points(rep(0.1, n_inputs), input_y, pch = 21, bg = "white", cex = 1.8)
  graphics::points(rep(0.55, n_hidden), hidden_y, pch = 21, bg = "grey90", cex = 2)
  graphics::points(0.95, output_y, pch = 21, bg = "white", cex = 2)
  graphics::text(rep(0.55, n_hidden), hidden_y, labels = paste0("h", seq_len(n_hidden)))
  graphics::text(0.95, output_y, labels = "y")

  edge_rows <- list()
  row <- 1L
  for (hidden in seq_len(n_hidden)) {
    for (input in seq_len(n_inputs)) {
      index <- (hidden - 1L) * (n_inputs + 1L) + input + 1L
      p_value <- weight_table$p_value[index]
      colour <- if (annotate == "p_value" && is.finite(p_value) && p_value < alpha) {
        "firebrick"
      } else {
        "grey45"
      }
      graphics::segments(0.12, input_y[input], 0.53, hidden_y[hidden], col = colour)
      label <- switch(
        annotate,
        estimate = formatC(x$weights[index], digits = digits, format = "fg"),
        p_value = if (is.finite(p_value)) {
          paste0("p=", formatC(p_value, digits = digits, format = "g"))
        } else {
          "p=NA"
        },
        none = ""
      )
      if (nzchar(label)) {
        graphics::text(
          0.325,
          (input_y[input] + hidden_y[hidden]) / 2,
          labels = label,
          cex = 0.65,
          col = colour
        )
      }
      edge_rows[[row]] <- data.frame(
        from = colnames(x$x)[input],
        to = paste0("h", hidden),
        weight_index = index,
        estimate = x$weights[index],
        p_value = p_value,
        stringsAsFactors = FALSE
      )
      row <- row + 1L
    }
  }
  output_start <- (n_inputs + 1L) * n_hidden
  for (hidden in seq_len(n_hidden)) {
    index <- output_start + hidden + 1L
    graphics::segments(0.57, hidden_y[hidden], 0.93, output_y, col = "grey45")
    edge_rows[[row]] <- data.frame(
      from = paste0("h", hidden),
      to = "y",
      weight_index = index,
      estimate = x$weights[index],
      p_value = weight_table$p_value[index],
      stringsAsFactors = FALSE
    )
    row <- row + 1L
  }
  graphics::title(
    main = sprintf(
      "statnnet %d-%d-1 network",
      n_inputs,
      n_hidden
    )
  )
  invisible(do.call(rbind, edge_rows))
}
