#' Plot partial covariate effects
#'
#' @param x A `statnnet` object.
#' @param variable One or more predictor names. The default plots all supported
#'   predictors.
#' @param uncertainty One of `"delta"`, `"simulation"`, or `"none"`.
#' @param level Pointwise confidence level. Curve intervals are not
#'   simultaneous confidence bands.
#' @param by Optional conditioning predictor. Two PCE curves are drawn at the
#'   observed mean minus and plus one standard deviation for a continuous
#'   predictor, or at both fitted levels for a binary predictor.
#' @param ... Additional arguments passed to [pce()].
#'
#' @return The plotted PCE data, invisibly.
#' @export
plot.statnnet <- function(x, variable = NULL,
                          uncertainty = c("delta", "simulation", "none"),
                          level = 0.95, by = NULL, ...) {
  uncertainty <- match.arg(uncertainty)
  if (is.null(variable)) {
    variable <- x$variables
  }
  if (!is.character(variable) || !length(variable) ||
      any(!variable %in% x$variables)) {
    stop("`variable` must contain predictor names used in the fitted formula.", call. = FALSE)
  }
  if (!is.null(by) && (!is.character(by) || length(by) != 1L || is.na(by) ||
      !by %in% x$variables)) {
    stop("`by` must name one predictor used in the fitted formula.", call. = FALSE)
  }
  if (!is.null(by) && by %in% variable) {
    stop("`by` must be different from every plotted `variable`.", call. = FALSE)
  }

  plotted <- lapply(variable, function(name) {
    pce(
      x,
      variable = name,
      uncertainty = uncertainty,
      level = level,
      by = by,
      ...
    )
  })
  names(plotted) <- variable

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  if (length(plotted) > 1L) {
    columns <- ceiling(sqrt(length(plotted)))
    rows <- ceiling(length(plotted) / columns)
    graphics::par(mfrow = c(rows, columns))
  }

  for (name in variable) {
    values <- plotted[[name]]
    categorical <- all(is.na(values$value)) || any(!is.na(values$contrast))
    conditioned <- !is.null(by)
    if (conditioned) {
      group_labels <- unique(values$by_label)
      group_rows <- lapply(group_labels, function(label) which(values$by_label == label))
      colours <- c("#0072B2", "#D55E00")
      line_types <- c(1, 2)
      limits <- range(
        c(values$conf_low, values$conf_high, values$estimate),
        na.rm = TRUE
      )
      if (!all(is.finite(limits))) limits <- range(values$estimate, na.rm = TRUE)
      if (categorical) {
        first_group <- values[group_rows[[1L]], , drop = FALSE]
        positions <- seq_len(nrow(first_group))
        offsets <- seq(-0.12, 0.12, length.out = length(group_rows))
        graphics::plot(
          positions,
          first_group$estimate,
          type = "n",
          xaxt = "n",
          xlab = name,
          ylab = "Partial covariate effect",
          main = paste(name, "by", by),
          ylim = limits
        )
        graphics::axis(1, at = positions, labels = first_group$contrast)
        for (index in seq_along(group_rows)) {
          group <- values[group_rows[[index]], , drop = FALSE]
          group_positions <- positions + offsets[index]
          graphics::points(
            group_positions, group$estimate, pch = 19, col = colours[index]
          )
          if (any(is.finite(group$conf_low))) {
            graphics::arrows(
              group_positions,
              group$conf_low,
              group_positions,
              group$conf_high,
              angle = 90,
              code = 3,
              length = 0.05,
              col = colours[index]
            )
          }
        }
      } else {
        first_group <- values[group_rows[[1L]], , drop = FALSE]
        graphics::plot(
          first_group$value,
          first_group$estimate,
          type = "l",
          xlab = name,
          ylab = "Partial covariate effect",
          main = paste(name, "by", by),
          ylim = limits,
          col = colours[1L],
          lty = line_types[1L]
        )
        for (index in seq_along(group_rows)) {
          group <- values[group_rows[[index]], , drop = FALSE]
          graphics::lines(
            group$value,
            group$estimate,
            col = colours[index],
            lty = line_types[index]
          )
          if (any(is.finite(group$conf_low))) {
            graphics::lines(group$value, group$conf_low, lty = 3, col = colours[index])
            graphics::lines(group$value, group$conf_high, lty = 3, col = colours[index])
          }
        }
      }
      graphics::legend(
        "topright",
        legend = group_labels,
        col = colours[seq_along(group_labels)],
        lty = line_types[seq_along(group_labels)],
        pch = if (categorical) 19 else NA_integer_,
        bty = "n"
      )
    } else if (categorical) {
      positions <- seq_len(nrow(values))
      limits <- range(c(values$conf_low, values$conf_high, values$estimate), na.rm = TRUE)
      if (!all(is.finite(limits))) limits <- range(values$estimate, na.rm = TRUE)
      graphics::plot(
        positions,
        values$estimate,
        xaxt = "n",
        xlab = name,
        ylab = "Partial covariate effect",
        main = name,
        ylim = limits,
        pch = 19
      )
      graphics::axis(1, at = positions, labels = values$contrast)
      if (any(is.finite(values$conf_low))) {
        graphics::arrows(
          positions,
          values$conf_low,
          positions,
          values$conf_high,
          angle = 90,
          code = 3,
          length = 0.05
        )
      }
    } else {
      limits <- range(c(values$conf_low, values$conf_high, values$estimate), na.rm = TRUE)
      if (!all(is.finite(limits))) limits <- range(values$estimate, na.rm = TRUE)
      graphics::plot(
        values$value,
        values$estimate,
        type = "l",
        xlab = name,
        ylab = "Partial covariate effect",
        main = name,
        ylim = limits
      )
      if (any(is.finite(values$conf_low))) {
        graphics::lines(values$value, values$conf_low, lty = 2)
        graphics::lines(values$value, values$conf_high, lty = 2)
      }
    }
    graphics::abline(h = 0, lty = 3, col = "grey60")
  }
  invisible(if (length(plotted) == 1L) plotted[[1L]] else plotted)
}
