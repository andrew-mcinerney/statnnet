#' statnnet
#'
#'
#' @param nnet nnet object
#' @return statnnet object
#' @export
statnnet <- function(nn, X) {

  if (class(nn) != "nnet") {
    stop("Error: Argument must be of class nnet")
  }

  if (is.null(colnames(X))) {
    colnames(X) <- colnames(X, do.NULL = FALSE, prefix = deparse(substitute(X)))
  }

  wald <- wald_test(X, nn$fitted.values + nn$residuals, nn$wts, nn$n[2])

  class(nn) <- "statnnet"

  nn$cl <- match.call()

  nn$wald_p <- wald$p_value

  nn$X <- X

  return(nn)

}

#' @export
print.statnnet <- function(x, ...) {
  cat("Call (nnet):\n")
  print(x$call)
  cat("Call (statnnet):\n")
  print(x$cl)
  cat("\n")
  cat("Model Architecture: ", x$n[1], "-", x$n[2], "-", "1", " network",
      sep = ""
  )
  cat(" with", (x$n[1] + 2) * x$n[2] + 1, "weights\n")
}

#' @export
coef.statnnet <- function(object, ...) {
  wts <- object$wts
  p <- object$n[1]
  q <- object$n[2]

  wm <- c(
    "b",
    paste(colnames(object$X), sep = ""),
    paste("h", seq_len(q), sep = ""),
    as.character(object$call$y)
  )

  conn <- c(rep(0:p, times = q), 0, (p + 1):(p + q))

  nunits <- p + q + 2

  nconn <- c(
    rep(0, times = p + 2),
    seq(p + 1, q * (p + 1), by = (p + 1)),
    q * (p + 2) + 1
  )

  names(wts) <- apply(
    cbind(
      wm[1 + conn],
      wm[1 + rep(1:nunits - 1, diff(nconn))]
    ),
    1,
    function(x) paste(x, collapse = "->")
  )
  return(wts)
}

#' @export
summary.statnnet <- function(object, ...) {
  p <- object$n[1]
  q <- object$n[2]

  nconn <- c(
    rep(0, times = p + 2),
    seq(p + 1, q * (p + 1), by = (p + 1)),
    q * (p + 2) + 1
  )

  object$nconn <- nconn

  object$BIC <- - 2 * nn_loglike(object) + 2 * log(nrow(nn$residuals))

  covariates <- colnames(object$X)

  eff <- covariate_eff(object$X, object$wts, object$n[2])

  coefdf <- data.frame(
    Covariate = covariates,
    Estimate = eff,
    Wald.p.value = object$wald_p
  )

  colnames(coefdf)[1] <- ""

  object$coefdf <- coefdf

  Signif <- symnum(object$wald_p, corr = FALSE, na = FALSE,
                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " "))
  object$coefdf$Wald.p.value <- paste(
    formatC(object$coefdf$Wald.p.value, format = "e", digits = 2),
    format(Signif))

  class(object) <- c("summary.statnnet", class(object))
  return(object)
}

#' @export
print.summary.statnnet <- function(x, ...) {
  cat("Call (nnet):\n")
  print(x$call)
  cat("Call (statnnet):\n")
  print(x$cl)
  cat("\n")
  cat("Number of input nodes:", x$n[1], "\n")
  cat("Number of hidden nodes:", x$n[2], "\n")
  cat("\n")
  cat("BIC:", x$BIC, "\n")
  cat("\n")
  cat("Coefficients:\n")

  print(x$coefdf, right = TRUE, na.print = "NA", digits = 2, row.names = FALSE)

  Signif <- symnum(x$wald_p, corr = FALSE, na = FALSE,
                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " "))
  if ((w <- getOption("width")) < nchar(sleg <- attr(Signif,
                                                     "legend")))
    sleg <- strwrap(sleg, width = w - 2, prefix = "  ")
  cat("---\nSignif. codes:  ", sleg, sep = "", fill = w +
        4 + max(nchar(sleg, "bytes") - nchar(sleg)))
  cat("\n")
  cat("Weights:\n")
  wts <- format(round(coef.statnnet(x), 2))
  lapply(
    split(wts, rep(1:(x$n[1] + x$n[2] + 2), diff(x$nconn))),
    function(x) print(x, quote = FALSE)
  )
}
