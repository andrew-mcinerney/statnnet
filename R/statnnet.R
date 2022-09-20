#' statnnet
#'
#'
#' @param nn nnet object
#' @param X matrix of values for nnet
#' @param B number of bootstrap replicates
#' @return statnnet object
#' @export
statnnet <- function(nn, X, B = 1000) {

  if (class(nn) != "nnet") {
    stop("Error: Argument must be of class nnet")
  }

  if (is.null(colnames(X))) {
    colnames(X) <- colnames(X, do.NULL = FALSE, prefix = deparse(substitute(X)))
  }

  n <- nrow(nn$residuals)

  y <- nn$fitted.values + nn$residuals

  class(nn) <- "statnnet"

  nn$BIC <- - 2 * nn_loglike(nn) + 2 * log(n)

  wald <- wald_test(X, y, nn$wts, nn$n[2])

  # Covariate effect and bootstrapped std. error
  nn$eff <- covariate_eff(X, nn$wts, nn$n[2])

  nn$eff_se <- apply(replicate(B,
                               covariate_eff(X[sample(n, size = n, replace = TRUE), ],
                                             W = nn$wts,
                                             q = nn$n[2])),
                     1, sd)


  nn$cl <- match.call()

  nn$wald_p <- wald$p_value
  nn$wald_chi <- wald$chisq

  nn$X <- X
  nn$y <- y
  nn$B <- B

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

  covariates <- colnames(object$X)

  coefdf <- data.frame(
    Covariate = covariates,
    Estimate = object$eff,
    Std.Error = object$eff_se,
    Wald.chi = object$wald_chi,
    Wald.p.value = object$wald_p
  )

  colnames(coefdf)[1] <- ""
  colnames(coefdf)[3] <- "Std. Error"
  colnames(coefdf)[4] <- "  X^2"
  colnames(coefdf)[5] <- "Pr(> X^2)"

  object$coefdf <- coefdf

  Signif <- symnum(object$wald_p, corr = FALSE, na = FALSE,
                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " "))
  object$coefdf$`Pr(> X^2)` <- paste(
    formatC(object$coefdf$`Pr(> X^2)`, format = "e", digits = 2),
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

  print(x$coefdf, right = TRUE, na.print = "NA",
        digits =  max(3L, getOption("digits") - 2L), row.names = FALSE)

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

#' @export
plot.statnnet <-
  function (x, which = c(1L:ncol(x$X)), x_axis_r = c(-3, 3), x_axis_l = 101,
            conf_int = FALSE, alpha = 0.05, B = x$B,
            caption = lapply(1:ncol(x$X),
                             function(iter) paste0("Covariate-Effect Plot for ",
                                                   colnames(x$X)[iter])),
            sub.caption = NULL, main = "",
            ask = prod(par("mfcol")) < length(which) && dev.interactive(), ...,
            label.pos = c(4,2), cex.caption = 1, cex.oma.main = 1.25){

    if (!inherits(x, "statnnet"))
      stop("use only with \"statnnet\" objects")

    if(!is.numeric(which) || any(which < 1) || any(which > ncol(x$X)))
      stop(sprintf("'which' must be in 1:%s", ncol(x$X)))

    if (conf_int == TRUE && is.null(alpha)) {
      stop("'alpha' must be not be NULL when 'conf_int == TRUE'")
    } else if (conf_int == TRUE && (alpha < 0 || alpha > 1 || length(alpha) != 1)) {
      stop("'alpha' must be a value between 0 and 1")
    }

    getCaption <- function(k) # allow caption = "" , plotmath etc
      if(length(caption) < k) NA_character_ else as.graphicsAnnot(caption[[k]])

    show <- rep(FALSE, ncol(x$X))
    show[which] <- TRUE

    cov_effs <- lapply(1:ncol(x$X),
                       function(iter) pdp_effect(x$wts, x$X, x$n[2],
                                                 iter,
                                                 x_r = x_axis_r,
                                                 len = x_axis_l))

    conf_val <- vector("list", length = ncol(x$X))
    for (i in 1:ncol(x$X)) {
      if (conf_int == TRUE && show[i]) {
        conf_val[[i]] <- mlesim(W = x$wts, X = x$X, y = x$y, q = x$n[2], ind =  i,
                                FUN = pdp_effect, B = x$B,
                                x_r = x_axis_r,
                                len = x_axis_l)
      }
    }


    xaxis <- seq(x_axis_r[1], x_axis_r[2], length.out = x_axis_l)

    labs <- colnames(x$X)

    one.fig <- prod(par("mfcol")) == 1
    if (ask) {
      oask <- devAskNewPage(TRUE)
      on.exit(devAskNewPage(oask))
    }
    ##---------- Do the individual plots : ----------
    for (i in 1:ncol(x$X)) {
      if (show[i]) {

        ylim <- range(c(cov_effs[[i]], conf_val[[i]]), na.rm = TRUE)


        if (ylim[1] > 0) ylim[1] = 0 else if (ylim[2] < 0) ylim[2] = 0
        dev.hold()
        plot(xaxis, cov_effs[[i]], xlab = labs[i], ylab = "Effect", main = main,
             ylim = ylim, type = "n", ...)
        lines(xaxis, cov_effs[[i]], ...)
        if (conf_int == TRUE) {
          lines(xaxis, conf_val[[i]]$upper, lty = 2, col = 2, ...)
          lines(xaxis, conf_val[[i]]$lower, lty = 2, col = 2, ...)
        }
        if (one.fig)
          title(sub = sub.caption, ...)
        mtext(getCaption(i), 3, 0.25, cex = cex.caption)
        abline(h = 0, lty = 3, col = "gray")
        dev.flush()
      }
    }
    if (!one.fig && par("oma")[3L] >= 1)
      mtext(sub.caption, outer = TRUE, cex = cex.oma.main)
    invisible()
  }
