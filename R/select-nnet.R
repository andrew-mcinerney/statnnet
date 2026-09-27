#' Select a parsimonious nnet model
#'
#' `select_nnet()` implements the information-criterion procedure of
#' McInerney and Burke: hidden-layer size selection, backward input selection,
#' and alternating fine-tuning of the hidden and input layers. Candidate models
#' are fitted from several random initialisations using [nnet::nnet()].
#'
#' The published procedure is currently implemented only for unpenalised
#' Gaussian regression. In the formula interface, complete right-hand-side
#' terms are selected together, so all design-matrix columns for a factor or a
#' transformed term enter or leave the model as a group. In the matrix
#' interface, each column is a selectable input.
#'
#' @param x A model formula or a numeric design matrix.
#' @param ... Arguments passed to methods and then to [nnet::nnet()]. Arguments
#'   that would change the supported architecture, response model, initial
#'   weights, or decay value are not allowed.
#' @param data A data frame used with the formula interface.
#' @param y A numeric response used with the matrix interface.
#' @param max_size Maximum number of hidden units considered in the initial
#'   hidden-layer search.
#' @param n_starts Number of random initialisations fitted for each candidate.
#' @param criterion Information criterion used to compare candidates.
#' @param start_range Numeric vector of length two giving the lower and upper
#'   bounds for uniformly generated starting weights.
#' @param maxit Maximum number of `nnet` optimiser iterations per start.
#' @param seed Optional random seed for reproducible starting weights.
#' @param progress Whether to print progress as candidate models are fitted.
#'
#' @return An object of class `"statnnet_selection"` containing the selected
#'   native `nnet` fit, selected terms or columns, the full search path, and
#'   start-level convergence diagnostics.
#' @references
#' McInerney, A. and Burke, K. (2025). A statistical modelling approach to
#' feedforward neural network model selection. *Statistical Modelling*, 25(4),
#' 323--342. \doi{10.1177/1471082X241258261}.
#' @examples
#' set.seed(1)
#' example_data <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
#' example_data$y <- 1 + 2 * example_data$x1 + rnorm(60, sd = 0.4)
#'
#' selection <- select_nnet(
#'   y ~ x1 + x2,
#'   data = example_data,
#'   max_size = 1,
#'   n_starts = 2,
#'   maxit = 500,
#'   seed = 2,
#'   progress = FALSE
#' )
#' selection
#' @export
select_nnet <- function(x, ...) {
  UseMethod("select_nnet")
}

#' @rdname select_nnet
#' @export
select_nnet.formula <- function(x, data, ..., max_size, n_starts = 5L,
                                criterion = c("BIC", "AIC", "AICc"),
                                start_range = c(-3, 3), maxit = 1000L,
                                seed = NULL, progress = TRUE) {
  if (missing(data)) {
    stop("`data` is required for the formula interface.", call. = FALSE)
  }
  call <- match.call()
  call[[1L]] <- quote(select_nnet)
  prepared <- .prepare_selection_formula(x, data, parent.frame())
  .select_nnet_impl(
    prepared = prepared,
    max_size = max_size,
    n_starts = n_starts,
    criterion = criterion,
    start_range = start_range,
    maxit = maxit,
    seed = seed,
    progress = progress,
    nnet_args = list(...),
    call = call
  )
}

#' @rdname select_nnet
#' @export
select_nnet.default <- function(x, y, ..., max_size, n_starts = 5L,
                                criterion = c("BIC", "AIC", "AICc"),
                                start_range = c(-3, 3), maxit = 1000L,
                                seed = NULL, progress = TRUE) {
  if (missing(y)) {
    stop("`y` is required for the matrix interface.", call. = FALSE)
  }
  call <- match.call()
  call[[1L]] <- quote(select_nnet)
  prepared <- .prepare_selection_matrix(x, y)
  .select_nnet_impl(
    prepared = prepared,
    max_size = max_size,
    n_starts = n_starts,
    criterion = criterion,
    start_range = start_range,
    maxit = maxit,
    seed = seed,
    progress = progress,
    nnet_args = list(...),
    call = call
  )
}

.prepare_selection_formula <- function(formula, data, environment) {
  if (!inherits(formula, "formula")) {
    stop("`x` must be a formula or numeric matrix.", call. = FALSE)
  }
  data <- as.data.frame(data)
  terms_full <- stats::terms(formula, data = data)
  if (attr(terms_full, "response") != 1L) {
    stop("The selection formula must contain a response.", call. = FALSE)
  }
  if (attr(terms_full, "intercept") != 1L) {
    stop(
      "The selection formula must include an intercept.",
      call. = FALSE
    )
  }
  if (any(attr(terms_full, "order") > 1L)) {
    stop(
      "Interaction terms are not currently supported by `select_nnet()`.",
      call. = FALSE
    )
  }
  term_labels <- attr(terms_full, "term.labels")
  if (!length(term_labels)) {
    stop("The selection formula must contain at least one predictor term.", call. = FALSE)
  }

  model_frame_all <- stats::model.frame(
    terms_full,
    data = data,
    na.action = stats::na.pass,
    drop.unused.levels = TRUE
  )
  complete <- stats::complete.cases(model_frame_all)
  model_frame <- model_frame_all[complete, , drop = FALSE]
  training_data <- data[complete, , drop = FALSE]
  if (!nrow(model_frame)) {
    stop("No complete observations are available for selection.", call. = FALSE)
  }
  response <- stats::model.response(model_frame)
  .validate_selection_response(response)

  model_matrix_full <- stats::model.matrix(terms_full, data = model_frame)
  intercept <- match("(Intercept)", colnames(model_matrix_full), nomatch = 0L)
  keep <- if (intercept) {
    seq_len(ncol(model_matrix_full)) != intercept
  } else {
    rep(TRUE, ncol(model_matrix_full))
  }
  design <- model_matrix_full[, keep, drop = FALSE]
  assign <- attr(model_matrix_full, "assign")[keep]
  groups <- lapply(seq_along(term_labels), function(index) {
    which(assign == index)
  })
  names(groups) <- term_labels
  if (any(lengths(groups) == 0L)) {
    stop("At least one formula term produced no model-matrix columns.", call. = FALSE)
  }

  list(
    interface = "formula",
    x = design,
    y = as.numeric(response),
    groups = groups,
    group_names = term_labels,
    full_columns = colnames(design),
    response_expression = deparse(formula[[2L]]),
    formula_environment = environment,
    full_formula = formula,
    training_data = training_data
  )
}

.prepare_selection_matrix <- function(x, y) {
  x <- as.matrix(x)
  if (!is.numeric(x) || length(dim(x)) != 2L) {
    stop("`x` must be a numeric matrix or data frame.", call. = FALSE)
  }
  if (!nrow(x) || !ncol(x)) {
    stop("`x` must contain at least one row and one column.", call. = FALSE)
  }
  if (length(y) != nrow(x)) {
    stop("`y` must have one value for every row of `x`.", call. = FALSE)
  }
  complete <- stats::complete.cases(x, y)
  x <- x[complete, , drop = FALSE]
  y <- y[complete]
  .validate_selection_response(y)
  if (any(!is.finite(x))) {
    stop("`x` must contain only finite values after removing incomplete rows.", call. = FALSE)
  }
  column_names <- colnames(x)
  if (is.null(column_names)) {
    column_names <- paste0("x", seq_len(ncol(x)))
  }
  if (anyDuplicated(column_names)) {
    stop("The columns of `x` must have unique names.", call. = FALSE)
  }
  colnames(x) <- column_names

  training_data <- data.frame(.response = as.numeric(y), check.names = FALSE)
  training_data <- cbind(
    training_data,
    as.data.frame(x, check.names = FALSE)
  )
  if (!is.null(rownames(x))) {
    rownames(training_data) <- rownames(x)
  }
  groups <- as.list(seq_len(ncol(x)))
  names(groups) <- column_names

  list(
    interface = "matrix",
    x = x,
    y = as.numeric(y),
    groups = groups,
    group_names = column_names,
    full_columns = column_names,
    response_expression = ".response",
    formula_environment = parent.frame(),
    full_formula = stats::reformulate(column_names, response = ".response"),
    training_data = training_data
  )
}

.validate_selection_response <- function(y) {
  if (!is.numeric(y) || is.matrix(y) || is.factor(y)) {
    stop(
      "`select_nnet()` currently supports only a numeric Gaussian response.",
      call. = FALSE
    )
  }
  if (!length(y) || any(!is.finite(y))) {
    stop("The response must contain finite numeric values.", call. = FALSE)
  }
  if (length(unique(y)) < 2L) {
    stop("The response must contain at least two distinct values.", call. = FALSE)
  }
  invisible(TRUE)
}

.selection_formula <- function(response, terms, environment) {
  formula <- stats::reformulate(terms, response = response)
  environment(formula) <- environment
  formula
}

.validate_selection_arguments <- function(max_size, n_starts, criterion,
                                          start_range, maxit, seed, progress,
                                          nnet_args, n_observations) {
  .selection_positive_integer(max_size, "max_size")
  .selection_positive_integer(n_starts, "n_starts")
  .selection_positive_integer(maxit, "maxit")
  criterion <- match.arg(criterion, c("BIC", "AIC", "AICc"))
  if (!is.numeric(start_range) || length(start_range) != 2L ||
      any(!is.finite(start_range)) || start_range[1L] >= start_range[2L]) {
    stop(
      "`start_range` must contain two finite increasing numeric bounds.",
      call. = FALSE
    )
  }
  if (!is.null(seed) &&
      (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed))) {
    stop("`seed` must be `NULL` or one finite number.", call. = FALSE)
  }
  seed_integer <- if (is.null(seed)) NULL else suppressWarnings(as.integer(seed))
  if (!is.null(seed) && is.na(seed_integer)) {
    stop("`seed` is outside the supported integer range.", call. = FALSE)
  }
  if (!is.logical(progress) || length(progress) != 1L || is.na(progress)) {
    stop("`progress` must be `TRUE` or `FALSE`.", call. = FALSE)
  }
  if (n_observations < 3L) {
    stop("At least three complete observations are required.", call. = FALSE)
  }
  reserved <- c(
    "x", "y", "size", "linout", "entropy", "softmax", "censored",
    "skip", "Wts", "trace", "decay", "Hess", "maxit", "formula", "data",
    "weights", "mask", "summ"
  )
  if (length(nnet_args) &&
      (is.null(names(nnet_args)) || any(!nzchar(names(nnet_args))))) {
    stop("All additional `nnet` arguments must be named.", call. = FALSE)
  }
  supplied <- names(nnet_args)
  supplied <- supplied[!is.na(supplied) & nzchar(supplied)]
  conflicting <- intersect(supplied, reserved)
  if (length(conflicting)) {
    stop(
      sprintf(
        "The following `nnet` arguments are controlled by `select_nnet()`: %s.",
        paste(conflicting, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  list(
    max_size = as.integer(max_size),
    n_starts = as.integer(n_starts),
    criterion = criterion,
    start_range = as.numeric(start_range),
    maxit = as.integer(maxit),
    seed = seed_integer,
    progress = progress
  )
}

.selection_positive_integer <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < 1 || x != as.integer(x)) {
    stop(sprintf("`%s` must be one positive integer.", name), call. = FALSE)
  }
  invisible(TRUE)
}

.select_nnet_impl <- function(prepared, max_size, n_starts, criterion,
                              start_range, maxit, seed, progress, nnet_args,
                              call) {
  control <- .validate_selection_arguments(
    max_size = max_size,
    n_starts = n_starts,
    criterion = criterion,
    start_range = start_range,
    maxit = maxit,
    seed = seed,
    progress = progress,
    nnet_args = nnet_args,
    n_observations = nrow(prepared$x)
  )
  if (!is.null(control$seed)) {
    set.seed(control$seed)
  }

  state <- new.env(parent = emptyenv())
  state$evaluation <- 0L
  state$path <- list()
  state$starts <- list()
  state$failures <- character()

  evaluate <- function(active_groups, size, phase, step, candidate) {
    state$evaluation <- state$evaluation + 1L
    indices <- unlist(prepared$groups[active_groups], use.names = FALSE)
    x_candidate <- prepared$x[, indices, drop = FALSE]
    if (isTRUE(control$progress)) {
      message(sprintf(
        "[%s %d] candidate %d: size = %d, inputs = %s",
        phase,
        step,
        state$evaluation,
        size,
        paste(prepared$group_names[active_groups], collapse = ", ")
      ))
    }
    fitted <- .fit_nnet_multistart(
      x = x_candidate,
      y = prepared$y,
      size = size,
      n_starts = control$n_starts,
      criterion = control$criterion,
      start_range = control$start_range,
      maxit = control$maxit,
      nnet_args = nnet_args
    )
    state$starts[[length(state$starts) + 1L]] <- transform(
      fitted$starts,
      evaluation = state$evaluation,
      phase = phase,
      step = step,
      candidate = candidate
    )
    state$path[[length(state$path) + 1L]] <- data.frame(
      evaluation = state$evaluation,
      phase = phase,
      step = as.integer(step),
      candidate = candidate,
      included_terms = paste(prepared$group_names[active_groups], collapse = " + "),
      n_terms = length(active_groups),
      n_inputs = ncol(x_candidate),
      size = as.integer(size),
      n_parameters = .nnet_parameter_count(ncol(x_candidate), size) + 1L,
      criterion = control$criterion,
      criterion_value = fitted$value,
      converged_starts = sum(fitted$starts$converged),
      n_starts = control$n_starts,
      best_start = fitted$best_start,
      selected = FALSE,
      stringsAsFactors = FALSE
    )
    list(
      active_groups = active_groups,
      size = as.integer(size),
      fit = fitted$fit,
      value = fitted$value,
      evaluation = state$evaluation
    )
  }

  all_groups <- seq_along(prepared$groups)
  hidden_candidates <- lapply(seq_len(control$max_size), function(size) {
    evaluate(all_groups, size, "hidden", 1L, paste0("size_", size))
  })
  current <- .best_selection_candidate(hidden_candidates, "initial hidden-layer search")

  step <- 0L
  repeat {
    step <- step + 1L
    if (length(current$active_groups) <= 1L) break
    removal_candidates <- lapply(current$active_groups, function(group) {
      active <- setdiff(current$active_groups, group)
      evaluate(
        active,
        current$size,
        "input",
        step,
        paste0("remove_", prepared$group_names[group])
      )
    })
    proposed <- .best_selection_candidate(removal_candidates, "input-removal search")
    if (!.criterion_improves(proposed$value, current$value)) break
    current <- proposed
  }

  fine_step <- 0L
  repeat {
    fine_step <- fine_step + 1L
    sizes <- unique(pmax(1L, current$size + c(-1L, 0L, 1L)))
    hidden_neighbours <- lapply(sizes, function(size) {
      evaluate(
        current$active_groups,
        size,
        "fine_hidden",
        fine_step,
        paste0("size_", size)
      )
    })
    proposed_hidden <- .best_selection_candidate(
      hidden_neighbours,
      "fine-tuning hidden-layer search"
    )
    if (proposed_hidden$size == current$size ||
        !.criterion_improves(proposed_hidden$value, current$value)) {
      break
    }
    current <- proposed_hidden

    moves <- list()
    if (length(current$active_groups) > 1L) {
      moves <- c(moves, lapply(current$active_groups, function(group) {
        evaluate(
          setdiff(current$active_groups, group),
          current$size,
          "fine_input",
          fine_step,
          paste0("remove_", prepared$group_names[group])
        )
      }))
    }
    excluded <- setdiff(all_groups, current$active_groups)
    if (length(excluded)) {
      moves <- c(moves, lapply(excluded, function(group) {
        evaluate(
          sort(c(current$active_groups, group)),
          current$size,
          "fine_input",
          fine_step,
          paste0("add_", prepared$group_names[group])
        )
      }))
    }
    if (!length(moves)) break
    proposed_input <- .best_selection_candidate(moves, "fine-tuning input search")
    if (!.criterion_improves(proposed_input$value, current$value)) break
    current <- proposed_input
  }

  if (is.null(current$fit) || !is.finite(current$value)) {
    stop("No candidate model converged successfully.", call. = FALSE)
  }
  path <- do.call(rbind, state$path)
  path$selected[path$evaluation == current$evaluation] <- TRUE
  rownames(path) <- NULL
  starts <- do.call(rbind, state$starts)
  rownames(starts) <- NULL

  selected_groups <- current$active_groups
  selected_names <- prepared$group_names[selected_groups]
  excluded_names <- prepared$group_names[setdiff(all_groups, selected_groups)]
  selected_columns <- unlist(prepared$groups[selected_groups], use.names = FALSE)
  selected_formula <- .selection_formula(
    prepared$response_expression,
    selected_names,
    prepared$formula_environment
  )
  selected_terms_object <- stats::terms(
    selected_formula,
    data = prepared$training_data
  )
  selected_frame <- stats::model.frame(
    selected_terms_object,
    data = prepared$training_data
  )
  selected_terms_object <- attr(selected_frame, "terms")
  selected_matrix <- stats::model.matrix(
    selected_terms_object,
    data = selected_frame
  )

  result <- list(
    fit = current$fit,
    call = call,
    interface = prepared$interface,
    formula = selected_formula,
    full_formula = prepared$full_formula,
    terms = selected_terms_object,
    contrasts = attr(selected_matrix, "contrasts"),
    xlevels = .model_xlevels(selected_frame),
    training_data = prepared$training_data,
    selected_terms = selected_names,
    excluded_terms = excluded_names,
    selected_columns = prepared$full_columns[selected_columns],
    full_columns = prepared$full_columns,
    size = current$size,
    max_size = control$max_size,
    criterion = control$criterion,
    criterion_value = current$value,
    search_path = path,
    start_diagnostics = starts,
    n_starts = control$n_starts,
    start_range = control$start_range,
    maxit = control$maxit,
    seed = control$seed,
    response = "continuous",
    diagnostics = list(
      n_candidates = nrow(path),
      n_failed_candidates = sum(path$converged_starts == 0L),
      n_converged_starts = sum(starts$converged),
      n_failed_starts = sum(!starts$converged),
      selected_evaluation = current$evaluation
    )
  )
  class(result) <- "statnnet_selection"
  result
}

.fit_nnet_multistart <- function(x, y, size, n_starts, criterion,
                                 start_range, maxit, nnet_args) {
  n_weights <- .nnet_parameter_count(ncol(x), size)
  starts <- vector("list", n_starts)
  fits <- vector("list", n_starts)
  values <- rep(Inf, n_starts)

  for (start in seq_len(n_starts)) {
    initial_weights <- stats::runif(
      n_weights,
      min = start_range[1L],
      max = start_range[2L]
    )
    warnings <- character()
    fit <- tryCatch(
      withCallingHandlers(
        do.call(
          nnet::nnet,
          c(
            list(
              x = x,
              y = y,
              size = as.integer(size),
              linout = TRUE,
              entropy = FALSE,
              softmax = FALSE,
              censored = FALSE,
              skip = FALSE,
              Wts = initial_weights,
              trace = FALSE,
              decay = 0,
              Hess = FALSE,
              maxit = as.integer(maxit)
            ),
            nnet_args
          )
        ),
        warning = function(warning) {
          warnings <<- c(warnings, conditionMessage(warning))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(error) error
    )
    error_message <- NA_character_
    converged <- FALSE
    rss <- NA_real_
    log_likelihood <- NA_real_
    value <- Inf
    iterations <- NA_integer_

    if (inherits(fit, "error")) {
      error_message <- conditionMessage(fit)
      fit <- NULL
    } else {
      converged <- identical(as.integer(fit$convergence), 0L) &&
        all(is.finite(fit$wts)) && all(is.finite(fit$fitted.values))
      iterations <- if (!is.null(fit$counts)) as.integer(fit$counts[1L]) else NA_integer_
      if (converged) {
        residual <- y - as.numeric(fit$fitted.values)
        rss <- sum(residual^2)
        sigma2 <- rss / length(y)
        if (is.finite(sigma2) && sigma2 > .Machine$double.eps) {
          log_likelihood <- -length(y) / 2 * (log(2 * pi * sigma2) + 1)
          value <- .nnet_information_criterion(
            log_likelihood = log_likelihood,
            n_parameters = n_weights + 1L,
            n_observations = length(y),
            criterion = criterion
          )
        } else {
          error_message <- "residual variance is numerically zero"
          converged <- FALSE
        }
      } else {
        error_message <- sprintf(
          "nnet convergence code %s",
          if (is.null(fit$convergence)) "unavailable" else fit$convergence
        )
      }
    }
    fits[[start]] <- fit
    values[start] <- value
    starts[[start]] <- data.frame(
      start = start,
      converged = converged,
      convergence_code = if (is.null(fit)) NA_integer_ else as.integer(fit$convergence),
      iterations = iterations,
      rss = rss,
      log_likelihood = log_likelihood,
      criterion_value = value,
      warning = if (length(warnings)) paste(unique(warnings), collapse = "; ") else NA_character_,
      error = error_message,
      stringsAsFactors = FALSE
    )
  }
  starts <- do.call(rbind, starts)
  best_start <- if (any(is.finite(values))) which.min(values) else NA_integer_
  list(
    fit = if (is.na(best_start)) NULL else fits[[best_start]],
    value = if (is.na(best_start)) Inf else values[best_start],
    best_start = best_start,
    starts = starts
  )
}

.nnet_parameter_count <- function(n_inputs, size) {
  as.integer((n_inputs + 1L) * size + size + 1L)
}

.nnet_information_criterion <- function(log_likelihood, n_parameters,
                                        n_observations, criterion) {
  if (criterion == "BIC") {
    return(-2 * log_likelihood + log(n_observations) * n_parameters)
  }
  if (criterion == "AIC") {
    return(-2 * log_likelihood + 2 * n_parameters)
  }
  denominator <- n_observations - n_parameters - 1L
  if (denominator <= 0L) return(Inf)
  -2 * log_likelihood + 2 * n_parameters * n_observations / denominator
}

.best_selection_candidate <- function(candidates, context) {
  values <- vapply(candidates, `[[`, numeric(1), "value")
  if (!any(is.finite(values))) {
    stop(sprintf("No candidate converged during the %s.", context), call. = FALSE)
  }
  candidates[[which.min(values)]]
}

.criterion_improves <- function(candidate, current) {
  is.finite(candidate) && candidate < current - sqrt(.Machine$double.eps) * max(1, abs(current))
}
