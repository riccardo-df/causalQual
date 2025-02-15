#' Summary Method for causalQual Objects
#'
#' Summarizes an \code{causalQual} object.
#'
#' @param object An \code{causalQual} object.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return
#' Summarizes an \code{causalQual} object.
#'
#' @examples
#' \donttest{## Generate synthetic data.
#' set.seed(1986)
#'
#' data <- generate_qualitative_data_soo(1000, assignment = "observational",
#'                                       outcome_type = "ordered")
#' Y <- data$Y
#' D <- data$D
#' X <- data$X
#'
#' ## Estimate probabilities of shifts.
#' fit <- causalQual_soo(Y = Y, D = D, X = X, outcome_type = "ordered")
#' summary(fit)}
#'
#' @import cli
#'
#' @author Riccardo Di Francesco
#'
#' @seealso causalQual
#'
#' @export
summary.causalQual <- function(object, ...) {
  if (object$identification == "selection_on_observables") {
    identification <- "Selection-on-Observables"
    estimand <- "Probability shifts"
  } else if (object$identification == "iv") {
    identification <- "Instrumental Variables"
    estimand <- "Local Probability Shifts"
  } else if (object$identification == "rd") {
    identification <- "Regression Discontinuity"
    estimand <- "Probability Shifts at the Cutoff"
  } else if (object$identification == "diff_in_diff") {
    identification <- "Difference-in-Differences"
    estimand <- "Probability Shifts on the Treated"
  }

  cli::cli_h1("CAUSAL INFERENCE WITH QUALITATIVE OUTCOMES")
  cli::cli_h2("Research design")
  cat("Identification:         ", identification, "\n")
  cat("Estimand:               ", estimand, "\n")
  cat("Outcome type:           ", object$outcome_type, "\n")
  cat("Classes:                ", if (object$identification == "diff_in_diff") sort(unique(object$data$Y_pre)) else sort(unique(object$data$Y)), "\n")
  cat("N. units:               ", length(object$data$D), "\n")
  cat("Fraction treated units: ", mean(object$data$D), "\n\n")

  cli::cli_h2("Point estimates and 95\\% confidence intervals")
  estimates <- object$estimates
  estimates <- format(round(estimates, 3), nsmall = 3)

  ci_lower <- object$confidence_intervals$lower
  ci_upper <- object$confidence_intervals$upper

  ci_lower <- format(round(ci_lower, 3), nsmall = 3)
  ci_upper <- format(round(ci_upper, 3), nsmall = 3)

  formatted_cis <- paste0("[", ci_lower, ", ", ci_upper, "]")

  for (i in seq_along(estimates)) {
    cat(paste0("Class ", i, ": ", estimates[i], "  ", formatted_cis[i], "\n"))
  }
}


#' Print Method for causalQual Objects
#'
#' Prints an \code{causalQual} object.
#'
#' @param x An \code{causalQual} object.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return
#' Prints an \code{causalQual} object.
#'
#' @examples
#' \donttest{## Generate synthetic data.
#' set.seed(1986)
#'
#' data <- generate_qualitative_data_soo(1000, assignment = "observational",
#'                                       outcome_type = "ordered")
#' Y <- data$Y
#' D <- data$D
#' X <- data$X
#'
#' ## Estimate probabilities of shifts.
#' fit <- causalQual_soo(Y = Y, D = D, X = X, outcome_type = "ordered")
#' print(fit)}
#'
#' @author Riccardo Di Francesco
#'
#' @seealso causalQual
#'
#' @export
print.causalQual <- function(x, ...) {
  summary.causalQual(x, ...)
}


#' Plot Method for causalQual Objects
#'
#' Plots an \code{causalQual} object.
#'
#' @param x An \code{causalQual} object.
#' @param hline Logical, whether to display an horizontal line at zero in the plot.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return
#' Plots an causalQual object.
#'
#' @examples
#' \donttest{## Generate synthetic data.
#' set.seed(1986)
#'
#' data <- generate_qualitative_data_soo(1000, assignment = "observational",
#'                                       outcome_type = "ordered")
#' Y <- data$Y
#' D <- data$D
#' X <- data$X
#'
#' ## Estimate probabilities of shifts.
#' fit <- causalQual_soo(Y = Y, D = D, X = X, outcome_type = "ordered")
#' plot(fit)}
#'
#' @import ggplot2 ggsci
#' @importFrom stats coef
#'
#' @author Riccardo Di Francesco
#'
#' @seealso causalQual
#'
#' @export
plot.causalQual <- function(x, hline = TRUE, ...) {
  ## 0.) Handle input and checks.
  if (!is.logical(hline)) stop("Invalid 'hline'. Must be logical.", call. = FALSE)

  estimate <- NULL
  ci_lower <- NULL
  ci_upper <- NULL

  ## 1.) Prepare data for plot.
  plot_data <- data.frame(class = factor(seq_along(x$estimates), labels = paste0("Class ", seq_along(x$estimates))),
                          estimate = x$estimates,
                          ci_lower = x$estimates - 1.96 * x$standard_errors,
                          ci_upper = x$estimates + 1.96 * x$standard_errors)

  ## 2.) Construct plot.
  plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = class, y = estimate)) +
    ggplot2::geom_point(size = 2, color = "#1F78B4") +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, color = "#1F78B4") +
    ggplot2::geom_hline(yintercept = if (hline) 0 else NULL, linetype = "dashed") +
    ggplot2::labs(x = "", y = "") +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(axis.text = ggplot2::element_text(size = 12))

  ## 3.) Plot
  plot
}
