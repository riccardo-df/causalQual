#' Summary Method for causalQual Objects
#' 
#' Summarizes an \code{\link{causalQual}} object.
#' 
#' @param object An \code{\link{causalQual}} object.
#' @param ... Further arguments passed to or from other methods.
#' 
#' @return 
#' Summarizes an \code{\link{causalQual}} object.
#' 
#' @examples 
#' \donttest{## Generate synthetic data.
#' set.seed(1986)
#' 
#' data <- generate_qualitative_data_soo(1000, assignment = "observational")
#' Y <- data$Y
#' D <- data$D
#' X <- data$X
#' 
#' ## Estimate probability shifts.
#' fit <- causalQual(Y, D, X, 
#'                   identification = "selection_on_observables",
#'                   outcome_type = "ordered")
#' summary(fit)}
#' 
#' @import cli
#' 
#' @author Riccardo Di Francesco
#'
#' @seealso \code{\link{causalQual}}
#' 
#' @export
summary.causalQual <- function(object, ...) {
  cli::cli_h1("CAUSAL INFERENCE WITH QUALITATIVE OUTCOMES")
  cli::cli_h2("Research design")
  cat("Identification:         ", object$identification, "\n")
  cat("Outcome type:           ", object$outcome_type, "\n")
  cat("Full sample size:       ", dim(object$data)[1], "\n")
  cat("N. covariates:          ", dim(object$data)[2] - 2, "\n")
  cat("Classes:                ", sort(unique(object$data$Y)), "\n")
  cat("Fraction treated units: ", mean(object$data$D), "\n\n")
  cli::cli_h2("Point estimates and 95\\% confidence intervals")
  estimates <- object$estimates
  estimates <- format(round(estimates, 3), nsmall = 3)
  
  ci_lower <- format(round(object$ci_lower, 3), nsmall = 3)
  ci_upper <- format(round(object$ci_upper, 3), nsmall = 3)
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
#' data <- generate_qualitative_data_soo(1000, assignment = "observational")
#' Y <- data$Y
#' D <- data$D
#' X <- data$X
#' 
#' ## Estimate probability shifts.
#' fit <- causalQual(Y, D, X, 
#'                   identification = "selection_on_observables",
#'                   outcome_type = "ordered")
#' print(fit)}
#' 
#' @author Riccardo Di Francesco
#'
#' @seealso \code{\link{causalQual}}
#' 
#' @export
print.causalQual <- function(x, ...) {
  summary.causalQual(x, ...)
}


#' Plot Method for causalQual Objects
#'
#' Plots an \code{\link{causalQual}} object.
#'
#' @param x An \code{\link{causalQual}} object.
#' @param hline Logical, whether to display an horizontal line at zero in the plot.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return
#' Plots an \code{\link{causalQual}} object.
#'
#' \donttest{## Generate synthetic data.
#' set.seed(1986)
#' 
#' data <- generate_qualitative_data_soo(1000, assignment = "observational")
#' Y <- data$Y
#' D <- data$D
#' X <- data$X
#' 
#' ## Estimate probability shifts.
#' fit <- causalQual(Y, D, X, 
#'                   identification = "selection_on_observables",
#'                   outcome_type = "ordered")
#' plot(fit)}
#' 
#' @import ggplot2 ggsci
#' @importFrom stats coef
#'
#' @author Riccardo Di Francesco
#'
#' @seealso \code{\link{causalQual}}
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
                          ci_lower = x$ci_lower,
                          ci_upper = x$ci_upper)
  
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