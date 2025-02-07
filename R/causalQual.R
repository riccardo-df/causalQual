#' Causal Inference with Qualitative Outcomes
#'
#' Estimates treatment effects on qualitative outcomes (e.g., multinomial or ordered non-numeric) under several research designs.
#'
#' @param Y Ordered non-numeric outcome. Must be labeled as \eqn{\{1, 2, \dots\}}.
#' @param D Binary treatment indicator.
#' @param X Covariate matrix (no intercept).
#' @param identification String controlling the employed identification strategy. Must be one of \code{"selection_on_observables"}.
#' @param outcome_type String controlling the outcome type. Must be either \code{"multinomial"} or \code{"ordered"}. Affects estimation of conditional class probabilities.
#' @param K Number of folds for nuisance functions estimation. Ignored if \code{identification != "selection_on_observables"}.
#'
#' @return An object of class \code{\link{causalQual}}.
#'
#' @examples
#' \donttest{## Generate synthetic data.
#' set.seed(1986)
#'
#' data <- generate_qualitative_data_soo(100, assignment = "observational")
#' Y <- data$Y
#' D <- data$D
#' X <- data$X
#'
#' ## Estimate probability shifts.
#' fit <- causalQual(Y, D, X,
#'                   identification = "selection_on_observables",
#'                   outcome_type = "ordered")
#'
#' summary(fit)
#' plot(fit)}
#'
#' @details
#' Depending on \code{identification}, \code{\link{causalQual}} calls the appropriate estimation function.
#'
#' ## Selection-on-observables
#' Under selection-on-observables, identification requires the treatment indicator to be (conditionally) independent of potential outcomes
#' (unconfoundedness), and that each unit has a non-zero probability of being treated (common support).\cr
#'
#' If \code{identification == "selection_on_observables"}, \code{\link{causalQual}} constructs and averages doubly robust scores for qualitative outcomes
#' to estimate the probability shifts \eqn{\delta_m := P(Y(1) = m) - P(Y(0) = m)}. For each class \eqn{m}, the doubly robust score is defined as:
#' \deqn{
#'     \hat{\Gamma}_{m, i} =
#'     \hat{P}(Y = m \mid D = 1, X_i) - \hat{P}(Y = m \mid D = 0, X_i) +
#'     D_i \frac{1\{Y = m\} - \hat{P}(Y = m \mid D = 1, X_i)}{\hat{P}(D_i = 1 | X_i)}
#'     - (1 - D_i) \frac{1\{Y = m\} - \hat{P}(Y = m \mid D = 0, X_i)}{1 - \hat{P}(D_i = 1 | X_i)}.
#' }
#'
#' The estimator for the probability shift is then the average of the scores:
#' \deqn{
#'     \hat{\delta}_m = \frac{1}{n} \sum_{i=1}^{n} \hat{\Gamma}_{m, i}.
#' }
#'
#' The variance of \eqn{\hat{\delta}_m} is estimated as:
#' \deqn{
#'     \widehat{\text{Var}} ( \hat{\delta}_m ) = \frac{1}{n} \sum_{i=1}^{n} ( \hat{\Gamma}_{m, i} - \hat{\delta}_m )^2.
#' }
#'
#' If \code{outcome_type == "multinomial"}, \eqn{\hat{P}(Y = m \mid D = 1, X)} and \eqn{\hat{P}(Y = m \mid D = 0, X)} are estimated using a \code{\link[ocf]{multinomial_ml}} strategy with regression forests
#' as base learner. Else, if \code{outcome_type == "ordered"}, \eqn{\hat{P}(Y = m \mid D = 1, X)} and \eqn{\hat{P}(Y = m \mid D = 0, X)} are estimated using the honest version of the \code{\link[ocf]{ocf}} estimator.
#' \eqn{\hat{P}(D_i = 1 | X_i)} is always estimated via a honest \code{\link[grf]{regression_forest}}. K-fold cross-fitting is employed.
#'
#' @seealso \code{\link{causalQual_soo}}
#'
#' @author Riccardo Di Francesco
#'
#' @export
causalQual <- function(Y, D, X, identification, outcome_type, K = 5) {
  ## 0.) Handling inputs and checks.
  if (!(identification %in% c("selection_on_observables"))) stop("Invalid 'identification'. Must be one of 'selection_on_observables'.", call. = FALSE)
  if (!(outcome_type %in% c("multinomial", "ordered"))) stop("Invalid 'outcome_type'. Must be either 'multinomial' or 'ordered'.", call. = FALSE)

  ## 1.) Call appropriate function.
  if (identification == "selection_on_observables") {
    fit <- causalQual_soo(Y, D, X, outcome_type, K)
  }

  ## 2.) Output.
  output <- list()
  output$identification <- identification
  output$outcome_type <- outcome_type
  output$estimates <- fit$estimates
  output$standard_errors <- fit$standard_errors
  output$ci_lower <- fit$ci_lower
  output$ci_upper <- fit$ci_upper
  output$data <- data.frame("Y" = Y, "D" = D, X)

  class(output) <- "causalQual"
  return(output)
}
