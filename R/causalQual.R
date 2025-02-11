#' Causal Inference with Qualitative Outcomes
#'
#' Estimates treatment effects on qualitative outcomes (e.g., multinomial or ordered non-numeric) under several research designs.
#'
#' @param Y Qualitative outcome. Must be labeled as \eqn{\{1, 2, \dots\}}. Ignored if \code{outcome_type == "diff_in_diff"}.
#' @param Y_pre Qualitative outcome before treatment. Must be labeled as \eqn{\{1, 2, \dots\}}. Ignored if \code{outcome_type != "diff_in_diff"}.
#' @param Y_post Qualitative outcome after treatment. Must be labeled as \eqn{\{1, 2, \dots\}}. Ignored if \code{outcome_type != "diff_in_diff"}.
#' @param D Binary treatment indicator.
#' @param X Covariate matrix (no intercept). Ignored if \code{outcome_type != "selection_on_observables"}.
#' @param Z Binary instrument. Ignored if \code{identification != "iv"}.
#' @param identification String controlling the employed identification strategy. Must be one of \code{"selection_on_observables"}, \code{"iv"}, \code{"diff_in_diff"}.
#' @param outcome_type String controlling the outcome type. Must be either \code{"multinomial"} or \code{"ordered"}. Affects estimation of conditional class probabilities. Ignored if \code{identification != "selection_on_observables"}.
#' @param K Number of folds for nuisance functions estimation. Ignored if \code{identification != "selection_on_observables"}.
#'
#' @return An object of class \code{\link{causalQual}}.
#'
#' @examples
#' \donttest{## Selection-on-observables.
#' # Generate synthetic data.
#' set.seed(1986)
#'
#' data <- generate_qualitative_data_soo(1000, assignment = "observational",
#'                                       outcome_type = "ordered")
#' Y <- data$Y
#' D <- data$D
#' X <- data$X
#'
#' # Estimate probabilities of shift.
#' fit <- causalQual(Y = Y, D = D, X = X,
#'                   identification = "selection_on_observables",
#'                   outcome_type = "ordered", K = 2)
#'
#' summary(fit)
#' plot(fit)
#'
#' ## Instrumental variables.
#' # Generate synthetic data.
#' set.seed(1986)
#'
#' data <- generate_qualitative_data_iv(1000,
#'                                      outcome_type = "ordered")
#' Y <- data$Y
#' D <- data$D
#' Z <- data$Z
#'
#' # Estimate local probabilities of shift.
#' fit <- causalQual(Y = Y, D = D, Z = Z,
#'                   identification = "iv",
#'                   outcome_type = "ordered")
#'
#' summary(fit)
#' plot(fit)
#'
#' ## Difference-in-differences.
#' # Generate synthetic data.
#' set.seed(1986)
#'
#' data <- generate_qualitative_data_did(1000, assignment = "observational",
#'                                       outcome_type = "ordered")
#' Y_pre <- data$Y_pre
#' Y_post <- data$Y_post
#' D <- data$D
#'
#' # Estimate probabilities of shift on the treated.
#' fit <- causalQual(Y_pre = Y_pre, Y_post = Y_post, D = D,
#'                   identification = "diff_in_diff",
#'                   outcome_type = "ordered")
#'
#' summary(fit)
#' plot(fit)}
#'
#' @details
#' ## Selection-on-observables
#' Under a selection-on-observables design, identification requires the treatment indicator to be (conditionally) independent of potential outcomes
#' (unconfoundedness), and that each unit has a non-zero probability of being treated (common support). If these hold, we can recover the probabilities of shift of all classes.\cr
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
#' as base learners. Else, if \code{outcome_type == "ordered"}, \eqn{\hat{P}(Y = m \mid D = 1, X)} and \eqn{\hat{P}(Y = m \mid D = 0, X)} are estimated using the honest version of the \code{\link[ocf]{ocf}} estimator.
#' \eqn{\hat{P}(D_i = 1 | X_i)} is always estimated via a honest \code{\link[grf]{regression_forest}}. K-fold cross-fitting is employed.
#'
#' ## Instrumental-variables
#' Under an instrumental-variables design, identification requires the treatment assignment (i.e., the instrument) to be independent of potential outcomes and potential treatments (exogeneity), that the
#' instrument influences the outcome solely through its effect on treatment (exclusion restriction), that the instrument has a nonzero effect on treatment probability (relevance), and that the instrument can only
#' increase/decrease the treatment probability (monotonicity). If these hold, we can recover the local probabilities of shift of all classes.\cr
#'
#' If \code{identification == "iv"}, \code{\link{causalQual}} fits one standard two-stage least squares model for each class \eqn{m} of \code{Y}. Each model uses a binary indicator \eqn{1(Y_{it} = m)} as an outcome.
#'
#' ## Difference-in-differences
#' Under a difference-in-difference design, identification requires that the probabilities time shift for class \eqn{m} evolve similarly for the treated and control groups  (parallel
#' trends on the probability mass functions of \code{Y}). If this holds, we can recover the probability of shift on the treated for class \eqn{m}.\cr
#'
#' If \code{identification == "selection_on_observables"}, \code{\link{causalQual}} fits one linear model for each class \eqn{m} of \code{Y}. Each model uses a binary indicator \eqn{1(Y_{it} = m)} as an outcome and a
#' time dummy, a treatment dummy (\code{D}), and an interaction thereof as regressors. The coefficient on the interaction identifies the probability of shift on the treated for class \eqn{m}. Standard errors are clustered
#' at the unit level.
#'
#' @author Riccardo Di Francesco
#'
#' @export
causalQual <- function(Y = NULL, Y_pre = NULL, Y_post = NULL, D = NULL, X = NULL, Z = NULL, identification = NULL, outcome_type = NULL, K = 5) {
  ## 0.) Handling inputs and checks.
  if (any(!(D %in% c(0, 1)))) stop("Invalid 'D'. Treatment must be binary {0, 1}.", call. = FALSE)
  if (!(identification %in% c("selection_on_observables", "iv", "diff_in_diff"))) stop("Invalid 'identification'. Must be one of 'selection_on_observables', 'iv', 'diff_in_diff'.", call. = FALSE)
  if (identification == "diff_in_diff" & !is.null(Y)) stop("When 'identification' is 'diff_in_diff', we do not need 'Y'. Rather, please provide 'Y_pre' and 'Y_post'.", call. = FALSE)
  if (identification != "diff_in_diff" & (!is.null(Y_pre) | !is.null(Y_post))) stop("When 'identification' is 'selection_on_observables', we do not need 'Y_pre' and 'Y_post'. Rather, please provide 'Y'.", call. = FALSE)
  if (identification != "selection_on_observables" & !is.null(X)) stop("When 'identification' is 'selection_on_observables', we do not need 'X'.", call. = FALSE)

  ## 1.) Call appropriate function.
  if (identification == "selection_on_observables") {
    fit <- causalQual_soo(Y, D, X, outcome_type, K)
  } else if (identification == "iv") {
    fit <- causalQual_iv(Y, D, Z)
  } else if (identification == "diff_in_diff") {
    fit <- causalQual_did(Y_pre, Y_post, D)
  }

  ## 2.) Output.
  output <- list()
  output$identification <- identification
  output$outcome_type <- outcome_type
  output$estimates <- fit$estimates
  output$standard_errors <- fit$standard_errors
  output$dr_scores <- if (identification == "selection_on_observables") fit$dr_scores else NULL
  output$fits <- if (identification != "diff_in_diff") fit$fits else NULL
  output$data <- list("Y" = Y, "Y_pre" = Y_pre, "Y_post" = Y_post, "D" = D, "X" = X, "Z" = Z)

  class(output) <- "causalQual"
  return(output)
}
