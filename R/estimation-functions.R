#' Causal Inference for Qualitative Outcomes under Selection-on-Observables
#'
#' Construct and average doubly robust scores for qualitative outcomes to estimate the probabilities of shift.
#'
#' @param Y Qualitative outcome. Must be labeled as \eqn{\{1, 2, \dots\}}.
#' @param D Binary treatment indicator.
#' @param X Covariate matrix (no intercept).
#' @param outcome_type String controlling the outcome type. Must be either \code{"multinomial"} or \code{"ordered"}. Affects estimation of conditional class probabilities.
#' @param K Number of folds for nuisance functions estimation.
#'
#' @return An object of class \code{causalQual}.
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
#' # Estimate probabilities of shift.
#' fit <- causalQual_soo(Y, D, X, outcome_type = "ordered", K = 2)
#'
#' summary(fit)
#' plot(fit)}
#'
#' @details
#' Under a selection-on-observables design, identification requires the treatment indicator to be (conditionally) independent of potential outcomes
#' (unconfoundedness), and that each unit has a non-zero probability of being treated (common support). If these assumptions hold, we can recover the
#' probabilities of shift of all classes:
#'
#' \deqn{\delta_m := P(Y_i(1) = m) - P(Y_i(0) = m), \, m = 1, \dots, M.}
#'
#' \code{\link{causalQual_soo}} constructs and averages doubly robust scores for qualitative outcomes
#' to estimate \eqn{\delta_m}. For each class \eqn{m}, the doubly robust score for unit \eqn{i} is defined as:
#'
#' \deqn{
#'     \hat{\Gamma}_{m, i} =
#'     \hat{P}(Y_i = m \mid D_i = 1, X_i) -
#'     \hat{P}(Y_i = m \mid D_i = 0, X_i) +
#' }
#' \deqn{
#'     D_i \frac{1\{Y_i = m\} - \hat{P}(Y_i = m \mid D_i = 1, X_i)}
#'     {\hat{P}(D_i = 1 | X_i)}
#'     - (1 - D_i) \frac{1\{Y_i = m\} - \hat{P}(Y_i = m \mid D_i = 0, X_i)}
#'     {1 - \hat{P}(D_i = 1 | X_i)}.
#' }
#'
#' The estimator for \eqn{\delta_m} is then the average of the scores:
#'
#' \deqn{\hat{\delta}_m = \frac{1}{n} \sum_{i=1}^{n} \hat{\Gamma}_{m, i},}
#'
#' with its variance estimated as:
#' \deqn{
#'     \widehat{\text{Var}} ( \hat{\delta}_m ) = \frac{1}{n} \sum_{i=1}^{n} ( \hat{\Gamma}_{m, i} - \hat{\delta}_m )^2.
#' }
#'
#' \code{\link{causalQual_soo}} uses these estimates to construct confidence intervals based on conventional normal approximations.\cr
#'
#' If \code{outcome_type == "multinomial"}, \eqn{\hat{P}(Y_i = m \mid D_i = 1, X_i)} and \eqn{\hat{P}(Y_i = m \mid D_i = 0, X_i)} are estimated using a \code{\link[ocf]{multinomial_ml}} strategy with regression forests
#' as base learners. Else, if \code{outcome_type == "ordered"}, \eqn{\hat{P}(Y_i = m \mid D_i = 1, X_i)} and \eqn{\hat{P}(Y_i = m \mid D_i = 0, X_i)} are estimated using the honest version of the \code{\link[ocf]{ocf}} estimator.
#' \eqn{\hat{P}(D_i = 1 | X_i)} is always estimated via a honest \code{\link[grf]{regression_forest}}. \code{K}-fold cross-fitting is employed for the estimation of all these functions.\cr
#'
#' Folds are created by random split. If some class of \code{Y} is not observed in one or more folds for one or both treatment groups, a new random partition is performed. This process is repeat until when all
#' classes are observed in all folds and for all treatment groups up to 1000 times, after which the routine raises an error.
#'
#' @import ocf grf
#' @importFrom caret createFolds
#' @importFrom magrittr %>%
#'
#' @author Riccardo Di Francesco
#'
#' @seealso \code{\link{causalQual_iv}} \code{\link{causalQual_rd}} \code{\link{causalQual_did}}
#'
#' @export
causalQual_soo <- function(Y, D, X, outcome_type, K = 5) {
  ## 0.) Handling inputs and checks.
  if (K < 1 | K %% 1 != 0) stop("Invalid 'K'. Number of folds must be a positive integer.", call. = FALSE)
  if (!(outcome_type %in% c("multinomial", "ordered"))) stop("Invalid 'outcome_type'. Must be either 'multinomial' or 'ordered'.", call. = FALSE)
  if (!all(sort(unique(Y)) == seq_len(max(Y)))) stop("Invalid 'Y'. Outcome must be coded with consecutive integers from 1 to max(Y).", call. = FALSE)

  classes <- sort(unique(Y))
  treated_idx <- which(D == 1)
  control_idx <- which(D == 0)

  ## 1.) Generate binary outcomes for each class. The m-th element of the list stores the indicator variables relative to the m-th class.
  indicators <- list()
  counter <- 1
  for (m in classes) {
    indicators[[counter]] <- ifelse(Y == m, 1, 0)
    counter <- counter + 1
  }
  names(indicators) <- paste0("Class", classes)

  ## 2.) Estimate conditional class probabilities and propensity score via cross-fitting. Repeat fold splits until all values of outcome are observed in all folds and groups. Stop after 1000 trials.
  folds <- caret::createFolds(seq_len(length(Y)), K, list = TRUE)

  counter <- 1
  while (any(sapply(folds, function(f) { as.numeric(table(Y[f], D[f])) } ) == 0)) {
    folds <- caret::createFolds(seq_len(length(Y)), K, list = TRUE)
    if (counter == 1000) stop("We cannot construct 'K' folds such that all values of 'Y' are observed in all folds. \nMaybe try with a smaller 'K'?", call. = FALSE)
    counter <- counter + 1
  }

  pclasses1 <- matrix(NA, nrow = length(Y), ncol = max(Y))
  pclasses0 <- matrix(NA, nrow = length(Y), ncol = max(Y))
  pscore_hat <- numeric(length(Y))

  for (k in seq_len(K)) {
    ## Select holdout fold and its complement.
    holdout_idx <- folds[[k]]
    train_idx <- setdiff(seq_len(length(Y)), holdout_idx)

    ## Identify treated and control groups in holdout fold and its complement.
    treated_idx_train <- intersect(train_idx, which(D == 1))
    control_idx_train <- intersect(train_idx, which(D == 0))

    ## Fit multinomial ML/OCF separately for treated and control groups for conditional class probabilities. Fit regression forest for propensity score.
    if (outcome_type == "multinomial") {
      forests1 <- ocf::multinomial_ml(Y[treated_idx_train], X[treated_idx_train, ], learner = "forest")
      forests0 <- ocf::multinomial_ml(Y[control_idx_train], X[control_idx_train, ], learner = "forest")
    } else if (outcome_type == "ordered") {
      forests1 <- ocf::ocf(Y[treated_idx_train], X[treated_idx_train, ], honesty = TRUE, alpha = 0.2)
      forests0 <- ocf::ocf(Y[control_idx_train], X[control_idx_train, ], honesty = TRUE, alpha = 0.2)
    }

    forest_pscore <- grf::regression_forest(as.matrix(X[train_idx, ]), D[train_idx])

    ## Predict on holdout.
    if (outcome_type == "multinomial") {
      pclasses1[holdout_idx, ] <- stats::predict(forests1, X[holdout_idx, ])
      pclasses0[holdout_idx, ] <- stats::predict(forests0, X[holdout_idx, ])
    } else if (outcome_type == "ordered") {
      pclasses1[holdout_idx, ] <- stats::predict(forests1, X[holdout_idx, ])$probabilities
      pclasses0[holdout_idx, ] <- stats::predict(forests0, X[holdout_idx, ])$probabilities
    }

    pscore_hat[holdout_idx] <- stats::predict(forest_pscore, as.matrix(X[holdout_idx, ]))$predictions
  }

  ## 3.) Construct doubly-robust scores.
  gammas <- list()
  counter <- 1
  for (m in sort(unique(Y))) {
    gammas[[counter]] <- pclasses1[, m] - pclasses0[, m] + D * (indicators[[m]] - pclasses1[, m]) / pscore_hat - (1 - D) * (indicators[[m]] - pclasses0[, m]) / (1 - pscore_hat)
    counter <- counter + 1
  }
  names(gammas) <- paste0("Class", sort(unique(Y)))

  ## 4.) Pick mean and standard deviation of doubly-robust scores. Construct conventional confidence intervals.
  pshifts_hat <- sapply(gammas, mean)
  pshifts_hat_se <- sapply(seq_along(gammas), function(m) mean((gammas[[m]] - pshifts_hat[m])^2) / length(Y)) %>%
    sqrt()
  cis_lower <- pshifts_hat - 1.96 * pshifts_hat_se
  cis_upper <- pshifts_hat + 1.96 * pshifts_hat_se

  ## 5.) Output.
  output <- list()
  output$identification <- "selection_on_observables"
  output$outcome_type <- outcome_type
  output$estimates <- pshifts_hat
  output$standard_errors <- pshifts_hat_se
  output$confidence_intervals <- list("lower" = cis_lower, "upper" = cis_upper)
  output$dr_scores <- gammas
  output$fits <- NULL
  output$data <- list("Y" = Y, "Y_pre" = NULL, "Y_post" = NULL, "D" = D, "X" = X, "Z" = NULL, "running_variable" = NULL, "cutoff" = NULL)

  class(output) <- "causalQual"
  return(output)
}


#' Causal Inference for Qualitative Outcomes under Instrumental Variables
#'
#' Fit two-stage least squares models for qualitative outcomes to estimate the local probabilities of shift.
#'
#' @param Y Qualitative outcome before treatment. Must be labeled as \eqn{\{1, 2, \dots\}}.
#' @param D Binary treatment indicator.
#' @param Z Binary instrument.
#'
#' @return An object of class \code{causalQual}.
#'
#' @examples
#' \donttest{## Generate synthetic data.
#' set.seed(1986)
#'
#' data <- generate_qualitative_data_iv(1000, outcome_type = "ordered")
#'
#' Y <- data$Y
#' D <- data$D
#' Z <- data$Z
#'
#' ## Estimate local probabilities of shift.
#' fit <- causalQual_iv(Y, D, Z)
#'
#' summary(fit)
#' plot(fit)}
#'
#' @details
#' Under an instrumental-variables design, identification requires the instrument to be independent of potential outcomes and potential treatments (exogeneity), that the
#' instrument influences the outcome solely through its effect on treatment (exclusion restriction), that the instrument has a nonzero effect on treatment probability (relevance), and that the instrument can only
#' increase/decrease the treatment probability (monotonicity). If these assumptions hold, we can recover the local probabilities of shift for all classes:\cr
#'
#' \deqn{\delta_{m, L} := P(Y_i(1) = m | i \, is \, complier) - P(Y_i(0) = m | i \, is \, complier), \, m = 1, \dots, M.}
#'
#' \code{\link{causalQual_iv}} applies, for each class \code{m}, the standard two-stage least squares method to the binary variable \eqn{1(Y_i = m)}. Specifically, the routine first estimates
#' the following first-stage regression model via OLS:
#'
#' \deqn{D_i = \gamma_0 + \gamma_1 Z_i + \nu_i,}
#'
#' and constructs the predicted values \eqn{\hat{D}_i}. It then uses these predicted values in the second-stage regressions:
#'
#' \deqn{1(Y_i = m) = \alpha_{m0} + \alpha_{m1} \hat{D}_i + \epsilon_{mi}, \quad m = 1, \dots, M.}
#'
#' The OLS estimate \eqn{\hat{\alpha}_{m1}} of \eqn{\alpha_{m1}} is then our estimate of \eqn{\delta_{m, L}}. Standard errors are computed using conventional procedures and used to construct
#' conventional confidence intervals. All of this is done by calling the \code{\link[AER]{ivreg}} function.
#'
#' @importFrom AER ivreg
#'
#' @author Riccardo Di Francesco
#'
#' @seealso \code{\link{causalQual_soo}} \code{\link{causalQual_rd}} \code{\link{causalQual_did}}
#'
#' @export
causalQual_iv <- function(Y, D, Z) {
  ## 0.) Handle inputs and checks.
  if (any(!(Z %in% c(0, 1)))) stop("Invalid 'Z'. Instrument must be binary {0, 1}.", call. = FALSE)
  if (!all(sort(unique(Y)) == seq_len(max(Y)))) stop("Invalid 'Y'. Outcome must be coded with consecutive integers from 1 to max(Y).", call. = FALSE)

  classes <- sort(unique(Y))
  data <- data.frame("Y" = Y, "D" = D, "Z" = Z)

  ## 1.) Fit models. Then construct conventional confidence intervals.
  fits <- list()
  local_pshifts_hat <- numeric()
  local_pshifts_hat_se <- numeric()

  for (m in classes) {
    data$indicator <- as.numeric(Y == m)
    fit <- AER::ivreg(indicator ~ D | Z, data = data)

    fits[[paste0("Class_", m)]] <- fit
    local_pshifts_hat[m] <- coef(fit)["D"]
    local_pshifts_hat_se[m] <-  summary(fit)$coefficients["D", 2]
  }

  cis_lower <- local_pshifts_hat - 1.96 * local_pshifts_hat_se
  cis_upper <- local_pshifts_hat + 1.96 * local_pshifts_hat_se

  ## 2.) Output.
  output <- list()
  output$identification <- "iv"
  output$outcome_type <- NULL
  output$estimates <- local_pshifts_hat
  output$standard_errors <- local_pshifts_hat_se
  output$confidence_intervals <- list("lower" = cis_lower, "upper" = cis_upper)
  output$dr_scores <- NULL
  output$fits <- fits
  output$data <- list("Y" = Y, "Y_pre" = NULL, "Y_post" = NULL, "D" = D, "X" = NULL, "Z" = Z, "running_variable" = NULL, "cutoff" = NULL)

  class(output) <- "causalQual"
  return(output)
}


#' Causal Inference for Qualitative Outcomes under Difference-in-Differences
#'
#' Fit two-group/two-period models for qualitative outcomes to estimate the probabilities of shift on the treated.
#'
#' @param Y_pre Qualitative outcome before treatment. Must be labeled as \eqn{\{1, 2, \dots\}}.
#' @param Y_post Qualitative outcome after treatment. Must be labeled as \eqn{\{1, 2, \dots\}}.
#' @param D Binary treatment indicator.
#'
#' @return An object of class \code{causalQual}.
#'
#' @examples
#' \donttest{## Generate synthetic data.
#' set.seed(1986)
#'
#' data <- generate_qualitative_data_did(1000, assignment = "observational",
#'                                       outcome_type = "ordered")
#'
#' Y_pre <- data$Y_pre
#' Y_post <- data$Y_post
#' D <- data$D
#'
#' ## Estimate probabilities of shift on the treated.
#' fit <- causalQual_did(Y_pre, Y_post, D)
#'
#' summary(fit)
#' plot(fit)}
#'
#' @details
#' Under a difference-in-difference design, identification requires that the probabilities time shift for \eqn{Y_{is} (0)} for class \eqn{m} evolve similarly for the treated and control groups (parallel
#' trends on the probability mass functions of \eqn{Y_{is}(0)}). If this assumption holds, we can recover the probability of shift on the treated for class \eqn{m}:
#'
#' \deqn{\delta_{m, T} := P(Y_{it} (1) = m | D_i = 1) - P(Y_{it}(0) = m | D_i = 1).}
#'
#' \code{\link{causalQual_did}} applies, for each class \eqn{m}, the canonical two-group/two-period method to the binary variable \eqn{1(Y_{is} = m)}. Specifically,
#' consider the following linear model:
#'
#' \deqn{1(Y_{is} = m) = D_i \beta_{m1} + 1(s = t) \beta_{m2} + D_i 1(s = t) \beta_{m3} + \epsilon_{mis}.}
#'
#' The OLS estimate \eqn{\hat{\beta}_{m3}} of \eqn{\beta_{m3}} is our estimate of the probability shift on the treated for class \code{m}. Standard errors are clustered at the unit level and used to construct
#' conventional confidence intervals.
#'
#' @importFrom lmtest coeftest
#' @importFrom sandwich vcovCL
#' @importFrom stats lm
#' @importFrom stats as.formula pnorm
#'
#' @author Riccardo Di Francesco
#'
#' @seealso \code{\link{causalQual_soo}} \code{\link{causalQual_iv}} \code{\link{causalQual_rd}}
#'
#' @export
causalQual_did <- function(Y_pre, Y_post, D) {
  ## 0.) Handle inputs and checks.
  if (!all(sort(unique(Y_pre)) == seq_len(max(Y_pre)))) stop("Invalid 'Y_pre'. Outcome must be coded with consecutive integers from 1 to max(Y).", call. = FALSE)
  if (!all(sort(unique(Y_post)) == seq_len(max(Y_post)))) stop("Invalid 'Y_post'. Outcome must be coded with consecutive integers from 1 to max(Y).", call. = FALSE)
  if (is.null(Y_pre) | is.null(Y_post)) stop("When 'identification' is 'diff_in_diff', we need both 'Y_pre' and 'Y_post'.", call. = FALSE)

  classes_pre <- sort(unique(Y_pre))
  classes_post <- sort(unique(Y_post))

  if (any(sort(unique(classes_pre) != sort(unique(classes_post))))) stop("'Y_pre' and 'Y_post' have different classes.", call. = FALSE)
  classes <- classes_pre

  Y <- c(Y_pre, Y_post)
  D_new <- c(rep(D, 2))
  time <- c(rep(0, length(Y_pre)), rep(1, length(Y_pre)))
  unit_id <- rep(seq_len(length(Y_pre)), 2)

  data <- data.frame("Y" = Y, "D" = D_new, "time" = time, "unit_id" = unit_id)

  ## 1.) Fit models. Then construct conventional confidence intervals.
  fits <- list()
  pshifts_treated_hat <- numeric()
  pshifts_treated_hat_se <- numeric()

  for (m in classes) {
    data$indicator <- as.numeric(Y == m)
    fit <- stats::lm(indicator ~ D*time, data = data)
    cluster_se <- lmtest::coeftest(fit, vcov = sandwich::vcovCL(fit, cluster = ~ unit_id, type = "HC0"))

    fits[[paste0("Class_", m)]] <- fit
    pshifts_treated_hat[m] <- coef(fit)["D:time"]
    pshifts_treated_hat_se[m] <- cluster_se["D:time", 2]
  }

  cis_lower <- pshifts_treated_hat - 1.96 * pshifts_treated_hat_se
  cis_upper <- pshifts_treated_hat + 1.96 * pshifts_treated_hat_se

  ## 2.) Output.
  output <- list()
  output$identification <- "diff_in_diff"
  output$outcome_type <- NULL
  output$estimates <- pshifts_treated_hat
  output$standard_errors <- pshifts_treated_hat_se
  output$confidence_intervals <- list("lower" = cis_lower, "upper" = cis_upper)
  output$dr_scores <- NULL
  output$fits <- fits
  output$data <- list("Y" = NULL, "Y_pre" = Y_pre, "Y_post" = Y_post, "D" = D, "X" = NULL, "Z" = NULL, "running_variable" = NULL, "cutoff" = NULL)

  class(output) <- "causalQual"
  return(output)
}


#' Causal Inference for Qualitative Outcomes under Regression Discontinuity
#'
#' Fit local polynomial regression models for qualitative outcomes to estimate the probabilities of shift at the cutoff.
#'
#' @param Y Qualitative outcome. Must be labeled as \eqn{\{1, 2, \dots\}}.
#' @param running_variable Running variable determining treatment assignment.
#' @param cutoff Cutoff or threshold. Units with \code{running_variable < cutoff} are considered controls, while units with \code{running_variable >= cutoff} are considered treated.
#'
#' @return An object of class \code{causalQual}.
#'
#' @examples
#' \donttest{## Generate synthetic data.
#' set.seed(1986)
#'
#' data <- generate_qualitative_data_rd(1000, outcome_type = "ordered")
#'
#' Y <- data$Y
#' running_variable <- data$running_variable
#' cutoff <- data$cutoff
#'
#' ## Estimate probabilities of shift at the cutoff.
#' fit <- causalQual_rd(Y, running_variable, cutoff)
#'
#' summary(fit)
#' plot(fit)}
#'
#' @details
#' Under a regression discontinuity design, identification requires that the probability mass functions for class \eqn{m} of potential outcomes are continuous in the running variable (continuity). If this assumption holds,
#' we can recover the probability shift at the cutoff for class \eqn{m}:
#'
#' \deqn{\delta_{m, C} := P(Y_i (1) = m | Running_i = cutoff) - P(Y_i(0) = m | Running_i = cutoff).}
#'
#' \code{\link{causalQual_rd}} applies, for each class \eqn{m}, standard local polynomial estimators to the binary variable \eqn{1(Y_i = m)}. Specifically, the ruotine implements the
#' robust bias-corrected inference procedure of Calonico et al. (2014) (see the \code{\link[rdrobust]{rdrobust}} function).
#'
#' @import rdrobust
#'
#' @author Riccardo Di Francesco
#'
#' @seealso \code{\link{causalQual_soo}} \code{\link{causalQual_iv}} \code{\link{causalQual_did}}
#'
#' @export
causalQual_rd <- function(Y, running_variable, cutoff) {
  ## 0.) Handle inputs and checks.
  if (!all(sort(unique(Y)) == seq_len(max(Y)))) stop("Invalid 'Y'. Outcome must be coded with consecutive integers from 1 to max(Y).", call. = FALSE)
  if (is.null(running_variable)) stop("When 'identification' is 'rd', we need 'running_variable'.", call. = FALSE)
  if (is.null(cutoff)) stop("When 'identification' is 'rd', we need 'cutoff'.", call. = FALSE)

  classes <- sort(unique(Y))

  ## 1.) Fit models.
  fits <- list()
  pshifts_cutoff_hat <- numeric()
  pshifts_cutoff_hat_se <- numeric()
  pshifts_cutoff_hat_ci_lower <- numeric()
  pshifts_cutoff_hat_ci_upper <- numeric()

  for (m in classes) {
    indicator <- as.numeric(Y == m)
    fit <- rdrobust::rdrobust(indicator, running_variable, cutoff)

    fits[[paste0("Class_", m)]] <- fit
    pshifts_cutoff_hat[m] <- fit$Estimate[2]
    pshifts_cutoff_hat_se[m] <- fit$se[3]
    pshifts_cutoff_hat_ci_lower[m] <- fit$ci["Robust", 1]
    pshifts_cutoff_hat_ci_upper[m] <- fit$ci["Robust", 2]
  }

  ## 2.) Output.
  output <- list()
  output$identification <- "rd"
  output$outcome_type <- NULL
  output$estimates <- pshifts_cutoff_hat
  output$standard_errors <- pshifts_cutoff_hat_se
  output$confidence_intervals <- list("lower" = pshifts_cutoff_hat_ci_lower, "upper" = pshifts_cutoff_hat_ci_upper)
  output$dr_scores <- NULL
  output$fits <- fits
  output$data <- list("Y" = Y, "Y_pre" = NULL, "Y_post" = NULL, "D" = NULL, "X" = NULL, "Z" = NULL, "running_variable" = running_variable, "cutoff" = cutoff)

  class(output) <- "causalQual"
  return(output)
}
