#' Double Machine Learning for Qualitative Outcomes
#' 
#' Construct and average doubly robust scores for qualitative outcomes to estimate the probability shifts. Identification requires the treatment indicator to be (conditionally) independent of potential outcomes 
#' (unconfoundedness), and that each unit has a non-zero probability of being treated (common support).
#' 
#' @param Y Ordered non-numeric outcome. Must be labeled as \eqn{\{1, 2, \dots\}}.
#' @param D Binary treatment indicator.
#' @param X Covariate matrix (no intercept).
#' @param outcome_type String controlling the outcome type. Must be either \code{"multinomial"} or \code{"ordered"}. Affects estimation of conditional class probabilities.
#' @param K Number of folds for nuisance functions estimation.
#' 
#' @return An object of class \code{causalQual}.
#'   
#' @details
#' For each class \eqn{m}, the doubly robust score is defined as:
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
#' @import ocf grf
#' @importFrom caret createFolds
#' @importFrom magrittr %>%
#' 
#' @seealso \code{\link{causalQual}}
#' 
#' @author Riccardo Di Francesco
#' 
#' @keywords internal
causalQual_soo <- function(Y, D, X, outcome_type, K = 5) {
  ## 0.) Handling inputs and checks.
  if (!all(sort(unique(Y)) == seq_len(max(Y)))) stop("Invalid 'Y'. Outcome must be coded with consecutive integers from 1 to max(Y).", call. = FALSE)
  if (any(!(D %in% c(0, 1)))) stop("Invalid 'D'. Treatment must be binary {0, 1}.", call. = FALSE)
  if (K < 1 | K %% 1 != 0) stop("Invalid 'K'. Number of folds must be a positive integer.", call. = FALSE)
  if (!(outcome_type %in% c("multinomial", "ordered"))) stop("Invalid 'outcome_type'. Must be either 'multinomial' or 'ordered'.", call. = FALSE)
  
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
  
  ## 2.) Estimate conditional class probabilities and propensity score via cross-fitting.
  folds <- caret::createFolds(seq_len(length(Y)), K, list = TRUE)
  
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
    } else if (outcom_type == "ordered") {
      forests1 <- ocf::ocf(Y[treated_idx_train], X[treated_idx_train, ], honesty = TRUE, alpha = 0.2)
      forests0 <- ocf::ocf(Y[control_idx_train], X[control_idx_train, ], honesty = TRUE, alpha = 0.2)
    }
    
    forest_pscore <- grf::regression_forest(as.matrix(X[train_idx, ]), D[train_idx])
    
    ## Predict on holdout.
    pclasses1[holdout_idx, ] <- stats::predict(forests1, X[holdout_idx, ])$probabilities
    pclasses0[holdout_idx, ] <- stats::predict(forests0, X[holdout_idx, ])$probabilities
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
  
  ## 4.) Pick mean and standard deviation of doubly-robust scores. Then construct confidence intervals.
  pshifts_hat <- sapply(gammas, mean)
  pshifts_hat_se <- sapply(seq_along(gammas), function(m) mean((gammas[[m]] - pshifts_hat[m])^2) / length(Y)) %>%
    sqrt()
  
  ci_lower <- pshifts_hat - 1.96 * pshifts_hat_se
  ci_upper <- pshifts_hat + 1.96 * pshifts_hat_se
  
  ## 5.) Output.
  return(list("estimates" = pshifts_hat, "standard_errors" = pshifts_hat_se, "ci_lower" = ci_lower, "ci_upper" = ci_upper))
}
