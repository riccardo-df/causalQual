#' Generate Qualitative Data (Selection-on-Observables)
#'
#' Generate a synthetic data set with qualitative outcomes under a selection-on-observables design. The data include a binary treatment indicator and a matrix of covariates. The treatment is either
#' independent or conditionally (on the covariates) independent of potential outcomes, depending on users' choices.
#'
#' @param n Sample size.
#' @param assignment String controlling treatment assignment. Must be either \code{"randomized"} (random assignment) or \code{"observational"} (random assigment conditional on the generated covariates).
#' @param outcome_type String controlling the outcome type. Must be either \code{"multinomial"} or \code{"ordered"}. Affects how potential outcomes are generated.
#'
#' @return A list storing a data frame with the observed data, the true propensity score, and the true probabilities of shift.
#'
#' @examples
#' \donttest{## Generate synthetic data.
#' set.seed(1986)
#'
#' data <- generate_qualitative_data_soo(100,
#'                                       assignment = "observational",
#'                                       outcome_type = "ordered")
#'
#' data$pshifts}
#'
#' @details
#' ## Outcome type
#' Potential outcomes are generated differently according to \code{outcome_type}. If \code{outcome_type == "multinomial"}, \code{\link{generate_qualitative_data_soo}} computes linear predictors for each class using the covariates:
#'
#' \deqn{\eta_{mi} (d) = \beta_{m1}^d X_{i1} + \beta_{m2}^d X_{i2} + \beta_{m3}^d X_{i3}, \quad d = 0, 1,}
#'
#' and then transforms \eqn{\eta_{mi} (d)} into valid probability distributions using the softmax function:
#'
#' \deqn{P(Y_i(d) = m | X_i) = \frac{\exp(\eta_{mi} (d))}{\sum_{m'} \exp(\eta_{m'i}(d))}, \quad d = 0, 1.}
#'
#' It then generates potential outcomes \eqn{Y_i(1)} and \eqn{Y_i(0)} by sampling from \{1, 2, 3\} using \eqn{P(Y_i(d) = m | X_i), \, d = 0, 1}.\cr
#'
#' If instead \code{outcome_type == "ordered"}, \code{\link{generate_qualitative_data_soo}} first generates latent potential outcomes:
#'
#' \deqn{Y_i^* (d) = \tau d + X_{i1} + X_{i2} + X_{i3} + N (0, 1), \quad d = 0, 1,}
#'
#' with \eqn{\tau = 2}. It then constructs \eqn{Y_i (d)} by discretizing \eqn{Y_i^* (d)} using threshold parameters \eqn{\zeta_1 = 2} and \eqn{\zeta_2 = 4}. Then,
#'
#' \deqn{P(Y_i(d) = m | X_i) = P(\zeta_{m-1} < Y_i^*(d) \leq \zeta_m | X_i) = \Phi (\zeta_m - \sum_j X_{ij} - \tau d) - \Phi (\zeta_{m-1} - \sum_j X_{ij} - \tau d), \quad d = 0, 1,}
#'
#' which allows us to analytically compute the probabilities of shift.
#'
#' ## Treatment assignment
#' Treatment is always assigned as \eqn{D_i \sim \text{Bernoulli}(\pi(X_i))}. If \code{assignment == "randomized"}, then the propensity score is specified as \eqn{\pi(X_i) = P ( D_i = 1 | X_i)) = 0.5}.
#' If instead \code{assignment == "observational"}, then \eqn{\pi(X_i) = (X_{i1} + X_{i3}) / 2}.
#'
#' ## Other details
#' The function always generates three independent covariates from \eqn{U(0,1)}. Observed outcomes \eqn{Y_i} are always constructed using the usual observational rule.\cr
#'
#' Controlling for \eqn{X_{i1}} and \eqn{X_{i3}} is sufficient for selection-on-observables to hold.
#'
#' @importFrom stats rbinom
#' @importFrom stats runif
#'
#' @author Riccardo Di Francesco
#'
#' @seealso \code{\link{generate_qualitative_data_iv}} \code{\link{generate_qualitative_data_rd}} \code{\link{generate_qualitative_data_did}}
#'
#' @export
generate_qualitative_data_soo <- function(n, assignment, outcome_type) {
  ## 0.) Handling inputs and checks.
  if (n < 1 | n %% 1 != 0) stop("Invalid 'n'. This must be an integer greater than or equal to one.", call. = FALSE)
  if (!(assignment %in% c("randomized", "observational"))) stop("Invalid 'assignment'. This must be either 'randomized' or 'observational'.", call. = FALSE)
  if (!(outcome_type %in% c("multinomial", "ordered"))) stop("Invalid 'outcome_type'. Must be either 'multinomial' or 'ordered'.", call. = FALSE)

  k <- 3
  classes <- seq_len(3)

  ## 1.) Generate covariates.
  X <- matrix(stats::runif(n * k), ncol = k)
  colnames(X) <- paste0("X", seq(k))

  ## 2.) Construct potential outcomes.
  if (outcome_type == "multinomial") {
    ## Construct linear predictors.
    beta1 <- matrix(c(0.5,  0.3, -0.2, # Class 1.
                      -0.2,  0.4,  0.1, # Class 2.
                      0.1, -0.3, 0.5), # Class 3.
                    ncol = 3, byrow = TRUE)
    beta0 <- matrix(c(0.7,  0.7, -0.2, # Class 1.
                      -0.2,  0.4,  0.1, # Class 2.
                      -0.2, -0.5,  0.1), # Class 3.
                    ncol = 3, byrow = TRUE)

    logits1 <- X %*% beta1
    logits0 <- X %*% beta0

    ## Softmax function to get valid probabilities.
    pmass1 <- softmax(logits1)
    pmass0 <- softmax(logits0)

    ## Sample potential outcomes.
    Y1 <- apply(pmass1, 1, function(p) sample(classes, size = 1, prob = p))
    Y0 <- apply(pmass0, 1, function(p) sample(classes, size = 1, prob = p))
  } else if (outcome_type == "ordered") {
    ## Generate latent potential outcomes with ATE tau.
    tau <- 2
    Ystar1 <- tau + rowSums(X) + rnorm(n)
    Ystar0 <- rowSums(X) + rnorm(n)

    ## Discretize to construct potential outcomes.
    zeta <- c(2, 3)

    Y1 <- cut(Ystar1, breaks = c(-Inf, zeta, Inf), labels = classes, right = TRUE)
    Y0 <- cut(Ystar0, breaks = c(-Inf, zeta, Inf), labels = classes, right = TRUE)

    Y1 <- as.integer(as.character(Y1))
    Y0 <- as.integer(as.character(Y0))

    ## Probability mass functions (needed to compute probabilities of shift later).
    pmass1 <- cbind(pnorm(zeta[1] - rowSums(X) - tau),
                    pnorm(zeta[2] - rowSums(X) - tau) - pnorm(zeta[1] - rowSums(X) - tau),
                    1 - pnorm(zeta[2] - rowSums(X) - tau))

    pmass0 <- cbind(pnorm(zeta[1] - rowSums(X)),
                    pnorm(zeta[2] - rowSums(X)) - pnorm(zeta[1] - rowSums(X)),
                    1 - pnorm(zeta[2] - rowSums(X)))
  }

  ## 3.) Define propensity scores and assign treatment.
  if (assignment == "randomized") {
    pscore <- rep(0.5, n)
  } else if (assignment == "observational") {
    pscore <- (X[, 1] + X[, 3]) / 2
  }

  D <- stats::rbinom(n, 1, pscore)

  ## 4.) Generate observed outcomes.
  Y <- D * Y1 + (1 - D) * Y0

  ## 5.) Probabilities of shift.
  pshifts <- colMeans(pmass1 - pmass0)

  ## 6.) Output.
  return(list("Y" = Y, "D" = D, "X" = X, "pscore" = pscore, "pshifts" = pshifts))
}


#' Generate Qualitative Data (Instrumental Variables)
#'
#' Generate a synthetic data set with qualitative outcomes under an instrumental variables design. The data include a binary treatment indicator and a binary instrument. Potential outcomes
#' and potential treatments are independent of the instrument. Moreover, the instrument does not directly impact potential outcomes, has an impact on treatment probability, and can only increase the probability
#' of treatment.
#'
#' @param n Sample size.
#' @param outcome_type String controlling the outcome type. Must be either \code{"multinomial"} or \code{"ordered"}. Affects how potential outcomes are generated.
#'
#' @return A list storing a data frame with the observed data, the true propensity score, the true instrument propensity score, and the true local probabilities of shift.
#'
#' @examples
#' \donttest{## Generate synthetic data.
#' set.seed(1986)
#'
#' data <- generate_qualitative_data_iv(100,
#'                                      outcome_type = "ordered")
#'
#' data$local_pshifts}
#'
#' @details
#' ## Outcome type
#' Potential outcomes are generated differently according to \code{outcome_type}. If \code{outcome_type == "multinomial"}, \code{\link{generate_qualitative_data_iv}} computes linear predictors for each class using the covariates:
#'
#' \deqn{\eta_{mi} (d) = \beta_{m1}^d X_{i1} + \beta_{m2}^d X_{i2} + \beta_{m3}^d X_{i3}, \quad d = 0, 1,}
#'
#' and then transforms \eqn{\eta_{mi} (d)} into valid probability distributions using the softmax function:
#'
#' \deqn{P(Y_i(d) = m | X_i) = \frac{\exp(\eta_{mi} (d))}{\sum_{m'} \exp(\eta_{m'i}(d))}, \quad d = 0, 1.}
#'
#' It then generates potential outcomes \eqn{Y_i(1)} and \eqn{Y_i(0)} by sampling from \{1, 2, 3\} using \eqn{P_i(Y(d) = m | X), \, d = 0, 1}.\cr
#'
#' If instead \code{outcome_type == "ordered"}, \code{\link{generate_qualitative_data_iv}} first generates latent potential outcomes:
#'
#' \deqn{Y_i^* (d) = \tau d + X_{i1} + X_{i2} + X_{i3} + N (0, 1), \quad d = 0, 1,}
#'
#' with \eqn{\tau = 2}. It then constructs \eqn{Y_i (d)} by discretizing \eqn{Y_i^* (d)} using threshold parameters \eqn{\zeta_1 = 2} and \eqn{\zeta_2 = 4}. Then,
#'
#' \deqn{P(Y_i(d) = m | X_i) = P(\zeta_{m-1} < Y_i^*(d) \leq \zeta_m | X_i) = \Phi (\zeta_m - \sum_j X_{ij} - \tau d) - \Phi (\zeta_{m-1} - \sum_j X_{ij} - \tau d), \quad d = 0, 1,}
#'
#' which allows us to analytically compute the local probabilities of shift.
#'
#' ## Treatment assignment and instrument
#' The instrument is always generated as \eqn{Z_i \sim \text{Bernoulli}(0.5)}. Treatment is always modeled as \eqn{D_i \sim \text{Bernoulli}(\pi(X_i, Z_i))}, with
#' \eqn{\pi(X_i, Z_i) = P ( D_i = 1 | X_i, Z_i)) = (X_{i1} + X_{i3} + Z_i) / 3}. Thus, \eqn{Z_i} can increase the probability of treatment intake but cannot decrease it.
#'
#' ## Other details
#' The function always generates three independent covariates from \eqn{U(0,1)}. Observed outcomes \eqn{Y_i} are always constructed using the usual observational rule.\cr
#'
#' The DGPs outlined above ensure that \eqn{Z_i} is a valid instrument for \code{D_i}.
#'
#' @importFrom stats rbinom
#' @importFrom stats runif
#'
#' @author Riccardo Di Francesco
#'
#' @seealso \code{\link{generate_qualitative_data_soo}} \code{\link{generate_qualitative_data_rd}} \code{\link{generate_qualitative_data_did}}
#'
#' @export
generate_qualitative_data_iv <- function(n, outcome_type) {
  ## 0.) Handling inputs and checks.
  if (n < 1 | n %% 1 != 0) stop("Invalid 'n'. This must be an integer greater than or equal to one.", call. = FALSE)
  if (!(outcome_type %in% c("multinomial", "ordered"))) stop("Invalid 'outcome_type'. Must be either 'multinomial' or 'ordered'.", call. = FALSE)

  k <- 3
  classes <- seq_len(3)

  ## 1.) Generate covariates.
  X <- matrix(runif(n * k), ncol = k)
  colnames(X) <- paste0("X", seq(k))

  ## 2.) Construct potential outcomes-
  if (outcome_type == "multinomial") {
    ## Construct linear predictors.
    beta1 <- matrix(c(0.5,  0.3, -0.2, # Class 1.
                      -0.2,  0.4,  0.1, # Class 2.
                      0.1, -0.3, 0.5), # Class 3.
                    ncol = 3, byrow = TRUE)
    beta0 <- matrix(c(0.7,  0.7, -0.2, # Class 1.
                      -0.2,  0.4,  0.1, # Class 2.
                      -0.2, -0.5,  0.1), # Class 3.
                    ncol = 3, byrow = TRUE)

    logits1 <- X %*% beta1
    logits0 <- X %*% beta0

    ## Softmax function to get valid probabilities.
    pmass1 <- exp(logits1) / rowSums(exp(logits1))
    pmass0 <- exp(logits0) / rowSums(exp(logits0))

    ## Sample potential outcomes.
    Y1 <- apply(pmass1, 1, function(p) sample(classes, size = 1, prob = p))
    Y0 <- apply(pmass0, 1, function(p) sample(classes, size = 1, prob = p))
  } else if (outcome_type == "ordered") {
    ## Generate latent potential outcomes with ATE tau.
    tau <- 2

    Ystar1 <- tau + rowSums(X) + rnorm(n)
    Ystar0 <- rowSums(X) + rnorm(n)

    ## Discretize to construct potential outcomes.
    zeta <- c(2, 3)

    Y1 <- cut(Ystar1, breaks = c(-Inf, zeta, Inf), labels = classes, right = TRUE)
    Y0 <- cut(Ystar0, breaks = c(-Inf, zeta, Inf), labels = classes, right = TRUE)

    Y1 <- as.integer(as.character(Y1))
    Y0 <- as.integer(as.character(Y0))

    ## Probability mass functions (needed to compute probabilities of shift later).
    pmass1 <- cbind(pnorm(zeta[1] - rowSums(X) - tau),
                    pnorm(zeta[2] - rowSums(X) - tau) - pnorm(zeta[1] - rowSums(X) - tau),
                    1 - pnorm(zeta[2] - rowSums(X) - tau))

    pmass0 <- cbind(pnorm(zeta[1] - rowSums(X)),
                    pnorm(zeta[2] - rowSums(X)) - pnorm(zeta[1] - rowSums(X)),
                    1 - pnorm(zeta[2] - rowSums(X)))
  }

  ## 3.) Define propensity scores and assign instrument and treatment.
  instrument_pscore <- rep(0.5, n)
  Z <- stats::rbinom(n, 1, instrument_pscore)

  pscore <- (X[, 1] + X[, 3] + Z) / 3
  D <- stats::rbinom(n, 1, pscore)

  ## 4.) Generate observed outcomes.
  Y <- D * Y1 + (1 - D) * Y0

  ## 5.) Probabilities of shift.
  local_pshifts <- colMeans(pmass1[D == Z, ] - pmass0[D == Z, ])

  ## 6.) Output.
  return(list("Y" = Y, "D" = D, "Z" = Z, "pscore" = pscore, "instrument_pscore" = instrument_pscore, "local_pshifts" = local_pshifts))
}


#' Generate Qualitative Data (Regression Discontinuity)
#'
#' Generate a synthetic data set with qualitative outcomes under a regression discontinuity design. The data include a binary treatment indicator and a single covariate (the running variable). The conditional
#' probability mass fuctions of potential outcomes are continuous in the running variable.
#'
#' @param n Sample size.
#' @param outcome_type String controlling the outcome type. Must be either \code{"multinomial"} or \code{"ordered"}. Affects how potential outcomes are generated.
#'
#' @return A list storing a data frame with the observed data, and the true probabilities of shift at the cutoff.
#'
#' @examples
#' \donttest{## Generate synthetic data.
#' set.seed(1986)
#'
#' data <- generate_qualitative_data_rd(100,
#'                                      outcome_type = "ordered")
#'
#' data$pshifts_cutoff}
#'
#' @details
#' ## Outcome type
#' Potential outcomes are generated differently according to \code{outcome_type}. If \code{outcome_type == "multinomial"}, \code{\link{generate_qualitative_data_rd}} computes linear predictors for each class using the covariates:
#'
#' \deqn{\eta_{mi} (d) = \beta_{m1}^d X_{i1} + \beta_{m2}^d X_{i2} + \beta_{m3}^d X_{i3}, \quad d = 0, 1,}
#'
#' and then transforms \eqn{\eta_{mi} (d)} into valid probability distributions using the softmax function:
#'
#' \deqn{P(Y_i(d) = m | X_i) = \frac{\exp(\eta_{mi} (d))}{\sum_{m'} \exp(\eta_{m'i}(d))}.}
#'
#' It then generates potential outcomes \eqn{Y_i(1)} and \eqn{Y_i(0)} by sampling from \{1, 2, 3\} using \eqn{P(Y_i(d) = m | X_i), \, d = 0, 1}.\cr
#'
#' If instead \code{outcome_type == "ordered"}, \code{\link{generate_qualitative_data_rd}} first generates latent potential outcomes:
#'
#' \deqn{Y_i^* (d) = \tau d + X_{i1} + X_{i2} + X_{i3} + N (0, 1), \quad d = 0, 1,}
#'
#' with \eqn{\tau = 2}. It then constructs \eqn{Y_i (d)} by discretizing \eqn{Y_i^* (d)} using threshold parameters \eqn{\zeta_1 = 2} and \eqn{\zeta_2 = 4}. Then,
#'
#' \deqn{P(Y_i(d) = m) = P(\zeta_{m-1} < Y_i^*(d) \leq \zeta_m) = \Phi (\zeta_m - \sum_j X_{ij} - \tau d) - \Phi (\zeta_{m-1} - \sum_j X_{ij} - \tau d), \quad d = 0, 1,}
#'
#' which allows us to analytically compute the probabilities of shift at the cutoff.
#'
#' ## Treatment assignment
#' Treatment is always assigned as \eqn{D_i = 1(X_i \geq 0.5)}.
#'
#' ## Other details
#' The function always generates three independent covariates from \eqn{U(0,1)}. Observed outcomes \eqn{Y_i} are always constructed using the usual observational rule.\cr
#'
#' The DGPs outlined above ensure identification under a standard regression discontinuity design.
#'
#' @importFrom stats runif
#'
#' @author Riccardo Di Francesco
#'
#' @seealso \code{\link{generate_qualitative_data_soo}} \code{\link{generate_qualitative_data_iv}} \code{\link{generate_qualitative_data_did}}
#'
#' @export
generate_qualitative_data_rd <- function(n, outcome_type) {
  ## 0.) Handling inputs and checks.
  if (n < 1 | n %% 1 != 0) stop("Invalid 'n'. This must be an integer greater than or equal to one.", call. = FALSE)
  if (!(outcome_type %in% c("multinomial", "ordered"))) stop("Invalid 'outcome_type'. Must be either 'multinomial' or 'ordered'.", call. = FALSE)

  k <- 3
  classes <- seq_len(3)
  cutoff <- 0.5

  ## 1.) Generate covariates.
  X <- matrix(stats::runif(n * k), ncol = k)
  colnames(X) <- paste0("X", seq(k))

  ## 2.) Construct potential outcomes.
  if (outcome_type == "multinomial") {
    ## Construct linear predictors.
    beta1 <- matrix(c(0.5,  0.3, -0.2, # Class 1.
                      -0.2,  0.4,  0.1, # Class 2.
                      0.1, -0.3, 0.5), # Class 3.
                    ncol = 3, byrow = TRUE)
    beta0 <- matrix(c(0.7,  0.7, -0.2, # Class 1.
                      -0.2,  0.4,  0.1, # Class 2.
                      -0.2, -0.5,  0.1), # Class 3.
                    ncol = 3, byrow = TRUE)

    logits1 <- X %*% beta1
    logits0 <- X %*% beta0

    ## Softmax function to get valid probabilities.
    pmass1 <- softmax(logits1)
    pmass0 <- softmax(logits0)

    ## Sample potential outcomes.
    Y1 <- apply(pmass1, 1, function(p) sample(classes, size = 1, prob = p))
    Y0 <- apply(pmass0, 1, function(p) sample(classes, size = 1, prob = p))
  } else if (outcome_type == "ordered") {
    ## Generate latent potential outcomes with ATE tau.
    tau <- 2
    Ystar1 <- tau + rowSums(X) + rnorm(n)
    Ystar0 <- rowSums(X) + rnorm(n)

    ## Discretize to construct potential outcomes.
    zeta <- c(2, 3)

    Y1 <- cut(Ystar1, breaks = c(-Inf, zeta, Inf), labels = classes, right = TRUE)
    Y0 <- cut(Ystar0, breaks = c(-Inf, zeta, Inf), labels = classes, right = TRUE)

    Y1 <- as.integer(as.character(Y1))
    Y0 <- as.integer(as.character(Y0))

    ## Probability mass functions (needed to compute probabilities of shift later).
    pmass1 <- cbind(pnorm(zeta[1] - rowSums(X) - tau),
                    pnorm(zeta[2] - rowSums(X) - tau) - pnorm(zeta[1] - rowSums(X) - tau),
                    1 - pnorm(zeta[2] - rowSums(X) - tau))

    pmass0 <- cbind(pnorm(zeta[1] - rowSums(X)),
                    pnorm(zeta[2] - rowSums(X)) - pnorm(zeta[1] - rowSums(X)),
                    1 - pnorm(zeta[2] - rowSums(X)))
  }

  ## 3.) Define propensity scores and assign treatment.
  D <- as.numeric(X[, 1] >= cutoff)

  ## 4.) Generate observed outcomes.
  Y <- D * Y1 + (1 - D) * Y0

  ## 5.) Probabilities of shift at cutoff. Because of continuous running variable, we focus on a really small neighborhood.
  tolerance <- 0.01
  pshifts_cutoff <- colMeans(pmass1[abs(X[, 1] - cutoff) < tolerance, ] - pmass0[abs(X[, 1] - cutoff) < tolerance, ])

  ## 6.) Output.
  return(list("Y" = Y, "running_variable" = X[, 1], "cutoff" = cutoff, "pshifts_cutoff" = pshifts_cutoff))
}


#' Generate Qualitative Data (Difference-in-Differences)
#'
#' Generate a synthetic data set with qualitative outcomes under a difference-in-differences design. The data include two time periods, a binary treatment indicator (applied only in the second period),
#' and a matrix of covariates. Probabilities time shift among the treated and control groups evolve similarly across the two time periods (parallel trends on the probability mass functions).
#'
#' @param n Sample size.
#' @param time_periods Controls how many time periods are generated.
#' @param assignment String controlling treatment assignment. Must be either \code{"randomized"} (random assignment) or \code{"observational"} (assignment based on covariates).
#' @param outcome_type String controlling the outcome type. Must be either \code{"multinomial"} or \code{"ordered"}.
#' @param assignment_type String controlling the type of assignment timing. Must be either \code{"block"} (all treated units get treatment at the same time period) or \code{"staggered"} (units receive treatment at different time periods).
#' @return A list storing a data frame with the observed data, the true propensity score, and the true probabilities of shift on the treated.
#'
#' @details
#' ## Outcome type
#' Potential outcomes are generated differently according to \code{outcome_type}. If \code{outcome_type == "multinomial"}, \code{\link{generate_qualitative_data_did}} computes linear predictors for each class using the covariates:
#'
#' \deqn{\eta_{mi} (d, s) = \beta_{m1}^d X_{i1} + \beta_{m2}^d X_{i2} + \beta_{m3}^d X_{i3}, \quad d = 0, 1, \quad s = t-1, t,}
#'
#' and then transforms \eqn{\eta_{mi} (d, s)} into valid probability distributions using the softmax function:
#'
#' \deqn{P(Y_{is}(d) = m | X_i) = \frac{\exp(\eta_{mi} (d, s))}{\sum_{m'} \exp(\eta_{m'i}(d, s))}, \quad d = 0, 1, \quad s = t-1, t.}
#'
#' It then generates potential outcomes \eqn{Y_{it-1}(1)}, \eqn{Y_{it}(1)}, \eqn{Y_{it-1}(0)}, and \eqn{Y_{it}(0)} by sampling from \{1, 2, 3\} using \eqn{P(Y(d, s) = m \mid X), \, d = 0, 1, \, s = t-1, t}.\cr
#'
#' If instead \code{outcome_type == "ordered"}, \code{\link{generate_qualitative_data_did}} first generates latent potential outcomes:
#'
#' \deqn{Y_i^* (d, s) = \tau d + X_{i1} + X_{i2} + X_{i3} + N (0, 1), \quad d = 0, 1, \quad s = t-1, t,}
#'
#' with \eqn{\tau = 2}. It then constructs \eqn{Y_i (d, s)} by discretizing \eqn{Y_i^* (d, s)} using threshold parameters \eqn{\zeta_1 = 2} and \eqn{\zeta_2 = 4}. Then,
#'
#' \deqn{P(Y_i(d, s) = m | X_i) = P(\zeta_{m-1} < Y_i^*(d, s) \leq \zeta_m | X_i) = \Phi (\zeta_m - \sum_j X_{ij} - \tau d) - \Phi (\zeta_{m-1} - \sum_j X_{ij} - \tau d), \quad d = 0, 1, \quad s = t-1, t,}
#'
#' which allows us to analytically compute the probabilities of shift on the treated for each group (determined by timing of treatment) at each time after treatment occurs for that group.\cr
#'
#' Notice that, by design, the probabilities of shift on the treated of each group are constant over time. However, effects can vary across groups due to variation in covariate distributions.
#'
#' ## Treatment assignment
#' Treatment is always assigned as \eqn{D_i \sim \text{Bernoulli}(\pi(X_i))}. If \code{assignment == "randomized"}, then the propensity score is specified as \eqn{\pi(X_i) = P ( D_i = 1 | X_i)) = 0.5}.
#' If instead \code{assignment == "observational"}, then \eqn{\pi(X_i) = (X_{i1} + X_{i3}) / 2}.\cr
#'
#' If \code{assignment_type == "block"}, then all treated units receive treatment at time \code{ceiling(time_periods / 2)}. If instead \code{assignment_type == "staggered"}, then the routine samples
#' the treatment date among \code{2:(time_periods-1)} for treated units (essentially, each treated unit is equally likely to be treated at one date between the second and the second-to-last time period).
#'
#' ## Other details
#' The function always generates three independent covariates from \eqn{U(0,1)}. Observed outcomes \eqn{Y_{is}} are always constructed using the usual observational rule.\cr
#'
#' The DGPs outlined above ensure that a parallel trends assumption holds for all classes.
#'
#' @examples
#' \donttest{## Generate synthetic data.
#' set.seed(1986)
#'
#' data <- generate_qualitative_data_did(100,
#'                                       time_periods = 5,
#'                                       assignment = "observational",
#'                                       assignment_type = "staggered",
#'                                       outcome_type = "ordered")
#'
#' data$pshifts_treated}
#'
#' @importFrom stats rbinom runif rnorm
#'
#' @author Riccardo Di Francesco
#'
#' @seealso \code{\link{generate_qualitative_data_soo}} \code{\link{generate_qualitative_data_iv}} \code{\link{generate_qualitative_data_rd}}
#'
#' @export
generate_qualitative_data_did <- function(n, time_periods = 2, assignment, assignment_type = "block", outcome_type) {
  ## 0.) Input validation
  if (n < 1 | n %% 1 != 0) stop("Invalid 'n'. Must be an integer greater than or equal to one.", call. = FALSE)
  if (time_periods < 2) stop("Invalid 'time_periods'. Must be at least 2.", call. = FALSE)
  if (time_periods %% 1 != 0) stop("Invalid 'time_periods'. Must be integer.", call. = FALSE)
  if (!(assignment %in% c("randomized", "observational"))) stop("Invalid 'assignment'. Must be either 'randomized' or 'observational'.", call. = FALSE)
  if (!(assignment_type %in% c("block", "staggered"))) stop("Invalid 'assignment_type'. Must be either 'block' or 'staggered'.", call. = FALSE)
  if (!(outcome_type %in% c("multinomial", "ordered"))) stop("Invalid 'outcome_type'. Must be either 'multinomial' or 'ordered'.", call. = FALSE)

  k <- 3
  classes <- seq_len(3)

  ## 1.) Generate covariates.
  X <- matrix(stats::runif(n * k), ncol = k) # Notice they are time-invariant.
  colnames(X) <- paste0("X", seq(k))

  ## 2.) Construct potential outcomes. Here we need to allocate matrices for potential outcomes and loop over time periods. At each loop iteration, we sample potential outcomes for that time period. Models are time-invarying.
  Y1 <- matrix(NA, nrow = n, ncol = time_periods) # Each row is a unit, each column is a time period.
  Y0 <- matrix(NA, nrow = n, ncol = time_periods)

  for (t in seq_len(time_periods)) {
    if (outcome_type == "multinomial") {
      ## Compute linear predictors for probabilities.
      beta1 <- matrix(c(0.5,  0.3, -0.2, # Class 1.
                        -0.2,  0.4,  0.1, # Class 2.
                        0.1, -0.3, 0.5), # Class 3.
                      ncol = 3, byrow = TRUE)
      beta0 <- matrix(c(0.7,  0.7, -0.2, # Class 1.
                        -0.2,  0.4,  0.1, # Class 2.
                        -0.2, -0.5,  0.1), # Class 3.
                      ncol = 3, byrow = TRUE)

      logits1 <- X %*% beta1
      logits0 <- X %*% beta0

      ## Softmax function to get valid probabilities. Notice they are time-invarying; since it's fast computing them, it's easier to do this within the loop.
      pmass1 <- softmax(logits1)
      pmass0 <- softmax(logits0)

      ## Sample potential outcomes.
      Y1[, t] <- apply(pmass1, 1, function(p) sample(classes, 1, prob = p))
      Y0[, t] <- apply(pmass0, 1, function(p) sample(classes, 1, prob = p))
    } else if (outcome_type == "ordered") {
      ## Generate latent potential outcomes with ATT tau.
      tau <- 2

      Ystar1 <- tau + rowSums(X) + rnorm(n)
      Ystar0 <- rowSums(X) + rnorm(n)

      ## Discretize to construct potential outcomes.
      zeta <- c(2, 3)

      Y1t <- cut(Ystar1, breaks = c(-Inf, zeta, Inf), labels = classes, right = TRUE)
      Y0t <- cut(Ystar0, breaks = c(-Inf, zeta, Inf), labels = classes, right = TRUE)

      Y1[, t] <- as.integer(as.character(Y1t))
      Y0[, t] <- as.integer(as.character(Y0t))

      ## Probability mass functions (needed to compute probabilities of shift on the treated later on). Notice they are time-invarying; since it's fast computing them, it's easier to do this within the loop.
      pmass1 <- cbind(pnorm(zeta[1] - rowSums(X) - tau),
                      pnorm(zeta[2] - rowSums(X) - tau) - pnorm(zeta[1] - rowSums(X) - tau),
                      1 - pnorm(zeta[2] - rowSums(X) - tau))

      pmass0 <- cbind(pnorm(zeta[1] - rowSums(X)),
                      pnorm(zeta[2] - rowSums(X)) - pnorm(zeta[1] - rowSums(X)),
                      1 - pnorm(zeta[2] - rowSums(X)))
    }
  }

  ## 3.) Define propensity scores and assign treatment.
  # Propensity score based on whether assignment is random or observational.
  if (assignment == "randomized") {
    pscore <- rep(0.5, n)
  } else if (assignment == "observational") {
    pscore <- (X[, 1] + X[, 3]) / 2
  }

  # Set treated units.
  D <- stats::rbinom(n, 1, pscore)

  # Set when treatment occurs. Initialize all units as never treated by default (treatment date set beyond last time period for simplicity). Then modify treatment start date for treated units.
  treatment_start <- rep(time_periods + 1, n) # For each units, determines when treatment occurs.

  if (assignment_type == "block") {
    block_time <- ceiling(time_periods / 2)
    treatment_start[D == 1] <- block_time
  } else if (assignment_type == "staggered") {
    treatment_start[D == 1] <- sample(2:(time_periods-1), sum(D), replace = TRUE)
  }

  ## 4.) Generate observed outcomes. We need again to loop over time periods.
  Y <- matrix(NA, nrow = n, ncol = time_periods)

  for (t in seq_len(time_periods)) {
    is_treated <- treatment_start <= t # Select units treated at time t or before. We observe Y1 for them.
    Y[is_treated, t] <- Y1[is_treated, t]
    Y[!is_treated, t] <- Y0[!is_treated, t]
  }

  ## 5.) Probabilities of shift on the treated. For each group (defined by time of receiving treatment), compute average effects at each time t > g and save in list. If assignment is "block", we have only one group.
  # Pre-allocate memory.
  treated_groups <- sort(unique(treatment_start[treatment_start <= time_periods]))
  pshifts_treated_list <- list()

  # Loop over groups.
  for (g in treated_groups) {
    # Identify units in this group.
    units_in_group <- which(treatment_start == g)

    # Pre-allocate memory for this group.
    pshifts_treated <- matrix(NA, nrow = time_periods - g + 1, ncol = length(classes))
    rownames(pshifts_treated) <- paste0("t=", g:time_periods)
    colnames(pshifts_treated) <- paste0("Class_", classes)

    # Loop over "future" time periods and compute individual effects for units in this group at each time t > g.
    for (t in g:time_periods) {
      pshifts_treated[t - g + 1, ] <- colMeans(pmass1[units_in_group, ] - pmass0[units_in_group, ])
    }

    pshifts_treated_list[[paste0("G=", g)]] <- pshifts_treated
  }

  ## 6.) Output.
  unit_id <- rep(seq_len(n), time_periods)
  time <- rep(seq_len(time_periods), n)
  Y <- as.vector(Y)
  D <- rep(D, time_periods)
  X <- do.call(rbind, replicate(time_periods, X, simplify = FALSE))
  treatment_start[treatment_start == time_periods+1] <- 0 # Recode never-treated to have 0 as treatment date.

  return(list("Y" = Y, "D" = D, "unit_id" = unit_id, "time" = time, "treatment_start" = treatment_start, "X" = X, "pscore" = pscore, "pshifts_treated" = pshifts_treated_list))
}
