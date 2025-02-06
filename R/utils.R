#' Generate Qualitative Data (Selection-on-Observables)
#'
#' Generate a synthetic data set with a qualitative outcome, a binary treatment indicator, and a matrix of covariates. The treatment is either independent or conditionally independent of potential outcomes,
#' depending on users' choices.
#'
#' @param n Sample size.
#' @param assignment String controlling treatment assignment. Must be either \code{"randomized"} (random assignment) or \code{"observational"} (random assigment conditional on the generated covariates).
#' @param outcome_type String controlling the outcome type. Must be either \code{"multinomial"} or \code{"ordered"}. Affects how potential outcomes are constructed.
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
#' Potential outcomes are generated differently according to \code{outcome_type}. If \code{outcome_type == "multinomial"}, the function computes linear predictors for each class using the covariates:
#'
#' \deqn{\eta_{mi} (d) = \beta_{m1}^d X_{i1} + \beta_{m2}^d X_{i2} + \beta_{m3}^d X_{i3},}
#'
#' and then transforms \eqn{\eta_{mi} (d)} into valid probability distributions using the softmax function:
#'
#' \deqn{P(Y(d) = m \mid X) = \frac{\exp(\eta_{mi (d)})}{\sum_{m'} \exp(\eta_{m'i(d)})}.}
#'
#' It then generates potential outcomes \eqn{Y(1)} and \eqn{Y(0)} by sampling from \{1, 2, 3\} using \eqn{P(Y(d) = m \mid X)}.\cr
#'
#' If instead \code{outcome_type == "ordered"}, we first generate latent potential outcomes:
#'
#' \deqn{Y_i^* (d) = \tau d + X_{i1} + X_{i2} + X_{i3} + N (0, 1),}
#'
#' with \eqn{\tau = 2}. We then construct \eqn{Y_i (d)} by discretizing \eqn{Y_i^* (d)} using threshold parameters \eqn{\zeta_1 = 2} and \eqn{\zeta_2 = 4}. We then have
#'
#' \deqn{P(Y_i(d) = m) = P(\zeta_{m-1} < Y_i^*(d) \leq \zeta_m) = \Phi (\zeta_m - \sum_j X_{ij} - \tau d) - \Phi (\zeta_{m-1} - \sum_j X_{ij} - \tau d).}
#'
#' which allows us to analytically compute the probabilities of shift.
#'
#' ## Treatment assignment
#' Treatment is always assigned as \eqn{D_i \sim \text{Bernoulli}(\pi(X_i))}. If \code{assignment == "randomized"}, then the propensity score is specified as \eqn{P ( D_i = 1 | X_i)) = (X_{i1} + X_{i2}) / 2}.
#' If instead \code{assignment == "observational"}, then \eqn{P ( D_i = 1 | X_i)) = 0.5}.
#'
#' ## Other details
#' The function always generates three independent covariates from \eqn{U(0,1)}. Observed outcomes \eqn{Y_i} are always constructed using the usual observational rule.
#'
#' @importFrom stats rbinom
#' @importFrom stats runif
#'
#' @author Riccardo Di Francesco
#'
#' @seealso \code{\link{causalQual}}
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

    ## Probabilities of shift.
    pshifts <- colMeans(pmass1 - pmass0)
  } else if (outcome_type == "ordered") {
    ## Generate latent potential outcomes with ATE tau.
    tau <- 2
    Ystar1 <- tau + rowSums(X) + rnorm(n) # Essentially, coefficients set to 1.
    Ystar0 <- rowSums(X) + rnorm(n)

    ## Discretize to construct potential outcomes.
    zeta <- c(2, 4)

    Y1 <- cut(Ystar1, breaks = c(-Inf, zeta, Inf), labels = 1:3, right = TRUE)
    Y0 <- cut(Ystar0, breaks = c(-Inf, zeta, Inf), labels = 1:3, right = TRUE)

    Y0 <- as.integer(as.character(Y0))
    Y1 <- as.integer(as.character(Y1))

    ## Probabilities of shift.
    pmass1 <- cbind(pnorm(zeta[1] - rowSums(X) - tau),
                    pnorm(zeta[2] - rowSums(X) - tau) - pnorm(zeta[1] - rowSums(X) - tau),
                    1 - pnorm(zeta[2] - rowSums(X) - tau))

    pmass0 <- cbind(pnorm(zeta[1] - rowSums(X)),
                pnorm(zeta[2] - rowSums(X)) - pnorm(zeta[1] - rowSums(X)),
                1 - pnorm(zeta[2] - rowSums(X)))

    pshifts <- colMeans(pmass1 - pmass0)
  }

  ## 2.) Define propensity scores and assign treatment.
  if (assignment == "randomized") {
    pscore <- rep(0.5, n)
  } else if (assignment == "observational") {
    pscore <- (X[, 1] + X[, 3]) / 2
  }

  D <- stats::rbinom(n, 1, pscore)

  ## 3.) Generate observed outcomes.
  Y <- D * Y1 + (1 - D) * Y0

  ## 4.) Output.
  if (outcome_type == "multinomial") {
    pshifts <- colMeans(pmass1 - pmass0)
  } else if (outcome_type == "ordered") {
  }

  return(list("Y" = Y, "D" = D, "X" = X, "pscore" = pscore, "pshifts" = pshifts))
}


#' Softmax function.
#'
#' Implementation of the softmax function.
#'
#' @param logits An nxM matrix, with n the sample size and M the number of classes of an ordinal outcome.
#'
#' @return The same matrix where each row is normalized so that row sums equal 1; useful to construct conditional probability mass functions.
#'
#' @author Riccardo Di Francesco
#'
#' @keywords internal
softmax <- function(logits) {
  exp_logits <- exp(logits)
  return(exp_logits / rowSums(exp_logits))
}

#' Renaming Variables for LATEX Usage
#'
#' Renames variables where the character "_" is used, which causes clashes in LATEX. Useful for the \code{phased} print method.
#'
#' @param names string vector.
#'
#' @importFrom stringr str_split
#'
#' @return
#' The renamed string vector. Strings where "_" is not found are not modified by \code{rename_latex}.
#'
#' @keywords internal
rename_latex <- function(names) {
  ## Locating variables that need renaming.
  idx <- grepl("_", names, fixed = TRUE)

  if (sum(idx) == 0) return(names)

  ## Renaming variables.
  split_names <- stringr::str_split(string = names[idx], pattern = "_", simplify = TRUE)
  attach_names <- paste(split_names[, 1], split_names[, 2], sep = "\\_")

  ## Replacing.
  names[idx] <- attach_names

  ## Output.
  return(names)
}
