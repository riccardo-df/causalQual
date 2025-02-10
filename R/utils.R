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
