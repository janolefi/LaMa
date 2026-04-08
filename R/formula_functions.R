expand_cosinor <- function(formula) {
  
  # recursive, walk along the AST and replace any calls to cosinor() with the corresponding sin/cos terms
  .replace_cosinor <- function(lang) {
    if (!is.call(lang)) return(lang)
    
    if (identical(lang[[1]], quote(cosinor))) {
      # swap function name, evaluate to get string terms, rebuild as AST
      lang[[1]] <- quote(cosinor_int)
      terms <- eval(lang, envir = list(cosinor_int = cosinor_int))
      return(str2lang(paste(terms, collapse = " + ")))
    }
    
    # recurse into all sub-expressions
    for (i in seq_along(lang)) lang[[i]] <- .replace_cosinor(lang[[i]])
    lang
  }
  
  # Internal: returns character vector of term strings
  # Called via AST manipulation in expand_cosinor — not for direct use
  cosinor_int <- function(x, period = 24) {
    xname <- deparse(substitute(x))
    terms <- character(0)
    for (p in period) {
      terms <- c(terms,
                 sprintf("sin(2*pi*%s/%g)", xname, p),
                 sprintf("cos(2*pi*%s/%g)", xname, p))
    }
    terms
  }
  
  new_rhs <- .replace_cosinor(formula[[length(formula)]])
  
  if (length(formula) == 3) # has left side
    stats::as.formula(call("~", formula[[2]], new_rhs))
  else # doesn't have left side
    stats::as.formula(call("~", new_rhs))
}


#' Trigonometric basis expansion
#'
#' Builds a design matrix of \code{sin}/\code{cos} pairs for use in models
#' with periodic predictors. Can be used directly or inside formulas passed
#' to \code{make_matrices} (where expansion is handled automatically).
#'
#' The resulting columns form the basis for linear predictors of the form
#' \deqn{
#'   \eta_t = \beta_0 + \sum_k \Bigl(
#'     \beta_{1k} \sin\!\Bigl(\tfrac{2 \pi x_t}{\text{period}_k}\Bigr) +
#'     \beta_{2k} \cos\!\Bigl(\tfrac{2 \pi x_t}{\text{period}_k}\Bigr)
#'   \Bigr).
#' }
#'
#' @param x Numeric vector of the periodic variable.
#' @param period Numeric vector of period lengths, e.g. \code{24} for a daily cycle with hourly data or \code{c(24, 12)} for a daily + semi-daily cycle.
#'
#' @return A numeric matrix with \code{2 * length(period)} columns named
#'   \code{sin(2*pi*x/period)} / \code{cos(2*pi*x/period)}.
#' @export
#'
#' @examples
#' cosinor(1:24, period = 24)
#' cosinor(1:24, period = c(24, 12, 6))
#'
#' ## In model formulas (expand_cosinor handles the expansion):
#' form <- ~ x + temp * cosinor(hour, c(24, 12))
#' data <- data.frame(x = runif(24), temp = rnorm(24, 20), hour = 1:24)
#' modmat <- make_matrices(form, data = data)
cosinor <- function(x, period = 24) {
  xname <- deparse(substitute(x))
  out   <- matrix(NA_real_, nrow = length(x), ncol = 0)
  nms   <- character(0)
  for (p in period) {
    out <- cbind(out, sin(2*pi*x/p), cos(2*pi*x/p))
    nms <- c(nms,
             sprintf("sin(2*pi*%s/%g)", xname, p),
             sprintf("cos(2*pi*%s/%g)", xname, p))
  }
  colnames(out) <- nms
  out
}
