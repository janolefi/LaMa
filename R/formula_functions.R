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

