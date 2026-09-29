
###### Regression setting ------------------------------------------------------

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

## the shape of a LaMa_matrices object, in one place
new_LaMa_matrices <- function(Z, S, pardim, coef, sp0, formula, data, gam, gam0, knots) {
  out <- list(Z = Z, S = S, pardim = pardim, coef = coef, sp0 = sp0,
              formula = formula, data = data, gam = gam, gam0 = gam0, knots = knots)
  class(out) <- "LaMa_matrices"
  out
}

## combine per-component results into one object. Z, pardim, gam, gam0 and knots
## stay per component; S and coef are flattened, because the fitting code wants a
## single flat list of penalty matrices. 'prefix' qualifies their names with the
## component name, which is what distinguishes a nested formula list from a flat one.
combine_LaMa_matrices <- function(res, data, prefix = FALSE) {
  if(prefix) for(nm in names(res)){
    names(res[[nm]]$S) <- paste0(nm, ".", names(res[[nm]]$S))
    names(res[[nm]]$coef) <- paste0(nm, ".", names(res[[nm]]$coef))
    names(res[[nm]]$sp0) <- paste0(nm, ".", names(res[[nm]]$sp0))
  }

  # NULL components are dropped, because the accumulating loop this replaces used
  # res[[name]] <- NULL, which deletes rather than assigns. Matters for 'knots',
  # which is NULL whenever the caller supplied none
  per <- function(what) {
    v <- lapply(res, `[[`, what)
    keep <- !vapply(v, is.null, TRUE)
    if(any(keep)) v[keep] else list()
  }
  # unname() first: c() would otherwise qualify the names a second time
  flat <- function(what) do.call(c, unname(per(what)))

  new_LaMa_matrices(Z = per("Z"), S = flat("S"), pardim = per("pardim"),
                    coef = flat("coef"), sp0 = flat("sp0"),
                    formula = per("formula"), data = data,
                    gam = per("gam"), gam0 = per("gam0"), knots = per("knots"))
}

## response variable of a formula, used to name its block. NULL for a right-side
## only formula, unless an index is given to fall back on
get_name <- function(fml, idx = NULL) {
  if(length(fml) == 3 && !deparse(fml[[2]]) %in% c("", ".")) deparse(fml[[2]])
  else if(!is.null(idx)) paste0("par", idx)
  else NULL
}

make_matrices_flat <- function(formula, data, knots = NULL) {
  
  process_single <- function(fml, name, knots_sub) {
    fml_given <- fml # as the user wrote it, for print()
    fml <- expand_cosinor(fml)
    
    # prepare gam model setup
    gam_setup <- gam(update(fml, dummy ~ .), data = cbind(dummy = 1, data), 
                   knots = knots_sub, fit = FALSE)
    
    # also prepare prediction gam, only ever used by predict.gam() to build
    # prediction design matrices. gam(G = ) reuses the setup above instead of
    # constructing the whole basis a second time
    gam_setup0 <- gam(G = gam_setup, control = list(maxit = 1))
    
    Z <- gam_setup$X
    colnames(Z) <- gam_setup$term.names
    gam_setup$X <- NULL # returned as Z, and a second copy is the largest thing here
    
    # mgcv's initial guesses, one per penalty matrix and in the same order as
    # gam_setup$S, so the loop below can name them as it walks the smooths
    sp_init <- if(length(gam_setup$S)) {
      mgcv::initial.sp(Z, gam_setup$S, gam_setup$off)
    } else numeric(0)
    # started well above mgcv's guess: the inner fits are better behaved when
    # constrained initially, and the outer iteration relaxes the penalty from
    # there. Capped, so that a term whose basis carries a large covariate scale
    # does not start somewhere absurd
    sp_init <- pmin(100 * sp_init, 1e5)

    S <- list()
    coef <- list()
    sp0 <- numeric(0)
    pardim <- list(fixed_eff = gam_setup$nsdf)
    counter <- 1
    
    for (i in seq_along(gam_setup$smooth)) {
      sm <- gam_setup$smooth[[i]]
      if(is.null(name)){
        label <- sm$label
      } else {
        label <- paste0(name, ".", sm$label)
      }
      
      if (is.null(sm$margin)) {
        # Single-penalty smooth
        S[[label]] <- gam_setup$S[[counter]]
        pardim[[sm$label]] <- nrow(S[[label]])
        coef[[label]] <- rep(0, nrow(S[[label]]))
        sp0[label] <- sp_init[counter]
        counter <- counter + 1
      } else {
        # Tensor product smooth
        n_pen <- length(sm$margin)
        margin_names <- sapply(sm$margin, `[[`, "term")
        penalty_list <- gam_setup$S[counter:(counter + n_pen - 1)]
        names(penalty_list) <- margin_names
        S[[label]] <- penalty_list
        pardim[[sm$label]] <- nrow(penalty_list[[1]])
        coef[[label]] <- rep(0, nrow(penalty_list[[1]]))
        # one per margin, named <smooth>.<margin> to mirror the sub-list in S
        sp0[paste0(label, ".", margin_names)] <- sp_init[counter:(counter + n_pen - 1)]
        counter <- counter + n_pen
      }
    }
    
    list(
      Z = Z, S = S, pardim = pardim,
      gam = gam_setup, gam0 = gam_setup0,
      knots = knots_sub, coef = coef, sp0 = sp0, formula = fml_given
    )
  }
  
  if (!inherits(formula, "list")) {
    res <- process_single(formula, get_name(formula), knots) # NULL name -> unprefixed labels
    return(new_LaMa_matrices(res$Z, res$S, res$pardim, res$coef, res$sp0,
                             res$formula, data, res$gam, res$gam0, knots))
  }
  
  # get names
  form_names <- names(formula) # also allow for named lists with right-side-only formulas
  if(is.null(form_names)){
    form_names <- sapply(seq_along(formula), function(i) get_name(formula[[i]], i))
  }
  
  # check knots
  if(!is.null(knots)){
    if(!all(names(knots) %in% form_names)){
      stop("Names of 'knots' must match the names of the formulas in 'formula'.")
    }
  }
  
  # a list of formulas: the labels are already prefixed by process_single()
  res <- lapply(seq_along(formula), function(i)
    process_single(formula[[i]], form_names[i], knots[[form_names[i]]]))
  names(res) <- form_names
  combine_LaMa_matrices(res, data)
}

#' Build the design and the penalty matrix for models involving penalised splines based on a formula and a data set
#'
#' @param formula formula as used in \code{mgcv}. Formulas can be right-side only, or contain a response variable, which is just extracted for naming.
#'
#' Can also be a list of formulas, which are then processed separately. In that case, both a named list of right-side only formulas or a list of formulas with response variables can be provided.
#' @param data data frame containing all the variables on the right side of the formula(s)
#' @param knots optional list containing user specified knot values for each covariate to be used for basis construction.
#' For most bases the user simply supplies the \code{knots} to be used, which must match up with the \code{k} value supplied (note that the number of knots is not always just \code{k}).
#' See \code{mgcv} documentation for more details.
#'
#' If \code{formula} is a list, this needs to be a named (based on the response variables) list over such lists.
#' 
#' @seealso \code{\link{predict.LaMa_matrices}} for prediction design matrix construction based on the model matrices object created by this function.
#'
#' @return a list of class \code{LaMa_matrices} containing:
#' \item{\code{Z}}{design matrix (or list of such matrices if \code{formula} is a list))}
#' \item{\code{S}}{list of penalty matrices (with names based on the response terms of the formulas as well as the smooth terms and covariates). For tensorproduct smooths, corresponding entries are themselves lists, containing the \eqn{d} marginal penalty matrices if \eqn{d} is the dimension of the tensor product)}
#' \item{\code{pardim}}{list of parameter dimensions (fixed and penalised separately) for each formula, for ease of setting up initial parameters}
#' \item{\code{coef}}{list of coefficient vectors filled with zeros of the correct length for each formula, for ease of setting up initial parameters}
#' \item{\code{sp0}}{named vector of \strong{initial penalty strengths}, with one entry per penalty matrix and names matching \code{S} (a tensor product contributes one entry per margin, named \code{<smooth>.<margin>}).
#'
#' These are \code{\link[mgcv]{initial.sp}}'s relative guesses raised by a factor of 100 and capped at \code{1e5}. The factor is deliberate: for the models this package targets, the inner optimisation is better behaved when it starts constrained, and the outer iteration relaxes the penalty from there, whereas an under-penalised first fit can settle in a poor local optimum that later iterations inherit through the warm start. The cap keeps a term whose basis carries a large covariate scale, such as a \code{by =} smooth with a large covariate, from starting absurdly high.
#'
#' It can be passed straight to \code{\link{penalty}} or \code{\link{penalty2}} when each smooth appears once. 
#' When a smooth is replicated, e.g. one spline per state or per off-diagonal element of the transition probability matrix, the entries have to be repeated, and \code{lambda} runs over replicates \strong{within} each smooth. 
#' A single \code{rep()} over \code{sp0} is therefore only correct for some models: use \code{each = N} if every smooth carries one penalty, \code{times = N} if there is a single tensor product, and neither if the two are mixed. Repeating per smooth always works:
#' \preformatted{
#' n <- sapply(modmat$S, function(x) if(is.matrix(x)) 1 else length(x))
#' lambda <- unlist(lapply(split(modmat$sp0, rep(seq_along(n), n)),
#'                         rep, times = N))
#' }}
#' \item{\code{formula}}{the formula(s) as supplied, in the same nesting as \code{Z}}
#' \item{\code{data}}{the data frame used for the model(s)}
#' \item{\code{gam}}{unfitted \code{mgcv::gam} object used for construction of \code{Z} and \code{S} (or list of such objects if \code{formula} is a list). Its \code{X} element is not kept, as it is returned as \code{Z}}
#' \item{\code{gam0}}{fitted \code{mgcv::gam} which is used internally to create prediction design matrices (or list of such objects if \code{formula} is a list)}
#' \item{\code{knots}}{knot list used in the basis construction (or named list over such lists if \code{formula} is a list)}
#'
#' @export
#' 
#' @importFrom mgcv gam s initial.sp
#' @importFrom stats update
#'
#' @examples
#' data = data.frame(x = runif(100), 
#'                   y = runif(100),
#'                   g = factor(rep(1:10, each = 10)))
#'
#' # unvariate thin plate regression spline
#' modmat = make_matrices(~ s(x), data)
#' # univariate P-spline
#' modmat = make_matrices(~ s(x, bs = "ps"), data)
#' # adding random intercept
#' modmat = make_matrices(~ s(g, bs = "re") + s(x, bs = "ps"), data)
#' # tensorproduct of x and y
#' modmat = make_matrices(~ s(x) + s(y) + ti(x,y), data)
#' # multiple formulas at once
#' modmat = make_matrices(list(mu ~ s(x) + y, sigma ~ s(g, bs = "re")), data = data)
make_matrices <- function(formula, data, knots = NULL){
  ## 0 = a single formula, 1 = a list of formulas, 2 = a list of such lists
  is_formula_list <- function(x){
    is.list(x) && length(x) > 0 && all(vapply(x, inherits, TRUE, what = "formula"))
  }
  depth <- if(inherits(formula, "formula")) 0L
           else if(is_formula_list(formula)) 1L
           else if(is.list(formula) && length(formula) > 0 &&
                   all(vapply(formula, is_formula_list, TRUE))) 2L
           else stop("'formula' must be a (nested) list of formulas or a single formula.")

  if(depth == 2L){
    names_out <- names(formula)
    if(is.null(names_out)) names_out <- paste0("stream", seq_along(formula))
    res <- lapply(seq_along(formula),
                  function(i) make_matrices_flat(formula[[i]], data, knots[[i]]))
    names(res) <- names_out
    # prefix: S and coef of every stream are qualified by the stream name
    return(combine_LaMa_matrices(res, data, prefix = TRUE))
  }

  if(depth == 1L && !is.null(knots)){
    if(!is.list(knots)){
      stop("'knots' must be a list of knots for each formula in 'formula'.")
    } else if(!all(sapply(knots, is.list))){
      stop("'knots' must be a named list of knot lists, containing knots for the corresponding formula (top level) and variable (bottom level) in 'formula'.")
    }
  }
  make_matrices_flat(formula, data, knots)
}

#' Repeat initial penalty strengths for replicated smooths
#'
#' @description
#' Expands the \code{sp0} element of a \code{\link{make_matrices}} object into a
#' \code{lambda} vector for \code{\link{penalty}} or \code{\link{penalty2}}, when each smooth is
#' used more than once, for example one spline per state or per off-diagonal element of the
#' transition probability matrix.
#'
#' @details
#' \code{lambda} runs over replicates \strong{within} each smooth, so for a tensor product the
#' margins vary fastest. A single \code{rep()} over \code{sp0} therefore only gives the right
#' answer for some models, silently producing the right length but the wrong values for others.
#' This function repeats each smooth's penalty strengths separately, which is always correct.
#'
#' @param model_matrices object of class \code{LaMa_matrices} as returned by \code{\link{make_matrices}}
#' @param N number of replicates. Either a single number used for every smooth, or one number per
#' smooth, i.e. per element of \code{model_matrices$S}
#'
#' @return named numeric vector of initial penalty strengths, of the length and in the order
#' expected for \code{lambda}
#'
#' @seealso \code{\link{make_matrices}}, which computes \code{sp0}
#' @export
#'
#' @examples
#' modmat = make_matrices(~ s(x) + ti(x, y), data.frame(x = runif(100), y = runif(100)))
#' modmat$sp0         # one value per penalty
#' rep_sp(modmat, 2)  # ready to use as 'lambda' for two replicates of each smooth
rep_sp <- function(model_matrices, N) {
  if(!inherits(model_matrices, "LaMa_matrices")){
    stop("'model_matrices' needs to be an object returned by 'make_matrices()'")
  }
  sp0 <- model_matrices$sp0
  n_pen <- vapply(model_matrices$S, function(s) if(is.matrix(s)) 1L else length(s), 1L)
  if(length(n_pen) == 0) return(numeric(0))
  if(length(N) == 1) N <- rep(N, length(n_pen))
  if(length(N) != length(n_pen)){
    stop("'N' needs to be a single number or one number per smooth, i.e. of length ",
         length(n_pen))
  }
  if(any(N < 1) || any(N != round(N))) stop("'N' needs to contain positive whole numbers")

  grp <- rep(seq_along(n_pen), n_pen) # which smooth each penalty belongs to
  # unnamed list, so that unlist() does not qualify the names a second time
  unlist(lapply(seq_along(n_pen), function(i) {
    v <- sp0[grp == i]
    stats::setNames(rep(v, times = N[i]),
                    paste0(names(v), ".", rep(seq_len(N[i]), each = length(v))))
  }))
}


#' Print a \code{LaMa_matrices} object
#'
#' @param x object of class \code{LaMa_matrices} as returned by \code{\link{make_matrices}}
#' @param ... ignored
#'
#' @return \code{x}, invisibly. Called for its side effect of printing.
#'
#' @seealso \code{\link{make_matrices}}
#' @export
print.LaMa_matrices <- function(x, ...) {
  # flatten the formula tree into rows: a stream name has no right-hand side,
  # a formula has both, and depth drives the indentation
  rows <- list()
  walk <- function(f, name, depth) {
    if(inherits(f, "formula")){
      has_lhs <- length(f) == 3
      lhs <- if(has_lhs) deparse(f[[2]]) else if(!is.null(name)) name else ""
      rhs <- paste(deparse(f[[length(f)]]), collapse = " ")
      rows[[length(rows) + 1]] <<- list(depth = depth, lhs = lhs,
                                        rhs = gsub("\\s+", " ", rhs))
    } else {
      nms <- names(f)
      for(i in seq_along(f)){
        nm <- if(is.null(nms)) NULL else nms[i]
        if(inherits(f[[i]], "formula")){
          walk(f[[i]], nm, depth)
        } else { # a further level of nesting, e.g. one entry per data stream
          rows[[length(rows) + 1]] <<- list(depth = depth,
                                            lhs = if(is.null(nm)) paste0("stream", i) else nm,
                                            rhs = NULL)
          walk(f[[i]], NULL, depth + 1)
        }
      }
    }
  }
  walk(x$formula, NULL, 0)

  cat("Model matrices for\n\n")
  fml_rows <- vapply(rows, function(r) !is.null(r$rhs), TRUE)
  # line the tildes up, per indentation level
  width <- stats::setNames(rep(0L, 10), as.character(0:9))
  for(r in rows[fml_rows]) width[as.character(r$depth)] <-
    max(width[as.character(r$depth)], nchar(r$lhs))
  for(r in rows){
    pad <- strrep("  ", r$depth + 1)
    if(is.null(r$rhs)) cat(pad, r$lhs, "\n", sep = "")
    else if(!nzchar(r$lhs)) cat(pad, "~ ", r$rhs, "\n", sep = "") # right-side only
    else cat(pad, formatC(r$lhs, width = width[as.character(r$depth)], flag = "-"),
             " ~ ", r$rhs, "\n", sep = "")
  }

  cat("\n---\n\n")
  if(length(x$sp0) == 0){
    cat("No penalised terms, so no initial penalty strengths.\n")
  } else {
    cat("Initial penalty strengths:\n\n")
    print(signif(x$sp0, 4))
    cat("\nUse 'rep_sp()' if a smooth is used more than once.\n")
  }
  invisible(x)
}


#' Process and standardise formulas for the state process of hidden Markov models
#'
#' @param formulas formulas for the transition process of a hidden Markov model, either as a single formula, a list of formulas, or a matrix.
#' @param nStates number of states of the Markov chain
#' @param ref optional vector of reference categories for each state, defaults to \code{1:nStates}. 
#' If provided, must be of length \code{nStates} and contain valid state indices. 
#' If a formula matrix is provided, this cannot be specified because reference categries are specified by one \code{"."} entry in each row.
#'
#' @returns named list of formulas of length \code{nStates * (nStates - 1)}, where each formula corresponds to a transition from state \eqn{i} to state \eqn{j}, excluding transitions from \code{i} to \code{ref[i]}.
#' @export
#'
#' @examples
#' # single formula for all non-reference category elements
#' formulas = process_hid_formulas(~ s(x), nStates = 3)
#' # now a list of length 6 with names tr.ij, not including reference categories
#' 
#' # different reference categories
#' formulas = process_hid_formulas(~ s(x), nStates = 3, ref = c(1,1,1))
#' 
#' # different formulas for different entries (and only for 2 of 6)
#' formulas = list(tr.12 ~ s(x), tr.23 ~ s(y))
#' formulas = process_hid_formulas(formulas, nStates = 3, ref = c(1,1,1))
#' # also a list of length 6, remaining entries filled with tr.ij ~ 1
#' 
#' # matrix input with reference categories
#' formulas = matrix(c(".", "~ s(x)", "~ s(y)",
#'                     "~ g", ".", "~ I(x^2)",
#'                     "~ y", "~ 1", "."), 
#'                     nrow = 3, byrow = TRUE)
#' # dots define reference categories
#' formulas = process_hid_formulas(formulas, nStates = 3)
process_hid_formulas <- function(formulas, 
                                 nStates,
                                 ref = NULL){
  
  # extract formula names from response variable
  get_name <- function(fml) {
    if (length(fml) == 3 && !deparse(fml[[2]]) %in% c("", ".")) {
      return(deparse(fml[[2]]))
    } else {
      return(NULL)
    }
  }
  
  # initialise transition names without reference categories
  make_tr_names <- function(nStates, ref){
    names <- c()
    for(i in 1:nStates){
      for(j in 1:nStates){
        if(j != ref[i]){
          names <- c(names, paste0("tr.", i, j))
        }
      }
    }
    names
  }
  
  # NULL input → uniform ~1 formulas
  if(is.null(formulas)){
    formulas <- ~ 1
  }
  
  # Matrix input: extract ref directly from "."
  if(inherits(formulas, "matrix")){
    if (!all(dim(formulas) == c(nStates, nStates))) {
      stop("Formula matrix must be of size nStates x nStates.")
    }
    
    ref <- integer(nStates)
    for (i in 1:nStates) {
      row <- formulas[i, ]
      dot_count <- sum(trimws(row) == "." | trimws(row) == "~ .")
      if (dot_count != 1) {
        stop(sprintf("Row %d of formula matrix must contain exactly one '.' (reference transition), found %d.", i, dot_count))
      }
      ref[i] <- which(trimws(row) %in% c(".", "~ ."))
    }
  }
  
  # Default reference if not yet defined
  if(is.null(ref)){
    ref <- 1:nStates
  } else{
    if (length(ref) != nStates) stop("'ref' must be of length nStates.")
    if (!all(ref %in% 1:nStates)) stop("'ref' must contain valid state indices.")
  }
  
  # Generate transition names and empty formula list
  names_out <- make_tr_names(nStates, ref)
  formula_list <- stats::setNames(vector("list", length(names_out)), names_out)
  
  # Fill from matrix (now that ref is known)
  if(inherits(formulas, "matrix")){
    for (i in 1:nStates) {
      for (j in 1:nStates) {
        entry <- formulas[i, j]
        if (trimws(entry) %in% c(".", "~ .")) next
        if (j == ref[i]) next
        
        name <- paste0("tr.", i, j)
        if (is.character(entry)) entry <- stats::as.formula(entry)
        formula_list[[name]] <- entry
      }
    }
    return(formula_list)
  }
  
  # Uniform formula for all transitions
  if(inherits(formulas, "formula")){
    for(i in seq_along(formula_list)){
      formula_list[[i]] <- formulas
    }
    return(formula_list)
  }
  
  # Named or unnamed list
  if(inherits(formulas, "list")){
    form_names <- names(formulas)
    if(is.null(form_names)){
      form_names <- sapply(formulas, get_name)
    }
    if(is.list(form_names)){
      stop("If formula list is provided, it must be a named list with names 'tr.ij'.")
    }
    if(!all(form_names %in% names_out)){
      msg <- "A formula for a reference category was provided. Formulas should not be provided for:"
      trs <- paste0("tr.", 1:nStates, ref, collapse = ", ")
      msg <- paste(msg, trs)
      stop(msg)
    }
    
    for (n in names_out) formula_list[[n]] <- ~1
    for (i in seq_along(formulas)) {
      formula_list[[form_names[i]]] <- formulas[[i]]
    }
    return(formula_list)
  }
  
  stop("Unsupported input type for 'formulas'. Must be a formula, list, or matrix.")
}





#' Build the prediction design matrix based on new data and model_matrices object created by \code{\link{make_matrices}}
#'
#' @param object model matrices object as returned from \code{\link{make_matrices}}
#' @param newdata data frame containing the variables in the formula and new data for which to evaluate the basis
#' @param what optional character string specifying which formula to use for prediction if \code{object} contains multiple formulas.
#' @param ... passed on to \code{\link{pred_matrix}}, in particular \code{exclude} to set terms to zero in the predicted design matrix
#'
#' @seealso \code{\link{make_matrices}} for creating objects of class \code{LaMa_matrices} which can be used for prediction by this function.
#'
#' @return prediction design matrix for \code{newdata} with the same basis as used for \code{model_matrices}
#' @export
#' 
#'
#' @examples
#' # single formula
#' modmat = make_matrices(~ s(x), data.frame(x = 1:10))
#' Z_p = predict(modmat, data.frame(x = 1:10 - 0.5))
#' # with multiple formulas
#' modmat = make_matrices(list(mu ~ s(x), sigma ~ s(x, bs = "ps")), data = data.frame(x = 1:10))
#' Z_p = predict(modmat, data.frame(x = 1:10 - 0.5), what = "mu")
#' # nested formula list
#' form = list(stream1 = list(mu ~ s(x), sigma ~ s(x, bs = "ps")))
#' modmat = make_matrices(form, data = data.frame(x = 1:10))
#' Z_p = predict(modmat, data.frame(x = 1:10 - 0.5), what = c("stream1", "mu"))
predict.LaMa_matrices <- function(object, newdata, what = NULL, ...){
  pred_matrix_internal(object, newdata = newdata, what = what, ...) # ... carries 'exclude'
}


#' @importFrom mgcv gam
#' @importFrom mgcv s
#' @importFrom mgcv predict.gam
#' @noRd
## the worker behind predict.LaMa_matrices(); pred_matrix() is a deprecated
## wrapper around it, so the deprecation warning is not raised by predict()
pred_matrix_internal <- function(model_matrices, newdata, what = NULL, exclude = NULL) {
  
  if(is.null(model_matrices$gam0)){
    stop("'model_matrices' contains no 'gam0' object; rebuild it with make_matrices()")
  } else {
    if(inherits(model_matrices$gam0, "gam")){
      if(!is.null(what)){
        message("'model_matrics' only contains a single gam object, 'what' is not used.")
      }
      gam_setup0 <- model_matrices$gam0
    } else {
      if(!inherits(model_matrices$gam0[[1]], "gam")){ # nested list of formulas
        if(is.null(what)){
          stop("'what' must be specified and contain a top-level name of the formula list and a name/ response variable from your original formulas.")
        }
        if(length(what) != 2 | !is.character(what)){
          stop("'what' must be a character vector of length 2, containing the top-level name of the formula list and a name/ response variable from your original formulas.")
        }
        if(!what[1] %in% names(model_matrices$gam0)){
          stop("'what[1]' must be a top-level name of your formula list.")
        }
        if(!what[2] %in% names(model_matrices$gam0[[what[1]]])){
          stop("'what[2]' must be one of the names/ response variables of your original formulas.")
        }
        gam_setup0 <- model_matrices$gam0[[what[1]]][[what[2]]]
      } else { # flat list of formulas
        if(is.null(what)){
          stop("'what' must be specified and be one of the names/ response variables from your original formulas.")
        }
        if(!what %in% names(model_matrices$gam0)){
          stop("'what' must be one of the names of your original formulas.")
        }
        gam_setup0 <- model_matrices$gam0[[what]]
      }
    }
  }
  
  # if pi in formula, mgcv requires pi column in data
  if("pi" %in% all.names(gam_setup0$formula)){
    newdata = cbind(newdata, pi)
  }
  
  predict.gam(gam_setup0, 
              newdata = cbind(dummy = 1, newdata), 
              type = "lpmatrix",
              exclude = exclude)
}

#' Build the prediction design matrix based on new data and model_matrices object created by \code{\link{make_matrices}}
#'
#' @description
#' \strong{Deprecated.} Use \code{\link[=predict.LaMa_matrices]{predict}} on the object returned by \code{\link{make_matrices}} instead, which does the same thing:
#' \code{predict(modmat, newdata)} in place of \code{pred_matrix(modmat, newdata)}.
#'
#' @param model_matrices model_matrices object as returned from \code{\link{make_matrices}}
#' @param newdata data frame containing the variables in the formula and new data for which to evaluate the basis
#' @param what optional character string specifying which formula to use for prediction, if \code{object} contains multiple formulas. If \code{NULL}, the first formula is used.
#' @param exclude optional vector of terms to set to zero in the predicted design matrix. Useful for predicting main effects only when e.g. \code{sd(..., bs = "re")} terms are present. See \code{mgcv::predict.gam} for more details.
#' @return prediction design matrix for \code{newdata} with the same basis as used for \code{model_matrices}
#'
#' @seealso \code{\link[=predict.LaMa_matrices]{predict}}, which replaces this function
#' @export
#'
#' @examples
#' modmat = make_matrices(~ s(x), data.frame(x = 1:10))
#' # deprecated:
#' # Z_p = pred_matrix(modmat, data.frame(x = 1:10 - 0.5))
#' # use instead:
#' Z_p = predict(modmat, data.frame(x = 1:10 - 0.5))
pred_matrix = function(model_matrices, 
                       newdata,
                       what = NULL,
                       exclude = NULL) {
  .Deprecated("predict", package = "LaMa",
              msg = paste("'pred_matrix()' is deprecated; use 'predict()' on the object",
                          "returned by 'make_matrices()' instead."))
  pred_matrix_internal(model_matrices, newdata = newdata, what = what, exclude = exclude)
}



###### Density estimation setting ----------------------------------------------

#' Build a standardised P-Spline design matrix and the associated P-Spline penalty matrix
#' 
#' This function builds the B-spline design matrix for a given data vector. 
#' Importantly, the B-spline basis functions are normalised such that the integral of each basis function is 1, hence this basis can be used for spline-based density estimation, when the basis functions are weighted by non-negative weights summing to one.
#'
#' @param x data vector
#' @param k number of basis functions
#' @param type type of the data, either \code{"real"} for data on the reals, \code{"positive"} for data on the positive reals or \code{"circular"} for circular data like angles.
#' @param degree degree of the B-spline basis functions, defaults to cubic B-splines
#' @param knots optional vector of knots (including the boundary knots) to be used for basis construction. 
#' If not provided, the knots are placed equidistantly for \code{"real"} and \code{"circular"} and using polynomial spacing for \code{"positive"}.
#'
#' For \code{"real"} and \code{"positive"} \code{k - degree + 1} knots are needed, for \code{"circular"} \code{k + 1} knots are needed.
#' # @param quantile logical, if \code{TRUE} use quantile-based knot spacing (instead of equidistant or polynomial)
#' @param diff_order order of differencing used for the P-Spline penalty matrix for each data stream. Defaults to second-order differences.
#' @param pow power for polynomial knot spacing
#' @param npoints number of points used in the numerical integration for normalizing the B-spline basis functions
#' 
#' Such non-equidistant knot spacing is only used for \code{type = "positive"}.
#'
#' @return list containing the design matrix \code{Z}, the penalty matrix \code{S}, the prediction design matrix \code{Z_predict}, the prediction grid \code{xseq}, and details for the basis expansion.
#' @export
#' @importFrom splines2 bSpline
#' @importFrom mgcv cSplineDes
#'
#' @examples
#' set.seed(1)
#' # real-valued
#' x <- rnorm(100)
#' modmat <- make_matrices_dens(x, k = 20)
#' # positive-continuouos
#' x <- rgamma2(100, mean = 5, sd = 2)
#' modmat <- make_matrices_dens(x, k = 20, type = "positive")
#' # circular
#' x <- rvm(100, mu = 0, kappa = 2)
#' modmat <- make_matrices_dens(x, k = 20, type = "circular")
#' # bounded in an interval
#' x <- rbeta(100, 1, 2)
#' modmat <- make_matrices_dens(x, k = 20)
make_matrices_dens = function(x, # data vector
                              k, # number of basis functions
                              type = "real", # type of the data
                              degree = 3, # degree of the B-Spline basis
                              knots = NULL, # default to automatic knots spacing, if provided, need to be k - degree + 1
                              # quantile = FALSE, # if TRUE, use quantile-based knots
                              diff_order = 2, # order of the differences for the penalty matrix
                              pow = 0.5, # power for polynomial knot spacing for positive values
                              npoints = 1e4 # number of points for numerical integration
){
  nObs <- length(x)
  
  quantile <- FALSE # no quantile spacing if knots are not custom

  if(type != "circular"){
    
    # getting the data range
    rangex <- range(x, na.rm = TRUE)
    # slightly inflating the range
    rangex <- rangex + c(-1, 1) * diff(rangex) / 20
    # computing the number of knots
    nrknots <- k - degree + 1 
    
    if(type == "real"){ # real data, defaults to equidistant knots
      # grid for numerical integration (normalisation)
      # xseq <- seq(rangex[1] + 1e-3, rangex[2] - 1e-3, length = npoints)
      
      # if knots not supplied, default to equidistant knots
      if(is.null(knots)){
        if(quantile){ # quantile spacing
          knots <- quantile(x, probs = seq(0, 1, length = nrknots), na.rm = TRUE)
          knots[1] <- rangex[1]
          knots[nrknots] <- rangex[2]
        } else { # equidistant spacing
          knots <- seq(rangex[1], rangex[2], length = nrknots)
        }
      } else{
        if(length(knots) != nrknots){
          stop("Number of knots provided is wrong. It should be 'k - degree + 1'.")
        }
      }
      
      boundary_knots <- knots[c(1, length(knots))] # set boundary knots
      knots <- knots[2:(length(knots)-1)] # only keep interior knots
      
      xseq <- seq(boundary_knots[1] + 1e-3, boundary_knots[2] - 1e-3, 
                  length = npoints)
      
      # numerical integration for normalizing the B-spline basis functions
      B0 <- bSpline(
        x = xseq, 
        knots = knots, 
        Boundary.knots = boundary_knots,
        degree = degree, 
        intercept = TRUE
        ) # unnormalised spline design matrix

      # interval width for numerical integration
      h <- xseq[2] - xseq[1]
      # numerical integration of each basis function
      w <- (h * colSums(B0))^-1
      
      # normalising
      B0 <- t(t(B0) * w)
      
      # actual data design matrix
      B <- matrix(NA, nrow = nObs, ncol = k)
      ind <- which(!is.na(x))
      B[ind,] <- bSpline(
        x = x[ind], 
        knots = knots, 
        Boundary.knots = boundary_knots,
        degree = degree, 
        intercept = TRUE
      ) # unnormalised spline design matrix
      
      # normalising
      B[ind,] <- t(t(B[ind,]) * w)
 
      # basis positions for initial values later (expected value of each basis function)
      basis_pos = colSums(xseq * t(t(B0) / rowSums(t(B0))))
      # basis_pos = knots[(degree):(length(knots)-degree+1)]
      
      # second-order difference matrix
      L <- diff(diag(k), differences = diff_order)

    } else if(type == "positive") { # non-equidistant knots, no mass on < 0
      
      if(min(x, na.rm = TRUE) <= 0){
        stop("When type = 'positive', 'x' can only contain positive values")
      } 
      
      # square-root spacing
      rangex[1] <- 0 # always from 0 to max(x)
      
      # grid for numerical integration (normalisation)
      # xseq <- seq(rangex[1] + 1e-3, rangex[2] - 1e-3, length = npoints)
      
      # if knots not supplied, default to square-root-spaced knots
      if(is.null(knots)){
        if(quantile){ # quantile spacing
          knots <- quantile(x, probs = seq(0, 1, length = nrknots), na.rm = TRUE)
          knots[1] <- rangex[1]
          knots[nrknots] <- rangex[2]
        } else { # polynomial spacing
          knots <- seq(rangex[1]^pow, rangex[2]^pow, length = nrknots)^(1/pow)
        }
      } else{
        if(length(knots) != nrknots){
          stop("Number of knots provided is wrong. It should be 'k - degree + 1'.")
        }
      }
      
      boundary_knots <- knots[c(1, length(knots))] # set boundary knots
      knots <- knots[2:(length(knots)-1)] # only keep interior knots
      
      xseq <- seq(boundary_knots[1] + 1e-3, boundary_knots[2] - 1e-3, 
                  length = npoints)
      
      # numerical integration for normalizing the B-spline basis functions
      B0 <- bSpline(
        x = xseq, 
        knots = knots, 
        Boundary.knots = boundary_knots,
        degree = degree, 
        intercept = TRUE
      ) # unnormalised spline design matrix
      
      # interval width for numerical integration
      h <- xseq[2] - xseq[1]
      
      # numerical integration of each basis function
      w <- (h * colSums(B0))^-1
      
      # normalising
      B0 <- t(t(B0) * w)
      
      # basis positions for initial values later (expected value of each basis function)
      basis_pos <- colSums(xseq * t(t(B0) / rowSums(t(B0))))
      
      B <- matrix(NA, nrow = nObs, ncol = k)
      ind <- which(!is.na(x))
      B[ind,] <- bSpline(
        x = x[ind], 
        knots = knots, 
        Boundary.knots = boundary_knots,
        degree = degree, 
        intercept = TRUE
      ) # unnormalised spline design matrix
      
      # normalising
      B[ind,] <- t(t(B[ind,]) * w)
      
      # second-order difference matrix
      L <- diff(diag(k), differences = diff_order) 
      
    } else {
      stop("type must be 'real', 'positive' or 'circular'")
    }
    
  } else if(type == "circular"){
    
    # if knots not supplied, default to equidistant knots
    if(is.null(knots)){
      if(quantile){ # quantile spacing
        knots <- quantile(x, probs = seq(0, 1, length = k+1), na.rm = TRUE)
        knots[1] <- -pi
        knots[length(knots)] <- pi
      } else { # equidistant spacing
        knots <- seq(-pi, pi, length = k + 1)
      }
    } else{
      if(length(knots) != nrknots){
        stop("Number of knots provided is wrong. For 'type = circular' it should be 'k + 1' with '-pi' and 'pi' as first and last entries.")
      }
    }
    
    xseq <- seq(-pi, pi, length = npoints)
    B0 <- cSplineDes(
      xseq, 
      knots, 
      ord = degree + 1
      )
    
    # interval width for numerical integration
    h <- xseq[2] - xseq[1]
    
    # numerical integration of each basis function
    w <- (h * colSums(B0))^-1
    
    # normalising
    B0 <- t(t(B0) * w)
    
    # basis positions for initial values later (expected value of each basis function)
    basis_pos <- colSums(xseq * t(t(B0) / rowSums(t(B0))))
    
    # actual data design matrix
    B <- matrix(NA, nrow = nObs, ncol = k)
    ind <- which(!is.na(x))
    B[ind,] <- cSplineDes(
      x[ind], 
      knots, 
      ord = degree + 1
    )
    # normalising
    B[ind,] <- t(t(B[ind,]) * w)
    
    # difference matrix for circular P-Spline penalty
    L <- diff(rbind(diag(k), diag(k)[1:diff_order,]), differences = diff_order)
  }
  
  # constructing penalty matrix
  S <- crossprod(L[,-k], L[,-k]) # leaving out last column because parameter set to zero
  message("Leaving out last column of the penalty matrix, fix the last spline coefficient at zero for identifiability!")
  
  basis <- list(
    type = type, 
    knots = knots, 
    w = w, 
    degree = degree, 
    basis_pos = basis_pos
    )
  
  out <- list(
    Z = B,
    S = S,
    Z_predict = B0[(1:500) * (npoints/500) -(npoints/500)/2,],
    xseq = xseq[(1:500)* (npoints/500) -(npoints/500)/2],
    basis = basis
  )
  
  return(out)
}

# helper function, not exported
make_splinecoef = function(model_matrices, 
                           type = "real", 
                           par){
  basis_pos = model_matrices$basis$basis_pos
  k = length(basis_pos)
  
  # scaling
  dens1 = rowSums(model_matrices$Z_predict)
  scaling = sapply(1:k, function(i) dens1[which.min(abs(basis_pos[i] - model_matrices$xseq))])
  
  if(type == "real"){ # if density has support on the reals -> use normal distribution to initialize
    beta = sapply(basis_pos, dnorm, mean = par$mean, sd = par$sd, log = TRUE)
    if(is.vector(beta)) beta = matrix(beta, nrow = 1)
    # rescaling to account for non-equidistant knot spacing
    beta = t(t(beta) - log(apply(model_matrices$Z_predict, 2, max)))
  } else if(type == "positive") { # if density has support on positive continuous -> use gamma distribution
    # transformation to scale and shape
    beta = sapply(basis_pos, dgamma2, mean = par$mean, sd = par$sd, log = TRUE)
    if(is.vector(beta)) beta = matrix(beta, nrow = 1)
    # rescaling to account for non-equidistant knot spacing
    beta = t(t(beta) - log(apply(model_matrices$Z_predict, 2, max)))
  } else if(type == "circular") {
    beta = sapply(basis_pos, LaMa::dvm, mu = par$mean, kappa = par$concentration, log = TRUE)
    if(is.vector(beta)) beta = matrix(beta, nrow = 1)
  }

  beta = beta - beta[,k]
  # beta = beta - beta[,k-1]
  message("Parameter matrix excludes the last column. Add zero column using 'cbind(coef, 0)' in your loss function!\n")
  return(beta[,-k])
}

#' Build the design and penalty matrices for smooth density estimation
#' 
#' @description
#' This high-level function can be used to prepare objects needed to estimate mixture models of smooth densities using P-Splines.
#' 
#' @details
#' Under the hood, \code{\link{make_matrices_dens}} is used for the actual construction of the design and penalty matrices.
#'
#' You can provide one or multiple data streams of different types (real, positive, circular) and specify initial means and standard deviations/ concentrations for each data stream. This information is then converted into suitable spline coefficients.
#' \code{smooth_dens_construct} then constructs the design and penalty matrices for standardised B-splines basis functions (integrating to one) for each data stream.
#' For types \code{"real"} and \code{"circular"} the knots are placed equidistant in the range of the data, for type \code{"positive"} the knots are placed using polynomial spacing.
#'
#' @param data named data frame of 1 or multiple data streams
#' @param par nested named list of initial means and sds/concentrations for each data stream
#' @param type vector of length 1 or number of data streams containing the type of each data stream, either \code{"real"} for data on the reals, \code{"positive"} for data on the positive reals or \code{"circular"} for angular data.
#' @param k vector of length 1 or number of data streams containing the number of basis functions for each data stream
#' @param knots optional list of knots vectors (including the boundary knots) to be used for basis construction. 
#' If not provided, the knots are placed equidistantly for \code{"real"} and \code{"circular"} and using polynomial spacing for \code{"positive"}.
#'
#' For \code{"real"} and \code{"positive"} \code{k - degree + 1} knots are needed, for \code{"circular"} \code{k + 1} knots are needed.
#' @param degree degree of the B-spline basis functions for each data stream, defaults to cubic B-splines
#' @param diff_order order of differencing used for the P-Spline penalty matrix for each data stream. Defaults to second-order differences.
#'
#' @return a nested list containing the design matrices \code{Z}, the penalty matrices \code{S}, the initial coefficients \code{coef} the prediction design matrices \code{Z_predict}, the prediction grids \code{xseq}, and details for the basis expansion for each data stream.
#' @export
#'
#' @examples
#' ## 3 data streams, each with one distribution
#' # normal data with mean 0 and sd 1
#' x1 = rnorm(100, mean = 0, sd = 1)
#' # gamma data with mean 5 and sd 3
#' x2 = rgamma2(100, mean = 5, sd = 3)
#' # circular data
#' x3 = rvm(100, mu = 0, kappa = 2)
#' 
#' data = data.frame(x1 = x1, x2 = x2, x3 = x3)
#' 
#' par = list(x1 = list(mean = 0, sd = 1),
#'            x2 = list(mean = 5, sd = 3),
#'            x3 = list(mean = 0, concentration = 2))
#' 
#' SmoothDens = smooth_dens_construct(data, 
#'                                    par,
#'                                    type = c("real", "positive", "circular"))
#'                              
#' # extracting objects for x1
#' Z1 = SmoothDens$Z$x1
#' S1 = SmoothDens$S$x1
#' coefs1 = SmoothDens$coef$x1
#' 
#' ## one data stream, but mixture of two distributions
#' # normal data with mean 0 and sd 1
#' x = rnorm(100, mean = 0, sd = 1)
#' data = data.frame(x = x)
#' 
#' # now parameters for mixture of two normals
#' par = list(x = list(mean = c(0, 5), sd = c(1,1)))
#' 
#' SmoothDens = smooth_dens_construct(data, par = par)
#' 
#' # extracting objects 
#' Z = SmoothDens$Z$x
#' S = SmoothDens$S$x
#' coefs = SmoothDens$coef$x
smooth_dens_construct <- function(data,
                                  par,
                                  type = "real",
                                  k = 25,
                                  knots = NULL,
                                  # quantile = FALSE,
                                  degree = 3,
                                  diff_order = 2
){
  quantile <- FALSE # no quantile spacing if knots are not custom
  
  if(!is.data.frame(data)){
    stop("data must be a data frame")
  }
  
  nStreams <- ncol(data)
  nObs <- nrow(data)
  varnames <- colnames(data)
  
  if(length(par) == length(varnames)){
    if(is.null(names(par))){
      stop("'par' must be a named nested list with names corresponding to the datastreams")
    } else if(any(names(par) != varnames)){
      stop("'par' must be a named nested list with names corresponding to the datastreams")
    }
  } else{
    stop("'par' must be a named nested list with names corresponding to the datastreams")
  }
  
  # processing input arguments
  if(length(k) == 1){
    k = rep(k, nStreams)
  } else if(length(k) != nStreams){
    stop("'k' must be a scalar or a vector of length equal to the number of datastreams")
  }
  if(length(type) == 1){
    type = rep(type, nStreams)
  } else if(length(type) != nStreams){
    stop("'type' must be of length 1 or equal to the number of datastreams")
  }
  if(length(degree) == 1){
    degree = rep(degree, nStreams)
  } else if(length(degree) != nStreams){
    stop("'degree' must be of length 1 or equal to the number of datastreams")
  }
  # if(length(quantile) == 1){
  #   quantile = rep(quantile, nStreams)
  # } else if(length(quantile) != nStreams){
  #   stop("'quantile' must be of length 1 or equal to the number of datastreams")
  # }
  if(length(diff_order) == 1){
    diff_order = rep(diff_order, nStreams)
  } else if(length(diff_order) != nStreams){
    stop("'diff_order' must be of length 1 or equal to the number of datastreams")
  }
  
  listseed = vector("list", length = nStreams)
  names(listseed) = varnames
  Z = S = Z_predict = xseq = betastart = basis = listseed
  
  for(i in seq_len(nStreams)){
    thisname = varnames[i]
    
    msg <- paste0("Data stream: ", thisname)
    message(msg)
    
    modmat <- make_matrices_dens(x = data[[thisname]], 
                                type = type[i], 
                                k = k[i], 
                                knots = knots[[thisname]],
                                # quantile = quantile[i],
                                degree = degree[i], 
                                diff_order = diff_order[i])
    Z[[thisname]] = modmat$Z
    S[[thisname]] = modmat$S
    Z_predict[[thisname]] = modmat$Z_predict
    xseq[[thisname]] = modmat$xseq
    basis[[thisname]] = modmat$basis
    betastart[[thisname]] = make_splinecoef(modmat, type = type[i], par = par[[thisname]])
  }
  
  out <- list(
    Z = Z,
    S = S,
    coef = betastart,
    Z_predict = Z_predict,
    xseq = xseq,
    basis = basis
  )
  
  class(out) <- "SmoothDens"
  return(out)
}


#' Compute the design matrix for a trigonometric basis expansion
#'
#' Given a periodically varying variable such as time of day or day of year and the associated cycle length, this function performs a basis expansion to efficiently calculate a linear predictor of the form
#' \deqn{ 
#'  \eta^{(t)} = \beta_0 + \sum_{k=1}^K \bigl( \beta_{1k} \sin(\frac{2 \pi k t}{L}) + \beta_{2k} \cos(\frac{2 \pi k t}{L}) \bigr). 
#'  }
#'  This is relevant for modeling e.g. diurnal variation and the flexibility can be increased by adding smaller frequencies (i.e. increasing \eqn{K}).
#'  
#' @param tod equidistant sequence of a cyclic variable
#' 
#' For time of day and e.g. half-hourly data, this could be 1, ..., L and L = 48, or 0.5, 1, 1.5, ..., 24 and L = 24.
#' @param L length of one cycle on the scale of the time variable. For time of day, this would be 24.
#' @param degree degree K of the trigonometric link above. Increasing K increases the flexibility.
#'
#' @return design matrix (without intercept column), ordered as sin1, cos1, sin2, cos2, ...
#'
#' @examples
#' # no examples
trigBasisExp = function(tod, L = 24, degree = 1){
  n = length(tod)
  Z = matrix(nrow = n, ncol = 2*degree)
  inner = 2*pi*tod/L
  for(k in seq_len(degree)){
    Z[,2*(k-1)+1:2] = cbind(sin(inner*k), cos(inner*k))
  }
  colnames(Z) = paste0(c("sin_", "cos_"), rep(1:degree, each = 2))
  Z
}
