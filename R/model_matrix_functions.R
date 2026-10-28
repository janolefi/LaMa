
###### Regression setting ------------------------------------------------------

#' Trigonometric basis expansion
#'
#' Builds a design matrix of \code{sin}/\code{cos} pairs for use in models
#' with periodic predictors. Can be used directly or inside formulas passed
#' to \code{modmat} (where expansion is handled automatically).
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
#' modmat <- modmat(form, data = data)
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
  # drop NULL entries, e.g. knots when none were supplied
  per <- function(what) {
    v <- lapply(res, `[[`, what)
    keep <- !vapply(v, is.null, TRUE)
    if(any(keep)) v[keep] else list()
  }
  # c() qualifies the names by the component name, wanted only with 'prefix'
  flat <- function(what) do.call(c, if(prefix) per(what) else unname(per(what)))

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
    # right side as a string, reported by tpm() for summary()
    attr(Z, "formula") <- paste("~", gsub("\\s+", " ", paste(deparse(fml_given[[length(fml_given)]]), collapse = " ")))
    gam_setup$X <- NULL # returned as Z, and a second copy is the largest thing here
    
    # mgcv's initial guesses, one per penalty matrix and in the same order as
    # gam_setup$S, so the loop below can name them as it walks the smooths
    sp_init <- initial.sp(Z, gam_setup$S, gam_setup$off)
    # started well above mgcv's guess, as the inner fits are better behaved when
    # constrained initially; capped for terms with a large covariate scale
    sp_init <- pmin(100 * sp_init, 1e6)

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
      if(length(sm$S) == 0){
        stop("'", sm$label, "' is unpenalised (e.g. fx = TRUE), which is not supported.")
      }

      if (length(sm$S) == 1) {
        # Single-penalty smooth
        S[[label]] <- gam_setup$S[[counter]]
        pardim[[sm$label]] <- nrow(S[[label]])
        coef[[label]] <- rep(0, nrow(S[[label]]))
        sp0[label] <- sp_init[counter]
        counter <- counter + 1
      } else {
        # several penalties: one per margin for a tensor product, but also e.g. bs = "fs"
        n_pen <- length(sm$S)
        margin_names <- if(length(sm$margin) == n_pen) {
          vapply(sm$margin, function(m) paste(m$term, collapse = ","), "") # "x,y" for a 2-d margin
        } else seq_len(n_pen)
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
#' @details \code{make_matrices()} is the former name of \code{modmat()} and still works.
#'
#' @return a list of class \code{LaMa_matrices} containing:
#' \item{\code{Z}}{design matrix (or list of such matrices if \code{formula} is a list))}
#' \item{\code{S}}{list of penalty matrices (with names based on the response terms of the formulas as well as the smooth terms and covariates). For tensorproduct smooths, corresponding entries are themselves lists, containing the \eqn{d} marginal penalty matrices if \eqn{d} is the dimension of the tensor product)}
#' \item{\code{pardim}}{list of parameter dimensions (fixed and penalised separately) for each formula, for ease of setting up initial parameters}
#' \item{\code{coef}}{list of coefficient vectors filled with zeros of the correct length for each formula, for ease of setting up initial parameters}
#' \item{\code{sp0}}{named vector of \strong{initial penalty strengths}, with one entry per penalty matrix and names matching \code{S} (a smooth with several penalties, e.g. a tensor product, contributes one entry each, named \code{<smooth>.<margin>}).
#'
#' These are \code{\link[mgcv]{initial.sp}}'s guesses times 100, capped at \code{1e6}: the inner optimisation is better behaved when it starts constrained, and the outer iteration relaxes the penalty from there.
#'
#' It can be passed straight to \code{\link{penalty}} or \code{\link{penalty2}} when each smooth appears once. 
#' When a smooth is replicated, e.g. one spline per state, use \code{\link{rep_sp}}, as a single \code{rep()} over \code{sp0} is only correct for some models.}
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
#' modmat = modmat(~ s(x), data)
#' # univariate P-spline
#' modmat = modmat(~ s(x, bs = "ps"), data)
#' # adding random intercept
#' modmat = modmat(~ s(g, bs = "re") + s(x, bs = "ps"), data)
#' # tensorproduct of x and y
#' modmat = modmat(~ s(x) + s(y) + ti(x,y), data)
#' # multiple formulas at once
#' modmat = modmat(list(mu ~ s(x) + y, sigma ~ s(g, bs = "re")), data = data)
modmat <- function(formula, data, knots = NULL){
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

# former name, still works
#' @rdname modmat
#' @export
make_matrices <- modmat

#' Repeat initial penalty strengths for replicated smooths
#'
#' @description
#' Expands the \code{sp0} element of a \code{\link{modmat}} object into a
#' \code{sp} vector for \code{\link{penalty}} or \code{\link{penalty2}}, when each smooth is
#' used more than once, for example one spline per state or per off-diagonal element of the
#' transition probability matrix.
#'
#' @details
#' \code{sp} runs over replicates \strong{within} each smooth, so for a tensor product the
#' margins vary fastest. A single \code{rep()} over \code{sp0} therefore only gives the right
#' answer for some models, silently producing the right length but the wrong values for others.
#' This function repeats each smooth's penalty strengths separately, which is always correct.
#'
#' @param model_matrices object of class \code{LaMa_matrices} as returned by \code{\link{modmat}}
#' @param N number of replicates. Either a single number used for every smooth, or one number per
#' smooth, i.e. per element of \code{model_matrices$S}
#'
#' @return named numeric vector of initial penalty strengths, of the length and in the order
#' expected for \code{sp}
#'
#' @seealso \code{\link{modmat}}, which computes \code{sp0}
#' @export
#'
#' @examples
#' modmat = modmat(~ s(x) + ti(x, y), data.frame(x = runif(100), y = runif(100)))
#' modmat$sp0         # one value per penalty
#' rep_sp(modmat, 2)  # ready to use as 'sp' for two replicates of each smooth
rep_sp <- function(model_matrices, N) {
  if(!inherits(model_matrices, "LaMa_matrices")){
    stop("'model_matrices' needs to be an object returned by 'modmat()'")
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
#' @param x object of class \code{LaMa_matrices} as returned by \code{\link{modmat}}
#' @param ... ignored
#'
#' @return \code{x}, invisibly. Called for its side effect of printing.
#'
#' @seealso \code{\link{modmat}}
#' @export
print.LaMa_matrices <- function(x, ...) {
  # one line per formula, named as in the rest of the object, tildes lined up
  show <- function(f, nms, pad, width) {
    rhs <- vapply(f, function(g) paste(deparse(g[[length(g)]]), collapse = " "), "")
    lhs <- if(is.null(nms)) "" else paste0(formatC(nms, width = width, flag = "-"), " ")
    cat(paste0(pad, lhs, "~ ", gsub("\\s+", " ", rhs), "\n"), sep = "")
  }
  f <- x$formula
  cat("Model matrices for\n\n")
  if(inherits(f, "formula")){
    show(list(f), get_name(f), "  ", 0)
  } else if(inherits(f[[1]], "formula")){
    show(f, names(f), "  ", max(nchar(names(f))))
  } else { # one entry per data stream
    width <- max(0, nchar(unlist(lapply(f, names))))
    for(s in names(f)){
      cat("  ", s, "\n", sep = "")
      show(f[[s]], names(f[[s]]), "    ", width)
    }
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





#' Build the prediction design matrix based on new data and model_matrices object created by \code{\link{modmat}}
#'
#' @param object model matrices object as returned from \code{\link{modmat}}
#' @param newdata data frame containing the variables in the formula and new data for which to evaluate the basis
#' @param what optional character string specifying which formula to use for prediction if \code{object} contains multiple formulas.
#' @param exclude optional vector of terms to set to zero in the predicted design matrix. Useful for predicting main effects only when e.g. \code{s(..., bs = "re")} terms are present. See \code{mgcv::predict.gam} for more details.
#' @param ... for method consistency only
#'
#' @seealso \code{\link{modmat}} for creating objects of class \code{LaMa_matrices} which can be used for prediction by this function.
#'
#' @return prediction design matrix for \code{newdata} with the same basis as used for \code{object}
#' @export
#' 
#'
#' @examples
#' # single formula
#' modmat = modmat(~ s(x), data.frame(x = 1:10))
#' Z_p = predict(modmat, data.frame(x = 1:10 - 0.5))
#' # with multiple formulas
#' modmat = modmat(list(mu ~ s(x), sigma ~ s(x, bs = "ps")), data = data.frame(x = 1:10))
#' Z_p = predict(modmat, data.frame(x = 1:10 - 0.5), what = "mu")
#' # nested formula list
#' form = list(stream1 = list(mu ~ s(x), sigma ~ s(x, bs = "ps")))
#' modmat = modmat(form, data = data.frame(x = 1:10))
#' Z_p = predict(modmat, data.frame(x = 1:10 - 0.5), what = c("stream1", "mu"))
#' @importFrom mgcv gam
#' @importFrom mgcv s
#' @importFrom mgcv predict.gam
predict.LaMa_matrices <- function(object, newdata, what = NULL, exclude = NULL, ...){
  
  if(is.null(object$gam0)){
    stop("'object' contains no 'gam0' object; rebuild it with modmat()")
  } else {
    if(inherits(object$gam0, "gam")){
      if(!is.null(what)){
        message("'object' only contains a single gam object, 'what' is not used.")
      }
      gam_setup0 <- object$gam0
    } else {
      if(!inherits(object$gam0[[1]], "gam")){ # nested list of formulas
        if(is.null(what)){
          stop("'what' must be specified and contain a top-level name of the formula list and a name/ response variable from your original formulas.")
        }
        if(length(what) != 2 | !is.character(what)){
          stop("'what' must be a character vector of length 2, containing the top-level name of the formula list and a name/ response variable from your original formulas.")
        }
        if(!what[1] %in% names(object$gam0)){
          stop("'what[1]' must be a top-level name of your formula list.")
        }
        if(!what[2] %in% names(object$gam0[[what[1]]])){
          stop("'what[2]' must be one of the names/ response variables of your original formulas.")
        }
        gam_setup0 <- object$gam0[[what[1]]][[what[2]]]
      } else { # flat list of formulas
        if(is.null(what)){
          stop("'what' must be specified and be one of the names/ response variables from your original formulas.")
        }
        if(!what %in% names(object$gam0)){
          stop("'what' must be one of the names of your original formulas.")
        }
        gam_setup0 <- object$gam0[[what]]
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

#' Build the prediction design matrix based on new data and model_matrices object created by \code{\link{modmat}}
#'
#' @description
#' \strong{Deprecated.} Use \code{\link[=predict.LaMa_matrices]{predict}} on the object returned by \code{\link{modmat}} instead, which does the same thing:
#' \code{predict(modmat, newdata)} in place of \code{pred_matrix(modmat, newdata)}.
#'
#' @param model_matrices model_matrices object as returned from \code{\link{modmat}}
#' @param newdata data frame containing the variables in the formula and new data for which to evaluate the basis
#' @param what optional character string specifying which formula to use for prediction if \code{model_matrices} contains multiple formulas.
#' @param exclude optional vector of terms to set to zero in the predicted design matrix. Useful for predicting main effects only when e.g. \code{sd(..., bs = "re")} terms are present. See \code{mgcv::predict.gam} for more details.
#' @return prediction design matrix for \code{newdata} with the same basis as used for \code{model_matrices}
#'
#' @seealso \code{\link[=predict.LaMa_matrices]{predict}}, which replaces this function
#' @export
#'
#' @examples
#' modmat = modmat(~ s(x), data.frame(x = 1:10))
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
                          "returned by 'modmat()' instead."))
  predict.LaMa_matrices(model_matrices, newdata = newdata, what = what, exclude = exclude)
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
#' @param diff_order order of differencing used for the P-Spline penalty matrix. Defaults to second-order differences.
#' @param pow power for polynomial knot spacing, only used for \code{type = "positive"}
#' @param npoints number of points used in the numerical integration for normalizing the B-spline basis functions
#'
#' @return list containing the design matrix \code{Z}, the penalty matrix \code{S}, the prediction design matrix \code{Z_predict}, the prediction grid \code{xseq}, and details for the basis expansion.
#'
#' The penalty matrix leaves out the last coefficient, which is fixed at zero for identifiability.
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
                              diff_order = 2, # order of the differences for the penalty matrix
                              pow = 0.5, # power for polynomial knot spacing for positive values
                              npoints = 1e4 # number of points for numerical integration
){
  type <- match.arg(type, c("real", "positive", "circular"))

  if(type == "circular"){

    x <- (x + pi) %% (2 * pi) - pi # wrap angles to [-pi, pi)

    # if knots not supplied, default to equidistant knots
    if(is.null(knots)){
      knots <- seq(-pi, pi, length = k + 1)
    } else if(length(knots) != k + 1){
      stop("Number of knots provided is wrong. For 'type = circular' it should be 'k + 1' with '-pi' and 'pi' as first and last entries.")
    }

  } else {

    if(type == "positive" && min(x, na.rm = TRUE) <= 0){
      stop("When type = 'positive', 'x' can only contain positive values")
    }

    # computing the number of knots
    nrknots <- k - degree + 1

    if(is.null(knots)){
      # slightly inflated data range
      rangex <- range(x, na.rm = TRUE)
      rangex <- rangex + c(-1, 1) * diff(rangex) / 20

      if(type == "real"){ # equidistant knots
        knots <- seq(rangex[1], rangex[2], length = nrknots)
      } else { # polynomial spacing, always from 0 to max(x)
        knots <- seq(0, rangex[2]^pow, length = nrknots)^(1/pow)
      }
    } else if(length(knots) != nrknots){
      stop("Number of knots provided is wrong. It should be 'k - degree + 1'.")
    }
  }

  # midpoint grid for numerical integration
  h <- diff(range(knots)) / npoints
  xseq <- min(knots) + (seq_len(npoints) - 0.5) * h

  # numerical integration for normalizing the B-spline basis functions
  B0 <- dens_basis(xseq, type, knots, degree) # unnormalised spline design matrix
  w <- (h * colSums(B0))^-1
  B0 <- t(t(B0) * w)

  # basis positions for initial values later (expected value of each basis function)
  basis_pos <- colSums(xseq * t(t(B0) / rowSums(t(B0))))

  # actual data design matrix, NA rows stay NA
  B <- matrix(NA, nrow = length(x), ncol = k)
  ind <- which(!is.na(x))
  B[ind,] <- t(t(dens_basis(x[ind], type, knots, degree)) * w)

  # difference matrix, cyclic for circular
  if(type == "circular"){
    L <- diff(rbind(diag(k), diag(k)[1:diff_order,]), differences = diff_order)
  } else{
    L <- diff(diag(k), differences = diff_order)
  }

  # constructing penalty matrix
  S <- crossprod(L[,-k]) # leaving out last column because parameter set to zero

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
# unnormalised B-spline basis at x, knots including the boundary knots
dens_basis = function(x, type, knots, degree){
  if(type == "circular"){
    cSplineDes(x, knots, ord = degree + 1)
  } else{
    nk <- length(knots)
    bSpline(x, knots = knots[-c(1, nk)], Boundary.knots = knots[c(1, nk)],
            degree = degree, intercept = TRUE)
  }
}

# helper function, not exported
make_splinecoef = function(model_matrices,
                           type = "real",
                           par){
  basis_pos = model_matrices$basis$basis_pos
  k = length(basis_pos)

  if(type == "circular") {
    kappa = if(is.null(par$kappa)) par$concentration else par$kappa # accept both names
    beta = sapply(basis_pos, dvm, mu = par$mean, kappa = kappa, log = TRUE)
  } else { # normal on the reals, gamma on the positive reals
    ddist = if(type == "real") dnorm else dgamma2
    beta = sapply(basis_pos, ddist, mean = par$mean, sd = par$sd, log = TRUE)
  }
  if(is.vector(beta)) beta = matrix(beta, nrow = 1)

  # rescaling to account for non-equidistant knot spacing
  if(type != "circular") beta = t(t(beta) - log(apply(model_matrices$Z_predict, 2, max)))

  beta = beta - beta[,k]
  return(beta[,-k, drop = FALSE]) # stays a matrix for one state
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
#' \code{smoothdens} then constructs the design and penalty matrices for standardised B-splines basis functions (integrating to one) for each data stream.
#' For types \code{"real"} and \code{"circular"} the knots are placed equidistant in the range of the data, for type \code{"positive"} the knots are placed using polynomial spacing.
#'
#' The last spline coefficient of each state is fixed at zero for identifiability, hence \code{coef} has \code{k - 1} columns.
#' \code{density(beta)} takes care of this: it maps the coefficient matrix to the \code{allprobs} matrix, with rows of missing observations equal to one.
#'
#' @param x data vector for one data stream, or named data frame or list of several data streams
#' @param par list of initial means and sds (\code{"real"}, \code{"positive"}) or means and concentrations (\code{"circular"}, called \code{kappa} or \code{concentration}), one value per state.
#' For several data streams, a named list of such lists.
#' @param type type of each data stream, either \code{"real"} for data on the reals, \code{"positive"} for data on the positive reals or \code{"circular"} for angular data.
#' @param k number of basis functions for each data stream
#' @param knots optional vector of knots (including the boundary knots) to be used for basis construction, for several data streams a named list of such vectors.
#' If not provided, the knots are placed equidistantly for \code{"real"} and \code{"circular"} and using polynomial spacing for \code{"positive"}.
#'
#' For \code{"real"} and \code{"positive"} \code{k - degree + 1} knots are needed, for \code{"circular"} \code{k + 1} knots are needed.
#' @param degree degree of the B-spline basis functions for each data stream, defaults to cubic B-splines
#' @param diff_order order of differencing used for the P-Spline penalty matrix for each data stream. Defaults to second-order differences.
#' @param pow power for polynomial knot spacing, only used for \code{type = "positive"}
#' @param npoints number of points used in the numerical integration for normalizing the B-spline basis functions
#'
#' \code{type}, \code{k}, \code{degree}, \code{diff_order}, \code{pow} and \code{npoints} can be of length 1 or one per data stream (named or in the order of \code{x}).
#'
#' @return For one data stream, an object of class \code{"smoothdens"}: a list containing the design matrix \code{Z}, the penalty matrix \code{S}, the initial coefficients \code{coef}, the prediction design matrix \code{Z_predict}, the prediction grid \code{xseq}, details for the basis expansion \code{basis}, and the function \code{density(beta)}.
#' For several data streams, a named list of such objects.
#' @seealso \code{\link[=predict.smoothdens]{predict}} for the estimated densities
#' @export
#'
#' @examples
#' ## one data stream, mixture of two normals
#' x = c(rnorm(100, 0, 1), rnorm(100, 5, 1))
#' dens = smoothdens(x, par = list(mean = c(0, 5), sd = c(1, 1)))
#'
#' Z = dens$Z # design matrix
#' S = dens$S # penalty matrix
#' beta = dens$coef # initial coefficients
#' allprobs = dens$density(beta) # use this in the likelihood function
#'
#' ## 3 data streams
#' data = data.frame(x1 = rnorm(100, mean = 0, sd = 1),
#'                   x2 = rgamma2(100, mean = 5, sd = 3),
#'                   x3 = rvm(100, mu = 0, kappa = 2))
#'
#' par = list(x1 = list(mean = 0, sd = 1),
#'            x2 = list(mean = 5, sd = 3),
#'            x3 = list(mean = 0, kappa = 2))
#'
#' dens = smoothdens(data, par, type = c("real", "positive", "circular"))
#'
#' # objects for x1
#' Z1 = dens$x1$Z
#' allprobs = dens$x1$density(dens$x1$coef) *
#'   dens$x2$density(dens$x2$coef) *
#'   dens$x3$density(dens$x3$coef)
smoothdens <- function(x,
                       par,
                       type = "real",
                       k = 25,
                       knots = NULL,
                       degree = 3,
                       diff_order = 2,
                       pow = 0.5,
                       npoints = 1e4
){
  # several data streams: one object per stream
  if(is.list(x)){
    varnames <- names(x)
    if(is.null(varnames) || !all(varnames %in% names(par))){
      stop("'par' must be a named nested list with names corresponding to the datastreams")
    }
    # value of an argument for stream i: length 1, named or in order
    pick <- function(arg, i){
      if(length(arg) == 1) return(arg)
      if(length(arg) != length(varnames)){
        stop(paste0("'", deparse(substitute(arg)), "' must be of length 1 or equal to the number of datastreams"))
      }
      if(is.null(names(arg))) arg[[i]] else arg[[varnames[i]]]
    }
    out <- lapply(seq_along(varnames), function(i){
      thisname <- varnames[i]
      smoothdens(x[[thisname]], par[[thisname]],
                 type = pick(type, i), k = pick(k, i), knots = knots[[thisname]],
                 degree = pick(degree, i), diff_order = pick(diff_order, i),
                 pow = pick(pow, i), npoints = pick(npoints, i))
    })
    names(out) <- varnames
    return(out)
  }

  if(is.null(par$mean)){
    stop("'par' must contain initial means, e.g. list(mean = c(0, 5), sd = c(1, 1))")
  }

  modmat <- make_matrices_dens(x = x,
                               type = type,
                               k = k,
                               knots = knots,
                               degree = degree,
                               diff_order = diff_order,
                               pow = pow,
                               npoints = npoints)
  modmat$coef <- make_splinecoef(modmat, type = modmat$basis$type, par = par)
  modmat$density <- dens_fun(modmat$Z)

  out <- modmat[c("Z", "S", "coef", "Z_predict", "xseq", "basis", "density")]
  class(out) <- "smoothdens"
  return(out)
}

# helper function, not exported
# weights from coefficients, last one fixed at zero (multinomial logit)
dens_weights = function(beta){
  alpha <- exp(cbind(beta, 0))
  alpha / rowSums(alpha)
}

# helper function, not exported
# function mapping beta to allprobs, only keeps Z in its environment
dens_fun = function(Z){
  ind <- which(!is.na(Z[,1])) # only for non-NA obs.
  function(beta){
    alpha <- dens_weights(beta)
    allprobs <- matrix(1, nrow(Z), nrow(alpha))
    allprobs[ind,] <- Z[ind,] %*% t(alpha)
    allprobs
  }
}

#' Estimated smooth densities
#'
#' Evaluates the smooth densities built with \code{\link{smoothdens}} for given spline coefficients, e.g. to plot the estimated state-dependent densities.
#'
#' @param object object of class \code{"smoothdens"} for one data stream
#' @param coef coefficient matrix, one row per state, without the last column fixed at zero
#' @param x optional points to evaluate the densities at, defaults to the prediction grid \code{object$xseq}
#' @param ... unused
#'
#' @return list with the points \code{x} and the matrix \code{density} of densities, one row per point and one column per state
#' @export
#'
#' @examples
#' x = c(rnorm(100, 0, 1), rnorm(100, 5, 1))
#' dens = smoothdens(x, par = list(mean = c(0, 5), sd = c(1, 1)))
#' # densities at the initial coefficients
#' pred = predict(dens, coef = dens$coef)
#' matplot(pred$x, pred$density, type = "l", lty = 1)
predict.smoothdens = function(object, coef, x = NULL, ...){
  if(is.null(x)){
    x <- object$xseq
    Z <- object$Z_predict
  } else{
    b <- object$basis
    xw <- if(b$type == "circular") (x + pi) %% (2 * pi) - pi else x # wrap angles to [-pi, pi)
    Z <- t(t(dens_basis(xw, b$type, b$knots, b$degree)) * b$w)
  }
  list(x = x, density = Z %*% t(dens_weights(coef)))
}

#' Build the design and penalty matrices for smooth density estimation
#'
#' @description
#' \code{smooth_dens_construct()} is deprecated, use \code{\link{smoothdens}} instead.
#' It returns the same objects, but ordered by component first and data stream second.
#'
#' @param data named data frame of 1 or multiple data streams
#' @inheritParams smoothdens
#'
#' @return a nested list containing the design matrices \code{Z}, the penalty matrices \code{S}, the initial coefficients \code{coef} the prediction design matrices \code{Z_predict}, the prediction grids \code{xseq}, and details for the basis expansion for each data stream.
#' @export
#' @keywords internal
#'
#' @examples
#' data = data.frame(x = rnorm(100))
#' # deprecated:
#' # SmoothDens = smooth_dens_construct(data, par = list(x = list(mean = 0, sd = 1)))
#' # use instead:
#' dens = smoothdens(data$x, par = list(mean = 0, sd = 1))
smooth_dens_construct <- function(data,
                                  par,
                                  type = "real",
                                  k = 25,
                                  knots = NULL,
                                  degree = 3,
                                  diff_order = 2
){
  .Deprecated("smoothdens", package = "LaMa",
              msg = paste("'smooth_dens_construct()' is deprecated; use 'smoothdens()' instead."))

  if(!is.data.frame(data)){
    stop("data must be a data frame")
  }
  dens <- smoothdens(as.list(data), par, type = type, k = k, knots = knots,
                     degree = degree, diff_order = diff_order)

  # old layout: component first, then data stream
  comps <- c("Z", "S", "coef", "Z_predict", "xseq", "basis")
  out <- lapply(stats::setNames(comps, comps), function(n) lapply(dens, `[[`, n))

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
