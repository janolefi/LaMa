
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

#' Build the design and the penalty matrix for models involving penalised splines based on a formula and a data set
#'
#' @param formula right side of a formula as used in \code{mgcv}
#' @param data data frame containing the variables in the formula
#' @param knots optional list containing user specified knot values to be used for basis construction
#' 
#' For most bases the user simply supplies the \code{knots} to be used, which must match up with the \code{k} value supplied (note that the number of knots is not always just \code{k}).
#' See \code{mgcv} documentation for more details.
#'
#' @return a list containing the design matrix \code{Z}, a (potentially nested) list of penalty matrices \code{S}, the \code{formula}, the \code{data}, the \code{knots}, and the original \code{mod} object returned by \code{mgcv}.
#' Note that for tensorproduct smooths, the corresponding list entry is itself a list, containing the d marginal penalty matrices if d is the dimension of the tensor product.
#' @export
#' 
#' @importFrom mgcv gam s
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
make_matrices_old = function(formula, 
                         data, 
                         knots = NULL
                         ){
  
  ## Potenially expand cosinor terms
  formula = expand_cosinor(formula)
  
  ## setting up the model using mgcv
  gam_setup = gam(formula = update(formula, dummy ~ .),
                  data = cbind(dummy = 1, data), 
                  knots = knots,
                  fit = FALSE)
  
  ## assiging design matrix
  term_names = gam_setup$term.names
  Z = gam_setup$X
  colnames(Z) = term_names

  ## dealing with the penalty matrices
  term_labels = sapply(gam_setup$smooth, function(x) x$label)
  
  # second option: tensorproduct -> save blown-up marginal penalty matrices (with constraints baked in)
  S = list()
  counter = 1
  for(i in seq_along(gam_setup$smooth)){
    sm = gam_setup$smooth[[i]]
    if(is.null(sm$margin)){
      S[[i]] = gam_setup$S[[counter]]
      counter = counter + 1
    } else{
      nPenMat = length(sm$margin)
      S_sublist = gam_setup$S[counter:(counter + nPenMat - 1)]
      margin_names = sapply(sm$margin, function(y) y$term)
      names(S_sublist) = margin_names
      S[[i]] = S_sublist
      counter = counter + nPenMat
    }
  }
  names(S) = term_labels
  
  pardim <- list(fixed_eff = gam_setup$nsdf)
  
  pardim_smooth = sapply(S, function(x){
    if(is.matrix(x)){
      return(nrow(x))
    } else{
      return(nrow(x[[1]]))
    }
  })
  
  pardim = c(pardim, pardim_smooth)
  
  out = list(Z = Z, 
             S = S, 
             # S2 = S2,
             pardim = pardim,
             formula = gam_setup$formula, 
             data = data, 
             knots = knots,
             gam = gam_setup)
  
  class(out) <- "LaMa_matrices"
  return(out)
}

make_matrices_flat <- function(formula, data, knots = NULL) {
  
  get_name <- function(fml, idx) {
    if (length(fml) == 3 && !deparse(fml[[2]]) %in% c("", ".")) {
      deparse(fml[[2]])
    } else {
      paste0("par", idx)
    }
  }
  
  process_single <- function(fml, name, knots_sub) {
    fml <- expand_cosinor(fml)
    
    # prepare gam model setup
    gam_setup <- gam(update(fml, dummy ~ .), data = cbind(dummy = 1, data), 
                   knots = knots_sub, fit = FALSE)
    
    # also prepare prediction gam
    gam_setup0 <- gam(update(fml, dummy ~ .), data = cbind(dummy = 1, data), 
                   knots = knots_sub, control = list(maxit = 1))
    
    Z <- gam_setup$X
    colnames(Z) <- gam_setup$term.names
    
    S <- list()
    coef <- list()
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
        counter <- counter + n_pen
      }
    }
    
    coef_fixed <- rep(0, gam_setup$nsdf)
    names(coef_fixed) <- colnames(Z)[seq_len(gam_setup$nsdf)]
    # coef <- c(setNames(list(coef_fixed), paste0(name, ".fixed_eff")), coef)
    
    list(
      Z = Z, S = S, pardim = pardim,
      gam = gam_setup, gam0 = gam_setup0,
      knots = knots_sub, coef = coef
    )
  }
  
  if (!inherits(formula, "list")) {
    name <- get_name(formula, 1)
    if(name == "par1"){
      name <- NULL # default name for single formula
    }
    res <- process_single(formula, name, knots)
    out <- list(
      Z = res$Z, S = res$S, 
      pardim = res$pardim,
      coef = res$coef,
      data = data, 
      gam = res$gam, gam0 = res$gam0,
      knots = knots
    )
    class(out) <- "LaMa_matrices"
    return(out)
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
  
  # Case: list of formulas
  Z_list <- list()
  S_flat <- list()
  pardim_list <- list()
  formula_list <- list()
  knots_list <- list()
  gam_list <- list()
  gam0_list <- list()
  coef_list <- list()
  
  # formula_names <- names(formula)
  
  for (i in seq_along(formula)) {
    fml <- formula[[i]]
    name <- form_names[i]
    res <- process_single(fml, name, knots[[name]])
    
    Z_list[[name]] <- res$Z
    S_flat <- c(S_flat, res$S)
    pardim_list[[name]] <- res$pardim
    knots_list[[name]] <- res$knots
    gam_list[[name]] <- res$gam
    gam0_list[[name]] <- res$gam0
    coef_list <- c(coef_list, res$coef)
  }
  
  out <- list(
    Z = Z_list, S = S_flat, 
    pardim = pardim_list,
    coef = coef_list,
    data = data,
    gam = gam_list, gam0 = gam0_list,
    knots = knots_list
  )
  class(out) <- "LaMa_matrices"
  return(out)
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
#' \item{\code{data}}{the data frame used for the model(s)}
#' \item{\code{gam}}{unfitted \code{mgcv::gam} object used for construction of \code{Z} and \code{S} (or list of such objects if \code{formula} is a list)}
#' \item{\code{gam0}}{fitted \code{mgcv::gam} which is used internally to create prediction design matrices (or list of such objects if \code{formula} is a list)}
#' \item{\code{knots}}{knot list used in the basis construction (or named list over such lists if \code{formula} is a list)}
#'
#' @export
#' 
#' @importFrom mgcv gam s
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
  if(!is.list(formula)){
    # not a list -> check if formula is a single formula
    if(inherits(formula, "formula")){
      # check knots
      return(make_matrices_flat(formula, data, knots))
    } else{
      stop("'formula' must be a (nested) list of formulas or a single formula.")
    }
  } else if(all(sapply(formula, inherits, what = "formula"))){     # check if flat list
    # check knots
    if(!is.null(knots)){
      if(!is.list(knots)){
        stop("'knots' must be a list of knots for each formula in 'formula'.")
      } else if(!all(sapply(knots, is.list))){
        stop("'knots' must be a named list of knot lists, containing knots for the corresponding formula (top level) and variable (bottom level) in 'formula'.")
      }
    }
    return(make_matrices_flat(formula, data, knots))
  } else{
    # check if two level list
    nested_check <- all(sapply(formula, function(f){
      all(sapply(f, inherits, what = "formula"))
    }))
    
    if(nested_check){
      Z_list <- list()
      S_list <- list()
      pardim_list <- list()
      coef_list <- list()
      gam_list <- list()
      gam0_list <- list()
      knots_list <- list()
      
      names <- names(formula)
      if(is.null(names)){
        names <- paste0("stream", seq_along(formula))
      }
      
      for(i in seq_along(formula)){
        thisname <- names[i]
        res <- make_matrices_flat(formula[[i]], data, knots[[i]])
        
        # handle nameing of S and coef
        if(length(res$S) > 0){
          names(res$S) <- paste0(thisname, ".", names(res$S))
        }
        names(res$coef) <- paste0(thisname, ".", names(res$coef))
        
        Z_list[[thisname]] <- res$Z
        S_list <- c(S_list, res$S)
        pardim_list[[thisname]] <- res$pardim
        coef_list <- c(coef_list, res$coef)
        gam_list[[thisname]] <- res$gam
        gam0_list[[thisname]] <- res$gam0
        knots_list[[thisname]] <- res$knots
      }
      
      out <- list(
        Z = Z_list, 
        S = S_list, 
        pardim = pardim_list,
        coef = coef_list,
        data = data,
        gam = gam_list, 
        gam0 = gam0_list,
        knots = knots_list
      )
      class(out) <- "LaMa_matrices"
      return(out)
    } else{
      stop("'formula' must be a (nested) list of formulas or a single formula.")
    }
  }
}

# process_obs_formulas <- function(formulas, 
#                                  dists,
#                                  nStates){
# 
#   # get dists
#   dist_fns <- lapply(dists, function(d){
#     if(is.character("d")){
#       get(paste0("d", d))
#     } else if(is.function(d)){
#       d
#     }
#   })
#   # for each dist, find parameter names of dist
#   parnames <- lapply(dist_fns, function(d) formals(d)[-1])
#   
#   
#   # create nested named formula list with streami.par.i ~ 1 for all by default
#   # match with provided formulas
# }

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
#' @param ... for method consistency only
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
  # dots <- list(...)
  # if(!is.null(dots$newdata)){
  #   newdata <- dots$newdata
  # } else{
  #   newdata <- dots[[1]]
  # }
  
  pred_matrix(object, newdata = newdata, what = what)
}


#' Build the prediction design matrix based on new data and model_matrices object created by \code{\link{make_matrices}}
#'
#' @param model_matrices model_matrices object as returned from \code{\link{make_matrices}}
#' @param newdata data frame containing the variables in the formula and new data for which to evaluate the basis
#' @param what optional character string specifying which formula to use for prediction, if \code{object} contains multiple formulas. If \code{NULL}, the first formula is used.
#' @param exclude optional vector of terms to set to zero in the predicted design matrix. Useful for predicting main effects only when e.g. \code{sd(..., bs = "re")} terms are present. See \code{mgcv::predict.gam} for more details.
#' @return prediction design matrix for \code{newdata} with the same basis as used for \code{model_matrices}
#' @export
#' 
#' @importFrom mgcv gam
#' @importFrom mgcv s
#' @importFrom mgcv predict.gam
#'
#' @examples
#' # single formula
#' modmat = make_matrices(~ s(x), data.frame(x = 1:10))
#' Z_p = pred_matrix(modmat, data.frame(x = 1:10 - 0.5))
#' # with multiple formulas
#' modmat = make_matrices(list(mu ~ s(x), sigma ~ s(x, bs = "ps")), data = data.frame(x = 1:10))
#' Z_p = pred_matrix(modmat, data.frame(x = 1:10 - 0.5), what = "mu")
#' # nested formula list
#' form = list(stream1 = list(mu ~ s(x), sigma ~ s(x, bs = "ps")))
#' modmat = make_matrices(form, data = data.frame(x = 1:10))
#' Z_p = pred_matrix(modmat, data.frame(x = 1:10 - 0.5), what = c("stream1", "mu"))
pred_matrix = function(model_matrices, 
                       newdata,
                       what = NULL,
                       exclude = NULL) {
  
  if(is.null(model_matrices$gam0)){
    gam_setup0 = mgcv::gam(model_matrices$formula,
                           data = cbind(dummy = 1, model_matrices$data),
                           knots = model_matrices$knots)
  } else {
    if(inherits(model_matrices$gam0, "gam")){
      if(!is.null(what)){
        message("'model_matrics' only contains a single gam object, 'what' is not used.")
      }
      gam_setup0 <- model_matrices$gam0
    } else {
      if(!inherits(model_matrices$gam0[[1]], "gam")){
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
      } else {
        if(is.null(what)){
          stop("'what' must be specified and be one of the names/ response variables from your original formulas.")
          # what <- names(model_matrices$Z)[1]
        }
        if(!what %in% names(model_matrices$gam0)){
          stop("'what' must be one of the names of your original formulas.")
        }
      }
      
      gam_setup0 <- model_matrices$gam0[[what]]
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
