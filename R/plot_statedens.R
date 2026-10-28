# State-dependent densities from the streams recorded by statedens(), for plot(mod, type = "statedens")

# streams of mod, after checking that allprobs is their product (nothing else entered it)
get_streams <- function(mod) {
  st <- mod$streams
  if(is.null(st)) stop("No data streams found: build allprobs with statedens() to plot or simulate the state-dependent densities.")
  P <- Reduce(`*`, lapply(st, function(s) dens_matrix(s$x, s$dens, s$pars, s$nStates)))
  if(!isTRUE(all.equal(P, mod$allprobs, check.attributes = FALSE))) {
    stop("allprobs is not the product of the data streams from statedens() (",
         paste(names(st), collapse = ", "), "), e.g. a mixture or another density also enters it.")
  }
  st
}

# state weights for each row: local state probabilities, decoded states, stationary distribution or given
state_weights <- function(mod, weights) {
  n <- nrow(mod$allprobs)
  N <- ncol(mod$allprobs)
  if(is.numeric(weights)) {
    if(length(weights) != N) stop("Numeric weights need to have length nStates.")
    return(matrix(weights, n, N, byrow = TRUE))
  }
  switch(weights,
         stateprobs = stateprobs(mod),
         viterbi = outer(viterbi(mod), 1:N, "==") * 1,
         stationary = {
           s <- summary(mod) # averaged over tracks, or over the cycle for periodic models
           d <- if(is.null(s$stationary)) s$occupancy else s$stationary
           if(is.null(d)) stop("No stationary distribution for this state process, use weights = \"stateprobs\".")
           matrix(d, n, N, byrow = TRUE)
         },
         stop("weights must be \"stateprobs\", \"viterbi\", \"stationary\" or a vector of length nStates."))
}

# weighted state-dependent densities of stream st on a grid: column j is the mean over the observed rows
# of p_tj f_j(x; theta_tj), so the columns add up to the fitted marginal density
stream_curves <- function(st, W, n = 200, xlim = NULL) {
  obs <- which(!is.na(st$x))
  W <- W[obs, , drop = FALSE]
  fam <- sub("^d", "", st$name)
  discrete <- fam %in% .discrete_fams
  circular <- fam %in% c("vm", "wrpcauchy")
  r <- range(st$x[obs])
  if(!is.null(xlim) && !circular) r <- if(discrete) c(ceiling(xlim[1]), floor(xlim[2])) else xlim # grid on the plotted range
  x <- if(discrete) r[1]:r[2] else if(circular) seq(-pi, pi, length.out = n) else seq(r[1], r[2], length.out = n)

  dens <- sapply(1:st$nStates, function(j) {
    pj <- state_pars(st$pars, j, obs)
    if(all(lengths(pj) == 1)) return(mean(W[, j]) * do.call(st$dens, c(list(x), pj)))
    # parameters per observation: average over the rows for each grid point
    sapply(x, function(g) mean(W[, j] * do.call(st$dens, c(list(rep(g, length(obs))), pj))))
  })
  colnames(dens) <- paste0("S", 1:st$nStates)
  list(x = x, density = dens, weights = colMeans(W), discrete = discrete, circular = circular, data = st$x[obs])
}

# histogram of each stream with the weighted state-dependent densities and their sum, in the style of the vignettes
plot_statedens <- function(mod, stream = NULL, weights = "stateprobs", col = NULL, lwd = 2, n = 200,
                           legend = "topright", ...) {
  st <- get_streams(mod)
  if(is.null(stream)) stream <- names(st)
  if(!all(stream %in% names(st))) stop("stream needs to be one of ", paste(names(st), collapse = ", "), ".")
  N <- ncol(mod$allprobs)
  if(is.null(col)) col <- LaMaColors(N)
  W <- state_weights(mod, weights)

  # one panel per stream, unless the user has set par(mfrow)
  if(length(stream) > 1 && all(par("mfrow") == 1)) {
    old <- par(mfrow = c(1, length(stream)))
    on.exit(par(old))
  }

  out <- lapply(setNames(stream, stream), function(s) {
    cv <- stream_curves(st[[s]], W, n, list(...)$xlim)
    total <- rowSums(cv$density)

    # histogram, arguments in ... replace the defaults
    # default breaks: one bar per value if discrete, equal bins on [-pi, pi] if circular
    breaks <- if(cv$discrete) seq(min(cv$x) - 0.5, max(cv$x) + 0.5) else if(cv$circular) seq(-pi, pi, length.out = 25) else 40
    args <- modifyList(list(breaks = breaks,
                            border = "white", main = "", xlab = s), list(...))
    h <- hist(cv$data, breaks = args$breaks, plot = FALSE)
    args$breaks <- NULL
    if(is.null(args$xlim)) args$xlim <- range(h$breaks, cv$x)
    if(is.null(args$ylim)) args$ylim <- c(0, max(h$density, total))
    do.call(plot, c(list(h, freq = FALSE), args))

    # densities as lines, probabilities as points connected by lines
    type <- if(cv$discrete) "b" else "l"
    for(j in 1:N) lines(cv$x, cv$density[, j], col = col[j], lwd = lwd, type = type, pch = 20)
    lines(cv$x, total, lwd = lwd, lty = 2, type = type, pch = 20)
    if(!isFALSE(legend)) {
      graphics::legend(legend, legend = c(paste("state", 1:N), "marginal"), col = c(col[1:N], 1),
                       lwd = lwd, lty = c(rep(1, N), 2), bty = "n")
    }
    cv
  })
  if(length(out) == 1) out[[1]] else out
}
