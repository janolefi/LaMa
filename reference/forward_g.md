# General [forward algorithm](https://www.taylorfrancis.com/books/mono/10.1201/b20790/hidden-markov-models-time-series-walter-zucchini-iain-macdonald-roland-langrock) with time-varying transition probability matrix

Calculates the log-likelihood of a sequence of observations under a
hidden Markov model with time-varying transition probabilities using the
**forward algorithm**.

## Usage

``` r
forward_g(
  delta,
  Gamma,
  allprobs,
  trackID = NULL,
  logspace = FALSE,
  ad = NULL,
  bw = NULL,
  report = TRUE
)
```

## Arguments

- delta:

  initial or stationary distribution of length `N`, or matrix of
  dimension `c(k,N)` for `k` independent tracks, if `trackID` is
  provided

- Gamma:

  array of transition probability matrices of dimension `c(N,N,n-1)`, as
  in a time series of length `n`, there are only `n-1` transitions.

  If an array of dimension `c(N,N,n)` for a single track is provided,
  the first slice will be ignored.

  If the elements of \\\Gamma^{(t)}\\ depend on covariate values at t or
  covariates \\t+1\\ is your choice in the calculation of the array,
  prior to using this function. When conducting the calculation by using
  [`tpm_g()`](https://janolefi.github.io/LaMa/reference/tpm_g.md), the
  choice comes down to including the covariate matrix `Z[-1,]` oder
  `Z[-n,]`.

  If `trackID` is provided, Gamma needs to be an array of dimension
  `c(N,N,n)`, matching the number of rows of allprobs. For each track,
  the transition matrix at the beginning will be ignored. If the
  parameters for Gamma are pooled across tracks or not, depends on your
  calculation of Gamma. If pooled, you can use `tpm_g(Z, beta)` to
  calculate the entire array of transition matrices when `Z` is of
  dimension `c(n,p)`.  

  This function can also be used to fit continuous-time HMMs, where each
  array entry is the Markov semigroup \\\Gamma(\Delta t) = \exp(Q \Delta
  t)\\ and \\Q\\ is the generator of the continuous-time Markov chain.

- allprobs:

  matrix of state-dependent probabilities/ density values of dimension
  `c(n, N)`

- trackID:

  optional vector of length `n` containing IDs that separate tracks.

  If provided, the total log-likelihood will be the sum of each track's
  likelihood contribution. In this case, `Gamma` must be an array of
  dimension `c(N,N,n)`, matching the number of rows of allprobs. For
  each track, the transition matrix at the beginning of the track will
  be ignored (as there is no transition between tracks). Furthermore,
  instead of a single vector `delta` corresponding to the initial
  distribution, a `delta` matrix of initial distributions, of dimension
  `c(k,N)`, can be provided, such that each track starts with it's own
  initial distribution.

- logspace:

  logical, indicating whether the probabilities/ densities in the
  `allprobs` matrix are on log-scale. If so, internal computations are
  also done on log-scale which is numerically more robust when the
  entries are very small. Note that this is only supported when used in
  AD mode with `RTMB`.

- ad:

  optional logical, indicating whether automatic differentiation with
  `RTMB` should be used. By default, the function determines this
  itself.

- bw:

  optional integer, indicating the bandwidth for a banded approximation
  of the forward algorithm. This is for expert users only, if sparsity
  in the Hessian matrix w.r.t. observations is required.

- report:

  logical, indicating whether `delta`, `Gamma`, `allprobs`, and
  potentially `trackID` should be reported from the fitted model.
  Defaults to `TRUE`, but only works if `ad = TRUE`, as it uses the
  `RTMB` package.

  When there are multiple tracks, for compatibility with downstream
  functions like
  [`viterbi_g`](https://janolefi.github.io/LaMa/reference/viterbi_g.md),
  [`stateprobs_g`](https://janolefi.github.io/LaMa/reference/stateprobs_g.md)
  or
  [`pseudo_res`](https://janolefi.github.io/LaMa/reference/pseudo_res.md),
  `forward_g` should only be called **once** with a `trackID` argument.

## Value

log-likelihood for given data and parameters

## See also

Other forward algorithms:
[`forward()`](https://janolefi.github.io/LaMa/reference/forward.md),
[`forward2()`](https://janolefi.github.io/LaMa/reference/forward2.md),
[`forward_g2()`](https://janolefi.github.io/LaMa/reference/forward_g2.md),
[`forward_hsmm()`](https://janolefi.github.io/LaMa/reference/forward_hsmm.md),
[`forward_ihsmm()`](https://janolefi.github.io/LaMa/reference/forward_ihsmm.md),
[`forward_p()`](https://janolefi.github.io/LaMa/reference/forward_p.md),
[`forward_phsmm()`](https://janolefi.github.io/LaMa/reference/forward_phsmm.md)

## Examples

``` r
## Simple usage
Gamma = array(c(0.9, 0.2, 0.1, 0.8), dim = c(2,2,10))
delta = c(0.5, 0.5)
allprobs = matrix(0.5, 10, 2)
forward_g(delta, Gamma, allprobs)
#> [1] -6.931472
# \donttest{
## Full model fitting example
## negative log likelihood function
nll = function(par, step, Z) {
 # parameter transformations for unconstrained optimisation
 beta = matrix(par[1:6], nrow = 2)
 Gamma = tpm_g(Z, beta) # multinomial logit link for each time point
 delta = stationary(Gamma[,,1]) # stationary HMM
 mu = exp(par[7:8])
 sigma = exp(par[9:10])
 # calculate all state-dependent probabilities
 allprobs = matrix(1, length(step), 2)
 ind = which(!is.na(step))
 for(j in 1:2) allprobs[ind,j] = dgamma2(step[ind], mu[j], sigma[j])
 # simple forward algorithm to calculate log-likelihood
 -forward_g(delta, Gamma, allprobs)
}

## fitting an HMM to the trex data
par = c(-1.5,-1.5,        # initial tpm intercepts (logit-scale)
        rep(0, 4),        # initial tpm slopes
        log(c(0.3, 2.5)), # initial means for step length (log-transformed)
        log(c(0.2, 1.5))) # initial sds for step length (log-transformed)
mod = nlm(nll, par, step = trex$step[1:500], Z = cosinor(trex$tod[1:500]))
#> Warning: NA/NaN replaced by maximum positive value
# }
```
