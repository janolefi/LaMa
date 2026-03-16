# [Forward algorithm](https://www.taylorfrancis.com/books/mono/10.1201/b20790/hidden-markov-models-time-series-walter-zucchini-iain-macdonald-roland-langrock) with for periodically varying transition probability matrices

Calculates the log-likelihood of a sequence of observations under a
hidden Markov model with periodically varying transition probabilities
using the **forward algorithm**.

## Usage

``` r
forward_p(
  delta,
  Gamma,
  allprobs,
  tod,
  trackID = NULL,
  ad = NULL,
  report = TRUE,
  logspace = FALSE
)
```

## Arguments

- delta:

  initial or stationary distribution of length N, or matrix of dimension
  c(k,N) for k independent tracks, if `trackID` is provided

- Gamma:

  array of transition probability matrices of dimension c(N,N,L).

  Here we use the definition \\\Pr(S_t=j \mid S\_{t-1}=i) =
  \gamma\_{ij}^{(t)}\\ such that the transition probabilities between
  time point \\t-1\\ and \\t\\ are an element of \\\Gamma^{(t)}\\.

- allprobs:

  matrix of state-dependent probabilities/ density values of dimension
  c(n, N)

- tod:

  (Integer valued) variable for cycle indexing in 1, ..., L, mapping the
  data index to a generalised time of day (length n)

  For half-hourly data L = 48. It could, however, also be day of year
  for daily data and L = 365.

- trackID:

  optional vector of length n containing IDs

  If provided, the total log-likelihood will be the sum of each track's
  likelihood contribution. Instead of a single vector `delta`
  corresponding to the initial distribution, a `delta` matrix of initial
  distributions of dimension c(k,N), can be provided, such that each
  track starts with it's own initial distribution.

- ad:

  optional logical, indicating whether automatic differentiation with
  `RTMB` should be used. By default, the function determines this
  itself.

- report:

  logical, indicating whether `delta`, `Gamma`, `allprobs`, and
  potentially `trackID` should be reported from the fitted model.
  Defaults to `TRUE`, but only works if `ad = TRUE`, as it uses the
  `RTMB` package.

  **Caution:** When there are multiple tracks, for compatibility with
  downstream functions like
  [`viterbi_p`](https://janolefi.github.io/LaMa/reference/viterbi_p.md),
  [`stateprobs_p`](https://janolefi.github.io/LaMa/reference/stateprobs_p.md)
  or
  [`pseudo_res`](https://janolefi.github.io/LaMa/reference/pseudo_res.md),
  `forward_p` should only be called **once** with a `trackID` argument.

- logspace:

  logical, indicating whether the probabilities/ densities in the
  `allprobs` matrix are on log-scale. If so, internal computations are
  also done on log-scale which is numerically more robust when the
  entries are very small.

## Value

HMM log-likelihood for given data and parameters

## Details

When the transition probability matrix only varies periodically (e.g. as
a function of time of day), there are only \\L\\ unique matrices if
\\L\\ is the period length (e.g. \\L=24\\ for hourly data and
time-of-day variation). Thus, it is much more efficient to only
calculate these \\L\\ matrices and index them by a time variable (e.g.
time of day or day of year) instead of calculating such a matrix for
each index in the data set (which would be redundant). This function
allows for that by only expecting a transition probability matrix for
each time point in a period and an integer valued (\\1, \dots, L\\) time
variable that maps the data index to the according time.

## See also

Other forward algorithms:
[`forward()`](https://janolefi.github.io/LaMa/reference/forward.md),
[`forward_g()`](https://janolefi.github.io/LaMa/reference/forward_g.md),
[`forward_hsmm()`](https://janolefi.github.io/LaMa/reference/forward_hsmm.md),
[`forward_ihsmm()`](https://janolefi.github.io/LaMa/reference/forward_ihsmm.md),
[`forward_phsmm()`](https://janolefi.github.io/LaMa/reference/forward_phsmm.md)

## Examples

``` r
## negative log likelihood function
nll = function(par, step, tod) {
 # parameter transformations for unconstrained optimisation
 beta = matrix(par[1:6], nrow = 2)
 Gamma = tpm_p(1:24, beta = beta) # multinomial logit link for each time point
 delta = stationary_p(Gamma, tod[1]) # stationary HMM
 mu = exp(par[7:8])
 sigma = exp(par[9:10])
 # calculate all state-dependent probabilities
 allprobs = matrix(1, length(step), 2)
 ind = which(!is.na(step))
 for(j in 1:2) allprobs[ind,j] = dgamma2(step[ind], mu[j], sigma[j])
 # simple forward algorithm to calculate log-likelihood
 -forward_p(delta, Gamma, allprobs, tod)
}

## fitting an HMM to the nessi data
par = c(-2,-2,            # initial tpm intercepts (logit-scale)
        rep(0, 4),        # initial tpm slopes
        log(c(0.3, 2.5)), # initial means for step length (log-transformed)
        log(c(0.2, 1.5))) # initial sds for step length (log-transformed)
mod = nlm(nll, par, step = trex$step[1:500], tod = trex$tod[1:500])
```
