# Optimise an `RTMB` object from several starting values

Latent Markov models often have several local optima of the likelihood.
`multistart()` optimises `obj` from several starting values: the first
is `obj$par` (the initial parameters), the others are random
perturbations of it on the working scale. Afterwards,
[`lama`](https://janolefi.github.io/LaMa/reference/lama.md)`(obj)`
builds the model from the best fit, as `RTMB` keeps the parameters with
the lowest negative log-likelihood it evaluated.

## Usage

``` r
multistart(
  obj,
  nstarts = 20,
  sd = 1,
  starts = NULL,
  scale = c("fixed", "unit_info"),
  stop = c("bayes", "hits", "none"),
  hits = 3,
  sd_max = 2,
  tol = 1e-06,
  seed = NULL,
  optimizer = stats::nlminb,
  silent = FALSE,
  ...
)
```

## Arguments

- obj:

  `RTMB` object from `MakeADFun()`

- nstarts:

  maximum number of starts, including `obj$par`

- sd:

  standard deviation of the perturbations on the working scale: a
  number, a vector of the length of `obj$par`, or a named list with one
  value per parameter (names as in the parameter list)

- starts:

  optional matrix of starting values, one per row, with columns as
  `obj$par`; then `nstarts` and `sd` are ignored

- scale:

  `"fixed"` for perturbations of `obj$par` with `sd`, or `"unit_info"`
  for perturbations of the first fit scaled by the unit information, see
  ‘Details’

- stop:

  stopping rule: `"bayes"` (not used for given `starts`), `"hits"` or
  `"none"`, see ‘Details’

- hits:

  number of times the best optimum needs to be reached for
  `stop = "hits"`

- sd_max:

  largest standard deviation for `scale = "unit_info"`

- tol:

  relative tolerance for two runs to reach the same optimum

- seed:

  optional seed for the random number generator

- optimizer:

  optimisation function called as `optimizer(par, fn, gr, ...)`,
  returning `par` and `objective` (like
  [`nlminb`](https://rdrr.io/r/stats/nlminb.html)) or `value` (like
  [`optim`](https://rdrr.io/r/stats/optim.html))

- silent:

  logical; if `TRUE`, nothing is printed

- ...:

  further arguments passed to `optimizer`

## Value

output of `optimizer` for the best run, with two more elements: `runs`,
a data frame with the negative log-likelihood `nll`, `convergence` code
and maximum gradient component `mgc` of each run, and `starts`, the
matrix of starting values used.

## Details

**Stopping.** Two runs reach the same optimum if their negative
log-likelihoods differ by less than `tol` (relative), so permuted states
(label switching) count as the same optimum. Only converged runs
(convergence code 0 and maximum gradient component below 0.01) count for
stopping. In any case, at most `nstarts` starts are run.

- `stop = "bayes"` uses the Bayesian stopping rule of Boender and
  Rinnooy Kan (1987): after \\k\\ runs with \\w\\ distinct optima, the
  expected number of optima is \\w(k-1)/(k-w-2)\\, and `multistart()`
  stops when this is below \\w + 0.5\\, i.e. probably no optimum is
  missed. The rule assumes that the starts are independent draws from
  one distribution, so only the random starts count, not the first run
  from `obj$par`, and it is not used for given `starts`. It needs at
  least 8 converged random starts.

- `stop = "hits"` stops once the best optimum was reached `hits` times.

**Caution.** Stopping does not mean that the global optimum was found.
"No optimum missed" only refers to the distribution of the starts:
optima far away from it, or with a small basin of attraction, can be
missed. For example, with a very small `sd`, all starts converge to the
same optimum and the rule stops quickly, although other starts might
find a better one. So choose `sd` wide enough and still check the fit.

**Scatter.**

- With `scale = "fixed"`, the starts are `obj$par` plus normal noise
  with standard deviation `sd`. As most working parameters are on log or
  logit scale, `sd = 1` is a relative change, so one value usually fits
  all parameters. Parameters on identity scale (e.g. means of normal
  distributions) may need their own `sd`, given as a named list like
  `list(eta = 1, mu = 5)`.

- With `scale = "unit_info"`, the starts are drawn around the first fit
  from a normal distribution with covariance \\n H^{-1}\\, where \\H\\
  is the Hessian of the negative log-likelihood and \\n\\ the number of
  observations. This is the information of a single observation (unit
  information), which does not shrink with \\n\\ and takes correlations
  of the parameters into account. Standard deviations are capped at
  `sd_max`, for weakly identified parameters. If the Hessian is not
  positive definite (or `allprobs` is not reported, so \\n\\ is
  unknown), the fixed `sd` is used instead.

## References

Boender, C. G. E., & Rinnooy Kan, A. H. G. (1987). Bayesian stopping
rules for multistart global optimization methods. *Mathematical
Programming*, 37(1), 59-80.

## Examples

``` r
nll <- function(par) {
  getAll(par, dat)
  Gamma <- tpm(eta)
  mu <- exp(log_mu); REPORT(mu)
  allprobs <- cbind(dgamma2(step, mu[1], exp(log_sigma[1])),
                   dgamma2(step, mu[2], exp(log_sigma[2])))
  -forward(stationary(Gamma), Gamma, allprobs)
}
dat <- list(step = trex$step[1:500])
par <- list(eta = c(-2, -2), log_mu = log(c(0.3, 2.5)), log_sigma = log(c(0.2, 1.5)))
obj <- MakeADFun(nll, par, silent = TRUE)

opt <- multistart(obj, seed = 1)
#> start 1 - nll: 536.2127 
#> start 2 - nll: 536.2127 
#> start 3 - nll: 536.2127 
#> start 4 - nll: 656.563 
#> start 5 - nll: 536.2127 
#> start 6 - nll: 536.2127 
#> start 7 - nll: 536.2127 
#> start 8 - nll: 536.2127 
#> start 9 - nll: 536.2127 
#> start 10 - nll: 536.2127 
#> 10 starts, 9 converged, 1 distinct optima, best reached 9 times; stopped: Bayesian stopping rule 
opt$runs
#>         nll convergence          mgc
#> 1  536.2127           0 1.657455e-04
#> 2  536.2127           0 9.364391e-05
#> 3  536.2127           0 2.764750e-03
#> 4  656.5630           0 1.919732e-02
#> 5  536.2127           0 2.028025e-04
#> 6  536.2127           0 7.778351e-04
#> 7  536.2127           0 1.693041e-05
#> 8  536.2127           0 6.041584e-04
#> 9  536.2127           0 2.866955e-03
#> 10 536.2127           0 3.189349e-05
mod <- lama(obj) # best fit

# starts scaled by the first fit
opt <- multistart(obj, scale = "unit_info", seed = 1)
#> start 1 - nll: 536.2127 
#> start 2 - nll: 536.2127 
#> start 3 - nll: 536.2127 
#> start 4 - nll: 661.2706 
#> start 5 - nll: 536.2127 
#> start 6 - nll: 536.2127 
#> start 7 - nll: 536.2127 
#> start 8 - nll: 536.2127 
#> start 9 - nll: 536.2127 
#> start 10 - nll: 536.2127 
#> start 11 - nll: 662.1921 
#> start 12 - nll: 536.2127 
#> start 13 - nll: 536.2127 
#> start 14 - nll: 536.2127 
#> start 15 - nll: 536.2127 
#> start 16 - nll: 536.2127 
#> start 17 - nll: 536.2127 
#> start 18 - nll: 536.2127 
#> start 19 - nll: 536.2127 
#> start 20 - nll: 536.2127 
#> 20 starts, 20 converged, 3 distinct optima, best reached 18 times; stopped: nstarts reached 
```
