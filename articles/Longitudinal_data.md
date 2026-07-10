# 4 Longitudinal data

> Before diving into this vignette, we recommend reading the vignettes
> [**Introduction to
> LaMa**](https://janolefi.github.io/LaMa/articles/Intro_to_LaMa.html)
> and [**Automatic differentiation via
> RTMB**](https://janolefi.github.io/LaMa/articles/LaMa_and_RTMB.html).

``` r

library(LaMa)
#> Loading required package: RTMB
```

In real-data applications, one will often be faced with a data set
consisting of several measurement tracks that can reasonably be assumed
to be mutually independent. Examples of such a longitudinal structure
include GPS tracks of several individuals (or several tracks, e.g. days,
of one individual), or when analysing sports data, time series for
separate games. In such settings, the researcher must decide whether to
pool parameters across tracks or not. Here, we provide brief examples
for complete and partial pooling.

When the data consists of K independent tracks, the joint likelihood is
L(\theta) = \prod\_{k=1}^K L_k(\theta), where L_k(\theta) is the usual
HMM likelihood for the k-th track, so the log-likelihood is a sum over K
tracks. In `LaMa`, this is handled efficiently by passing a `trackID`
vector to
[`forward()`](https://janolefi.github.io/LaMa/reference/forward.md),
which specifies the track index of each observation. The forward
algorithm then resets to the initial distribution at the start of each
new track and sums up the individual log-likelihood contributions
automatically. Importantly, always using the `trackID` argument rather
than calling
[`forward()`](https://janolefi.github.io/LaMa/reference/forward.md)
separately for each track ensures that the automatic reporting mechanism
inside
[`forward()`](https://janolefi.github.io/LaMa/reference/forward.md) — on
which
[`viterbi()`](https://janolefi.github.io/LaMa/reference/viterbi.md),
[`stateprobs()`](https://janolefi.github.io/LaMa/reference/stateprobs.md),
and other downstream functions rely — works correctly.

## Complete pooling

With complete pooling, all parameters are shared across individuals.

### Generating data

We generate K separate tracks from the same model:

``` r

# parameters shared across individuals
mu = c(15, 60)    # state-dependent means
sigma = c(10, 40) # state-dependent standard deviations
Gamma = matrix(c(0.95, 0.05,
                 0.15, 0.85), nrow = 2, byrow = TRUE)
delta = stationary(Gamma)

set.seed(123)
K = 200 # number of individuals
N = 2   # number of states
n = 50  # observations per individual

s = x = rep(NA, n * K)
for(k in 1:K){
  sk = xk = rep(NA, n)
  sk[1] = sample(1:2, 1, prob = delta)
  xk[1] = rnorm(1, mu[sk[1]], sigma[sk[1]])
  for(t in 2:n){
    sk[t] = sample(1:2, 1, prob = Gamma[sk[t-1],])
    xk[t] = rnorm(1, mu[sk[t]], sigma[sk[t]])
  }
  s[(k-1)*n + 1:n] = sk
  x[(k-1)*n + 1:n] = xk
}

trackID = rep(1:K, each = n)
```

### Writing the negative log-likelihood function

The likelihood function is almost identical to the single-track case —
the only difference is passing `trackID` as an additional argument to
[`forward()`](https://janolefi.github.io/LaMa/reference/forward.md):

``` r

nll_pool = function(par) {
  getAll(par, dat)
  Gamma = tpm(eta)
  delta = stationary(Gamma)
  mu = exp(log_mu); REPORT(mu)
  sigma = exp(log_sigma); REPORT(sigma)
  allprobs = matrix(1, length(x), N)
  for(j in 1:N) allprobs[,j] = dnorm(x, mu[j], sigma[j])
  -forward(delta, Gamma, allprobs, trackID)
}
```

### Estimating the model

``` r

par = list(
  eta = rep(-2, N*(N-1)),    # initial t.p.m. parameters (logit scale)
  log_mu = log(c(15, 60)),   # initial state-dependent means (log scale)
  log_sigma = log(c(10, 40)) # initial state-dependent sds (log scale)
)

dat = list(
  x = x,
  trackID = trackID,
  N = N
)

obj_pool = MakeADFun(nll_pool, par, silent = TRUE)
system.time(
  opt_pool <- nlminb(obj_pool$par, obj_pool$fn, obj_pool$gr)
)
#>    user  system elapsed 
#>    0.17    0.00    0.17
mod_pool = report(obj_pool)
```

We can inspect the estimated parameters, which should be close to the
true values used in the simulation:

``` r

mod_pool$mu     # true: 15, 60
#> [1] 15.00258 59.10962
mod_pool$sigma  # true: 10, 40
#> [1]  9.913063 39.856396
mod_pool$Gamma  # true: c(0.95, 0.05, 0.15, 0.85)
#>           S1         S2
#> S1 0.9500705 0.04992952
#> S2 0.1450977 0.85490231
```

## Partial pooling

If some parameters are individual-specific while the rest are shared, we
speak of partial pooling. We demonstrate this for K = 5 individuals
whose transition probability matrices depend on an individual-level
covariate, while the state-dependent observation parameters remain
shared. We estimate the effect of this covariate on the transition
probabilities.

### Generating data

``` r

K = 5   # number of individuals
N = 2   # number of states
n = 200 # observations per individual

# state-dependent parameters shared across individuals
mu = c(15, 60)
sigma = c(10, 40)

# individual-specific tpms depending on a covariate (e.g. age)
set.seed(123)
z = rnorm(K)
beta = matrix(c(-2, -2, 1, -1), nrow = N*(N-1))

# compute K individual-specific tpms
Gamma = tpm(beta, matrix(z, ncol = 1))
# stationary() on an array returns a c(K, N) matrix of stationary distributions
Delta = stationary(Gamma)

# simulate all tracks
set.seed(123)
s = x = rep(NA, n * K)
for(k in 1:K){
  sk = xk = rep(NA, n)
  sk[1] = sample(1:2, 1, prob = Delta[k,])
  xk[1] = rnorm(1, mu[sk[1]], sigma[sk[1]])
  for(t in 2:n){
    sk[t] = sample(1:2, 1, prob = Gamma[sk[t-1],,k])
    xk[t] = rnorm(1, mu[sk[t]], sigma[sk[t]])
  }
  s[(k-1)*n + 1:n] = sk
  x[(k-1)*n + 1:n] = xk
}
```

### Writing the negative log-likelihood function

Each individual has its own constant transition probability matrix,
determined by its covariate value. We compute all K individual-specific
tpms at once using
[`tpm()`](https://janolefi.github.io/LaMa/reference/tpm.md). When passed
an array of transition matrices of dimension `c(N,N,K)`,
[`stationary()`](https://janolefi.github.io/LaMa/reference/stationary.md)
returns a matrix of dimension `c(K,N)` with the stationary distribution
for each tpm in its rows. Both are then passed to
[`forward()`](https://janolefi.github.io/LaMa/reference/forward.md)
together with `trackID`.

``` r

nll_partial = function(par) {
  getAll(par, dat)
  Gamma = tpm(beta, matrix(z, ncol = 1)) # array of K individual-specific tpms
  Delta = stationary(Gamma) # when passed an array, returns a matrix of dim c(K, N),
                            # each row being the stationary distribution for one tpm
  mu = exp(log_mu); REPORT(mu)
  sigma = exp(log_sigma); REPORT(sigma)
  REPORT(beta)
  allprobs = matrix(1, length(x), N)
  for(j in 1:N) allprobs[,j] = dnorm(x, mu[j], sigma[j])
  # Delta matrix and Gamma array supply individual-specific starts and tpms
  -forward(Delta, Gamma, allprobs, trackID)
}
```

### Estimating the model

``` r

trackID = rep(1:K, each = n)

par = list(
  beta = matrix(c(-2, -2, 0, 0), nrow = N*(N-1)), # regression coefficients for tpm
  log_mu = log(c(15, 60)),                         # initial state-dependent means
  log_sigma = log(c(10, 40))                       # initial state-dependent sds
)

dat = list(
  x = x,
  z = z,
  trackID = trackID,
  N = N
)

obj_partial = MakeADFun(nll_partial, par, silent = TRUE)
system.time(
  opt_partial <- nlminb(obj_partial$par, obj_partial$fn, obj_partial$gr)
)
#>    user  system elapsed 
#>   0.018   0.000   0.018
mod_partial = report(obj_partial)
```

The estimated regression coefficients for the transition probabilities
and the shared state-dependent parameters:

``` r

mod_partial$beta   # true: c(-2, -2, 1, -1)
#>           [,1]       [,2]
#> [1,] -2.055222  0.6029733
#> [2,] -2.088381 -1.2557279
mod_partial$mu     # true: 15, 60
#> [1] 15.70795 57.88094
mod_partial$sigma  # true: 10, 40
#> [1] 10.33312 41.92609
```

## Random effects

In the partial pooling example above, we estimated a separate transition
probability matrix for each individual based on a fixed covariate. An
alternative and often more parsimonious approach is to treat the
individual-level variation in transition probabilities as arising from
**Gaussian random effects**. Instead of estimating K independent
transition structures, we model each individual’s linear predictor as a
deviation from a shared population-level intercept:
\text{logit}(\gamma\_{ij}^{(k)}) = \beta_0^{(ij)} + u\_{ij}^{(k)}, \quad
u\_{ij}^{(k)} \overset{\text{iid}}{\sim} N\bigl(0,
\sigma_u^{(ij)}\bigr), for each off-diagonal element (i \neq j) and each
individual k = 1, \dots, K. The population-level intercept
\beta_0^{(ij)} captures the average switching tendency, while
\sigma_u^{(ij)} quantifies how much individuals vary around it. The
number of state-process parameters is 2 \times N(N-1) — two intercepts
and two standard deviations for a two-state model — regardless of how
many individuals are in the data.

To fit this model we need to integrate out the random effects from the
joint likelihood. `RTMB` handles this automatically when the random
effects are declared via the `random` argument of `MakeADFun()`.

### Generating data

``` r

set.seed(123)
K = 50    # number of individuals
N = 2     # number of states
n = 100   # observations per individual

# Shared state-dependent parameters
mu = c(15, 60)
sigma = c(10, 40)

# Population-level intercepts and individual random effects
beta_0 = c(-2, -2)  # one per off-diagonal element
sigma_u = c(1, 1)   # one per off-diagonal element

u_true = matrix(rnorm(K * N*(N-1), 0, rep(sigma_u, each = K)),
                nrow = K, ncol = N*(N-1))

# K individual-specific tpms and their stationary distributions
Gamma_true = array(NA, c(N, N, K))
for(k in 1:K) Gamma_true[,,k] = tpm(beta_0 + u_true[k,])
Delta_true = stationary(Gamma_true)

# Simulate all tracks
s = x = rep(NA, n * K)
for(k in 1:K){
  sk = xk = rep(NA, n)
  sk[1] = sample(1:2, 1, prob = Delta_true[k,])
  xk[1] = rnorm(1, mu[sk[1]], sigma[sk[1]])
  for(t in 2:n){
    sk[t] = sample(1:2, 1, prob = Gamma_true[sk[t-1],,k])
    xk[t] = rnorm(1, mu[sk[t]], sigma[sk[t]])
  }
  s[(k-1)*n + 1:n] = sk
  x[(k-1)*n + 1:n] = xk
}

trackID = rep(1:K, each = n)
```

### Writing the joint negative log-likelihood

The joint negative log-likelihood includes both the data likelihood (via
[`forward()`](https://janolefi.github.io/LaMa/reference/forward.md)) and
the Gaussian prior on the random effects:

``` r

jnll_re = function(par) {
  getAll(par, dat)

  # Random effect standard deviations
  sigma_u = exp(log_sigma_u); REPORT(sigma_u); ADREPORT(sigma_u)

  # Gaussian prior on u: u[,ij] ~ N(0, sigma_u[ij])
  jnll = 0
  for(ij in 1:(N*(N-1))) {
    jnll = jnll - sum(dnorm(u[,ij], 0, sigma_u[ij], log = TRUE))
  }

  # Individual-specific tpms from population intercept + random effect
  Gamma = AD(array(0, c(N, N, K)))
  for(k in 1:K) Gamma[,,k] = tpm(beta_0 + u[k,])
  Delta = stationary(Gamma)

  # Observation model
  mu = exp(log_mu); REPORT(mu)
  sigma = exp(log_sigma); REPORT(sigma)
  allprobs = matrix(1, length(x), N)
  for(j in 1:N) allprobs[,j] = dnorm(x, mu[j], sigma[j])

  jnll - forward(Delta, Gamma, allprobs, trackID)
}
```

### Estimating the model

The random effects matrix `u` is declared via `random = "u"`, telling
`RTMB` to integrate it out. The optimiser only sees the marginal
(random-effects-free) likelihood as a function of the fixed-effect
parameters `beta_0`, `log_sigma_u`, `log_mu`, and `log_sigma`.

``` r

par = list(
  beta_0 = rep(-2, N*(N-1)),      # population-level intercepts
  log_sigma_u = log(c(1, 1)),     # log-SDs of random effects
  u = matrix(0, K, N*(N-1)),      # random effects (integrated out)
  log_mu = log(c(15, 60)),        # initial state-dependent means
  log_sigma = log(c(10, 40))      # initial state-dependent sds
)

dat = list(
  x = x,
  trackID = trackID,
  K = K,
  N = N
)

obj_re = MakeADFun(jnll_re, par, random = "u", silent = TRUE)
system.time(
  opt_re <- nlminb(obj_re$par, obj_re$fn, obj_re$gr)
)
#>    user  system elapsed 
#>   4.036   0.017   4.053
mod_re = report(obj_re)
```

The estimated standard deviations of the random effects tell us how much
the transition probabilities vary across individuals. We can check the
key population-level estimates against the true values used in the
simulation:

``` r

mod_re$sigma_u  # true: 1, 1
#> [1] 1.073841 0.934637
mod_re$mu       # true: 15, 60
#> [1] 14.62212 60.00429
mod_re$sigma    # true: 10, 40
#> [1]  9.937914 39.891970
```

Standard errors for all fixed effects and
[`ADREPORT()`](https://rdrr.io/pkg/RTMB/man/TMB-interface.html)ed
quantities are available via `sdreport()`:

``` r

sdr_re = sdreport(obj_re)
summary(sdr_re, "report") # sigma_u with SEs
#>         Estimate Std. Error
#> sigma_u 1.073841  0.1693166
#> sigma_u 0.934637  0.1572146
```

> Continue reading with [**Penalised
> splines**](https://janolefi.github.io/LaMa/articles/Penalised_splines.html).
