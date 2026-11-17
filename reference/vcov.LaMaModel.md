# Variance-covariance matrix of a `LaMaModel` object

Computes the approximate covariance matrix of the estimated parameters
`coef(object)`, either as the inverse Hessian of the negative
log-likelihood at the estimate (`type = "hessian"`), or as a sandwich
that remains valid if the model is misspecified: heteroskedasticity and
autocorrelation consistent (`type = "HAC"`) or cluster-robust over
tracks (`type = "cluster"`).

## Usage

``` r
# S3 method for class 'LaMaModel'
vcov(
  object,
  type = c("hessian", "HAC", "cluster"),
  lag = NULL,
  cluster = NULL,
  cr = c("CR1", "CR0"),
  ...
)
```

## Arguments

- object:

  model object of class `"LaMaModel"`

- type:

  `"hessian"` (default) for the inverse Hessian, `"HAC"` for the HAC
  sandwich, or `"cluster"` for the cluster-robust sandwich

- lag:

  bandwidth \\L\\ of the HAC sandwich, i.e. the maximum lag of the
  autocovariances included. Defaults to the rule of thumb above.
  `lag = 0` gives the ordinary (heteroskedasticity robust) sandwich.

- cluster:

  cluster identifier for `type = "cluster"`, one value per observation
  (of any type, in the order of the observations in
  [`forward`](https://janolefi.github.io/LaMa/reference/forward.md)).
  Defaults to the tracks, i.e. the `trackID` reported by
  [`forward`](https://janolefi.github.io/LaMa/reference/forward.md).

- cr:

  small-sample adjustment for `type = "cluster"`: `"CR1"` (default)
  multiplies the meat by \\G / (G - 1)\\, `"CR0"` does not adjust.

- ...:

  additional arguments (not used)

## Value

covariance matrix with rows and columns named like `coef(object)`

## Details

**Hessian.** The Hessian \\H\\ is

- the exact Hessian `obj$he()` for models without random effects,

- for models with random effects, the Hessian of the
  Laplace-approximated marginal likelihood w.r.t. the fixed parameters,
  by finite differences of `obj$gr()` (as in `sdreport`),

- for models from
  [`qreml`](https://janolefi.github.io/LaMa/reference/qreml.md), the
  Hessian of the penalised negative log-likelihood, including the
  penalised coefficients. The covariance is then conditional on the
  estimated smoothing parameters.

If the Hessian is not positive definite, which indicates that the model
has not converged, a warning is issued.

**HAC sandwich.** With the per-observation score contributions \\u_t\\
from [`scores`](https://janolefi.github.io/LaMa/reference/scores.md),
the covariance is \\H^{-1} M H^{-1}\\, with the Newey-West (Bartlett
kernel) estimate \$\$M = \Gamma_0 + \sum\_{l=1}^{L} \Bigl(1 -
\frac{l}{L+1}\Bigr) (\Gamma_l + \Gamma_l^\top), \qquad \Gamma_l = \sum_t
u_t u\_{t-l}^\top,\$\$ where the sum in \\\Gamma_l\\ only runs over
pairs of observations from the same track (as reported by
[`forward`](https://janolefi.github.io/LaMa/reference/forward.md)), as
different tracks are independent. Missing observations contribute zero
scores, so lags remain in units of time. The default bandwidth is \\L =
\lfloor 4 (\bar T / 100)^{2/9} \rfloor\\, with \\\bar T\\ the mean track
length.

If the model is correctly specified, the \\u_t\\ are martingale
differences, so \\\Gamma_l \approx 0\\ for \\l \ge 1\\ and \\M \approx
H\\: both types then agree approximately, and marked differences point
to misspecification.

**Cluster-robust sandwich.** With the score sums \\s_k = \sum\_{t \in k}
u_t\\ of the clusters \\k = 1, \dots, G\\ (by default the tracks), the
meat is \$\$M = c \sum\_{k=1}^G s_k s_k^\top,\$\$ with \\c = 1\\
(`cr = "CR0"`) or \\c = G / (G - 1)\\ (`cr = "CR1"`, default). It makes
no assumption about the dependence within clusters, so no bandwidth is
needed, but requires independent clusters and is only reliable for many
of them. As the \\s_k\\ sum to the gradient, which is zero at the
estimate, \\M\\ has rank at most \\G - 1\\, so the covariance is
singular if \\G - 1\\ is less than the number of parameters. For
parameters informed by a single cluster only, such as track-specific
initial distributions, it is not meaningful. Clusters may combine
several tracks (e.g. several tracks of one individual), but a track must
not be split between clusters.

Both sandwiches require the likelihood function to call
[`forward`](https://janolefi.github.io/LaMa/reference/forward.md) once
(with `trackID` for several tracks) and are not available for models
with random effects or from
[`qreml`](https://janolefi.github.io/LaMa/reference/qreml.md).

## References

Newey, W. K. and West, K. D. (1987). A simple, positive semi-definite,
heteroskedasticity and autocorrelation consistent covariance matrix.
*Econometrica*, 55(3), 703–708.

Cameron, A. C. and Miller, D. L. (2015). A practitioner's guide to
cluster-robust inference. *Journal of Human Resources*, 50(2), 317–372.

## See also

[`coef.LaMaModel`](https://janolefi.github.io/LaMa/reference/coef.LaMaModel.md),
[`scores`](https://janolefi.github.io/LaMa/reference/scores.md)
