# Variance-covariance matrix of a `LaMaModel` object

Computes the approximate covariance matrix of the estimated parameters
`coef(object)`, either as the inverse Hessian of the negative
log-likelihood at the estimate (`type = "hessian"`), or as a sandwich
that remains valid if the model is misspecified: with the outer product
of the scores (`type = "sandwich"`), heteroskedasticity and
autocorrelation consistent (`type = "HAC"`) or cluster-robust over
tracks (`type = "cluster"`).

## Usage

``` r
# S3 method for class 'LaMaModel'
vcov(
  object,
  type = c("hessian", "sandwich", "HAC", "cluster"),
  kernel = c("QS", "Bartlett"),
  bw = NULL,
  cluster = NULL,
  cr = c("CR1", "CR0", "CR3", "CR3onestep"),
  ...
)
```

## Arguments

- object:

  model object of class `"LaMaModel"`

- type:

  `"hessian"` (default) for the inverse Hessian, `"sandwich"` for the
  sandwich with the outer product of the scores, `"HAC"` for the HAC
  sandwich, or `"cluster"` for the cluster-robust sandwich

- kernel:

  kernel of the HAC sandwich: `"QS"` (Quadratic Spectral, default) or
  `"Bartlett"` (Newey-West), see Details

- bw:

  bandwidth of the HAC sandwich, a non-negative number. Defaults to the
  data-driven bandwidth of Andrews (1991) for the chosen kernel, see
  Details. For the Bartlett kernel, `bw = L + 1` includes the lags up to
  \\L\\, and `bw` at most 1 gives the same as `type = "sandwich"`.

- cluster:

  cluster identifier for `type = "cluster"`, one value per observation
  (of any type, in the order of the observations in
  [`forward`](https://janolefi.github.io/LaMa/reference/forward.md)).
  Defaults to the tracks, i.e. the `trackID` reported by
  [`forward`](https://janolefi.github.io/LaMa/reference/forward.md).

- cr:

  small-sample adjustment for `type = "cluster"`: `"CR1"` (default)
  multiplies the meat by \\G / (G - 1)\\, `"CR0"` does not adjust,
  `"CR3"` is the jackknife, which refits the model without each cluster,
  and `"CR3onestep"` its approximation by one Newton step per cluster,
  see Details.

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

**Sandwich.** With the per-observation score contributions \\u_t\\ from
[`scores`](https://janolefi.github.io/LaMa/reference/scores.md), the
covariance is \\H^{-1} M H^{-1}\\ with the outer product of the scores
as meat, \\M = \sum_t u_t u_t^\top\\ (the outer product of gradients,
OPG). It is robust to a misspecified predictive distribution as long as
the \\u_t\\ are serially uncorrelated, which they are under correct
specification (martingale differences), but not to dependence the model
misses. It equals the HAC sandwich with the Bartlett kernel and `bw` at
most 1. Not to be confused with the inverse of the outer product alone,
\\M^{-1}\\ (the BHHH estimator), which, like \\H^{-1}\\, is only valid
under correct specification.

**HAC sandwich.** With the per-observation score contributions \\u_t\\
from [`scores`](https://janolefi.github.io/LaMa/reference/scores.md),
the covariance is \\H^{-1} M H^{-1}\\, with the kernel estimate of the
variance of the score sum \$\$M = \Gamma_0 + \sum\_{l \ge 1} k(l / b)
(\Gamma_l + \Gamma_l^\top), \qquad \Gamma_l = \sum_t u_t
u\_{t-l}^\top,\$\$ where the sum in \\\Gamma_l\\ only runs over pairs of
observations from the same track (as reported by
[`forward`](https://janolefi.github.io/LaMa/reference/forward.md)), as
different tracks are independent. Missing observations contribute zero
scores, so lags remain in units of time. The kernel \\k\\ down-weights
the autocovariances at longer lags, which are estimated with much noise,
and the bandwidth \\b\\ (`bw`) scales the lag axis. Two kernels are
available:

- `kernel = "QS"` (default), the Quadratic Spectral kernel \$\$k(x) =
  \frac{25}{12 \pi^2 x^2} \Bigl(\frac{\sin(6 \pi x / 5)}{6 \pi x / 5} -
  \cos(6 \pi x / 5)\Bigr).\$\$ Its weights are positive up to lag \\1.2
  b\\ and then oscillate with small amplitude, so all lags contribute.
  It is the most efficient kernel among those that guarantee a positive
  semi-definite \\M\\ (Andrews, 1991): it keeps the short lags at nearly
  full weight, which reduces the bias.

- `kernel = "Bartlett"`, \\k(x) = \max(1 - \|x\|, 0)\\, the Newey-West
  estimator: with \\b = L + 1\\, lags up to \\L\\ enter with weights
  \\1 - l / (L + 1)\\.

The default bandwidth is the data-driven choice of Andrews (1991), which
minimises the asymptotic mean squared error of \\M\\ for an AR(1)
approximation of the scores: with \\\rho_a\\ and \\\sigma_a^2\\ the
coefficient and innovation variance of an AR(1) with intercept fitted to
each score column \\a\\ (by least squares on consecutive observations
within tracks), \$\$b = 1.3221 \\ (\alpha_2 \bar T)^{1/5} \\
\text{(QS)}, \qquad b = 1.1447 \\ (\alpha_1 \bar T)^{1/3} \\
\text{(Bartlett)},\$\$ with \\\bar T\\ the mean track length, \\\alpha_2
= \sum_a 4 \rho_a^2 \sigma_a^4 (1 - \rho_a)^{-8} / D\\, \\\alpha_1 =
\sum_a 4 \rho_a^2 \sigma_a^4 (1 - \rho_a)^{-6} (1 + \rho_a)^{-2} / D\\
and \\D = \sum_a \sigma_a^4 (1 - \rho_a)^{-4}\\. For a single track,
this is
[`sandwich::bwAndrews()`](https://zeileis.codeberg.page/sandwich/reference/weightsAndrews.html)
without prewhitening and with equal weights. All score columns are
weighted equally, so the bandwidth is driven by the most persistent and
most variable ones. Stronger serial dependence of the scores gives a
larger bandwidth.

All lags are computed exactly, with the fast Fourier transform within
each track, which takes \\O(p T \log T + T p^2)\\ operations for \\T\\
observations and \\p\\ parameters instead of \\O(T^2 p^2)\\ for summing
all lagged products directly.

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

`cr = "CR3"` gives the jackknife estimator (MacKinnon, Nielsen and Webb,
2023) \$\$V = \frac{G - 1}{G} \sum\_{k=1}^G d_k d_k^\top, \qquad d_k =
\hat\theta - \hat\theta^{(k)},\$\$ where \\\hat\theta^{(k)}\\ is the
estimate without cluster \\k\\. It is found by minimising the weighted
negative log-likelihood with weight zero on the observations of cluster
\\k\\, started at \\\hat\theta\\, so the model is refitted once per
cluster, and a message says so. Compared to CR1, it inflates the
contributions of clusters with high leverage (e.g. long tracks), whose
score sums are shrunk the most by fitting. Refits that do not converge
give a warning. If the Hessian at the estimate of a refit is not
positive definite, as for parameters not or only weakly identified
without that cluster, an error is issued. That a refit has moved to
another mode, for example with a different order of the states, is not
checked; standard errors much larger than those of the one-step version
can be a sign of this.

`cr = "CR3onestep"` approximates the change \\d_k\\ by one Newton step
from the estimate, \$\$d_k = (H - H_k)^{-1} s_k,\$\$ where \\H - H_k\\
is the Hessian without cluster \\k\\, from the weighted likelihood with
weight zero on its observations. It costs one Hessian per cluster and no
refits. For linear least squares, the step is exact; otherwise, the
approximation is poorest if there are few clusters or if one of them
carries much of the information. With equal leverage \\H_k = H / G\\, it
coincides with CR1. The one-step deletion is a pure Newton step, so it
requires \\H - H_k\\ to be positive definite: if leaving out a cluster
makes it singular or indefinite, an error is issued.

All sandwiches require the likelihood function to call
[`forward`](https://janolefi.github.io/LaMa/reference/forward.md) once
(with `trackID` for several tracks) and are not available for models
with random effects or from
[`qreml`](https://janolefi.github.io/LaMa/reference/qreml.md).

## References

Andrews, D. W. K. (1991). Heteroskedasticity and autocorrelation
consistent covariance matrix estimation. *Econometrica*, 59(3), 817–858.

Newey, W. K. and West, K. D. (1987). A simple, positive semi-definite,
heteroskedasticity and autocorrelation consistent covariance matrix.
*Econometrica*, 55(3), 703–708.

MacKinnon, J. G., Nielsen, M. O. and Webb, M. D. (2023). Cluster-robust
inference: A guide to empirical practice. *Journal of Econometrics*,
232(2), 272–299.

Cameron, A. C. and Miller, D. L. (2015). A practitioner's guide to
cluster-robust inference. *Journal of Human Resources*, 50(2), 317–372.

## See also

[`coef.LaMaModel`](https://janolefi.github.io/LaMa/reference/coef.LaMaModel.md),
[`scores`](https://janolefi.github.io/LaMa/reference/scores.md)
