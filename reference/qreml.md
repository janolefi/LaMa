# Automatic smoothness selection for arbitrary statistical models with penalised splines, simple random effects, or tensorproducts

Efficiently fits models involving quadratic penalties of the form
\$\$\sum\_{i} \lambda_i b^T S_i b,\$\$ by combining automatic
differentiation via `RTMB` with a custom implementation of the extended
Fellner-Schall update.

Users only need to supply a custom penalised log-likelihood function
that calls
[`penalty()`](https://janolefi.github.io/LaMa/reference/penalty.md) or
[`penalty2()`](https://janolefi.github.io/LaMa/reference/penalty2.md).

## Usage

``` r
qreml(
  pnll,
  par,
  dat,
  random,
  map = NULL,
  silent = 1,
  spname = "lambda",
  tol_edf = 0.01,
  maxiter = 100,
  alpha = 0.1,
  tol = 0.01,
  smoothing = 1,
  method = "BFGS",
  lsp_max = 15,
  step_small = 0.05,
  max_halve = 6,
  control = list(),
  spHess = FALSE,
  joint_unc = FALSE,
  saveall = FALSE
)
```

## Arguments

- pnll:

  penalised negative log-likelihood function that is structured as
  dictated by `RTMB` and uses the
  [`penalty`](https://janolefi.github.io/LaMa/reference/penalty.md) or
  [`penalty2`](https://janolefi.github.io/LaMa/reference/penalty2.md)
  function to compute the penalty

- par:

  named list of initial parameters

- dat:

  initial data list that contains the data used in the likelihood
  function, hyperparameters, and the **initial penalty strength** vector

- random:

  vector of names of the random effects/ penalised parameters in `par`

  **Caution:** The ordering of `random` needs to match the order of the
  random effects passed to `penalty`.

- map:

  optional map argument, containing factor vectors to indicate parameter
  sharing or fixing

- silent:

  integer silencing level: 0 corresponds to full printing of inner and
  outer iterations, 1 to printing of outer iterations only, and 2 to no
  printing

- spname:

  optional name given to the penalty strength parameter in `dat`.
  Defaults to `"lambda"`

- tol_edf:

  convergence tolerance on the **effective degrees of freedom**.
  Defaults to 0.001.

  This is the primary convergence criterion. The iteration stops once no
  smooth's effective degrees of freedom has changed by more than
  `tol_edf` across the last four outer iterations, provided the step in
  `log(lambda)` is also small. The effective degrees of freedom are what
  say whether the **fitted smooth** is still changing, and unlike a
  tolerance on the restricted likelihood they mean the same thing on
  every model, being measured in effective parameters rather than in
  nats. The penalty strengths themselves are a poor proxy: `lambda` can
  slide along a flat ridge for many iterations, moving substantially in
  relative terms, while neither the criterion nor the fit changes
  appreciably.

- maxiter:

  maximum number of outer iterations

- alpha:

  smallest factor by which a penalty strength may **decrease** in one
  outer iteration, a number in \[0, 1). Defaults to 0.1.

  Penalty strengths are free to increase as fast as the update proposes,
  but cannot collapse faster than this per iteration. Reducing a penalty
  strength too quickly can push the inner optimisation into a local
  optimum or a numerically awkward region, which matters more here than
  in a GAM because the likelihood is user-written. Set to zero to remove
  the floor entirely. Step length is handled separately, by `max_halve`
  and the adaptive multiplier.

- tol:

  **fallback** convergence tolerance, on the restricted log-likelihood.
  Defaults to 0.01.

  Used only when the effective degrees of freedom are not trustworthy,
  i.e.\\ when they fall outside \\\[0, K_i\]\\, which happens when the
  data Hessian is indefinite (see `tol_edf` and the returned
  `edf_valid`). The iteration then stops once the restricted
  log-likelihood has changed by less than `tol` over the last four outer
  iterations. A tolerance in nats is not comparable across models, which
  is why it is the fallback rather than the primary criterion.

- smoothing:

  optional scaling factor for the final penalty strength parameters.
  Increasing this beyond one leads to a smoother final model

- method:

  inner optimisation method to be used by
  [`optim`](https://rdrr.io/r/stats/optim.html). Defaults to `"BFGS"`

- lsp_max:

  largest value allowed for `log(lambda)`. Defaults to 15, as in `mgcv`,
  i.e. penalty strengths saturate at roughly 3.3e6

- step_small:

  size of a step in `log(lambda)` below which the step multiplier is
  allowed to double. Defaults to 0.05, as in `mgcv`

- max_halve:

  maximum number of times a step that decreases the restricted
  likelihood is halved before it is accepted anyway. Defaults to 6.

  `mgcv` never shortens below the full Fellner-Schall step and accepts a
  worse one instead, which for a user-written likelihood can drift
  downhill for tens of iterations.

- control:

  list of control parameters for
  [`optim`](https://rdrr.io/r/stats/optim.html) to use in the inner
  optimisation

- spHess:

  logical, if `TRUE`, the sparse automatic differentiation Hessian is
  used for evaluation. The factorisation is dense either way

- joint_unc:

  logical, if `TRUE`, joint `RTMB` object is returned allowing for joint
  uncertainty quantification

- saveall:

  logical, if `TRUE`, then all model objects from each iteration are
  saved in the final model object

## Value

model object of class `"qremlModel"`, carrying the fitted quantities,
the smoothness selection diagnostics shown by
[`summary.qremlModel`](https://janolefi.github.io/LaMa/reference/summary.qremlModel.md),
and `outer_hessian()`, the outer Hessian that
[`sdreport_outer`](https://janolefi.github.io/LaMa/reference/sdreport_outer.md)
turns into standard errors

## Details

Step size control follows the extended Fellner-Schall implementation in
`mgcv` while convergence is judged on the effective degrees of freedom
(differing from `mgcv`).

**The criterion.** The penalty strengths are chosen to maximise the
Laplace-approximate restricted log-likelihood \$\$V(\lambda) =
\ell_p(\hat{b}; \lambda) + \tfrac{1}{2} \log \|S\_\lambda\|\_+ -
\tfrac{1}{2} \log \|J\|,\$\$ where \\\ell_p\\ is the penalised
log-likelihood at its mode \\\hat{b}\\, \\S\_\lambda = \sum_i \lambda_i
S_i\\, \\\|\cdot\|\_+\\ denotes the product of the positive eigenvalues,
and \\J\\ is the Hessian of the penalised negative log-likelihood at
\\\hat{b}\\. Every outer iteration refits the inner problem to obtain
\\\hat{b}\\ and \\J\\, so one call is a sequence of complete model fits,
each warm started at the previous one.

**The update.** The penalty strengths are updated multiplicatively,
\\\lambda_i \leftarrow \lambda_i r_i\\, by the extended Fellner-Schall
ratio \$\$r_i = \frac{\mathrm{tr}(S\_\lambda^- S_i) - \mathrm{tr}(J^{-1}
S_i)}{\hat{b}^T S_i \hat{b}},\$\$ which equals one exactly where
\\\partial V / \partial \lambda_i = 0\\ and otherwise moves
\\\lambda_i\\ in the direction of the gradient, without ever making it
negative. For a smooth carrying a single penalty the first trace is
\\\mathrm{rank}(S_i) / \lambda_i\\, obtained once from an
eigendecomposition instead of at every iteration. Only the diagonal
blocks of \\J^{-1}\\ enter, and these come from the Cholesky factor by
triangular solves rather than from a full inverse.

**Step control.** The step is taken in \\\log \lambda\\ and scaled by a
multiplier that persists across iterations: it doubles while the
criterion keeps improving and is halved, up to `max_halve` times, while
it does not. A step that still worsens the criterion after that many
halvings is accepted, which is the only way the iteration can move
downhill. In addition no penalty strength may fall by more than a factor
`alpha` in one iteration, and \\\log \lambda\\ is capped at `lsp_max`.
The fit returned is the best iterate seen, not necessarily the last one.

**Convergence.** The iteration stops once no smooth's effective degrees
of freedom \$\$\mathrm{edf}\_i = K_i - \mathrm{tr}(J^{-1} \lambda_i
S_i)\$\$ has moved by more than `tol_edf` over the last four iterations,
provided the step in \\\log \lambda\\ is also small; see `tol_edf`.
Being eigenvalue sums, the effective degrees of freedom must lie in
\\\[0, K_i\]\\, and they leave that range only if the data Hessian is
indefinite, in which case the traces driving the update are not the
quantities the method assumes. That is reported rather than repaired,
through the returned `edf_valid` and a warning, and the criterion then
falls back to `tol` on the restricted log-likelihood.

**Numerics.** \\J\\ is the exact automatic differentiation Hessian,
which matters because the finite differencing error of
[`optimHess`](https://rdrr.io/r/stats/optim.html) is of the same order
as the criterion differences the convergence test has to resolve. It is
factorised by a pivoted Cholesky decomposition after scaling to a unit
diagonal; should it not be positive definite, a ridge is added and
reported through the returned `hessian_repaired`.

**Uncertainty.** The returned `outer_hessian()` gives the Hessian of the
negative restricted log-likelihood in the penalty strengths, which
[`sdreport_outer`](https://janolefi.github.io/LaMa/reference/sdreport_outer.md)
turns into standard errors. Note that these describe the curvature of
\\V\\ only, and treat the approximation to it as exact.

## References

Wood, S. N., & Fasiolo, M. (2017). A generalized Fellner-Schall method
for smoothing parameter optimization with application to Tweedie
location, scale and shape models. Biometrics, 73(4), 1071-1081.

Koslik, J. O. (2024). Efficient smoothness selection for nonparametric
Markov-switching models via quasi restricted maximum likelihood. arXiv
preprint arXiv:2411.11498.

## See also

[`penalty`](https://janolefi.github.io/LaMa/reference/penalty.md) and
[`penalty2`](https://janolefi.github.io/LaMa/reference/penalty2.md) to
compute the penalty inside the likelihood function, and
[`qreml_old`](https://janolefi.github.io/LaMa/reference/qreml_old.md)
for the original algorithm

## Examples

``` r
data = trex[1:2000,] # subset of the data

# initial parameter list
par = list(logmu = log(c(0.3, 2.5)), # step mean
           logsigma = log(c(0.3, 1.5)), # step sd
           beta0 = c(-2,-2), # state process intercept
           beta_spline = matrix(rep(0, 18), nrow = 2)) # state process spline coefs
          
# data object with initial penalty strength lambda
dat = list(step = data$step, # step length
           tod = data$tod, # time of day covariate
           N = 2, # number of states
           lambda = rep(20,2)) # initial penalty strength

# building model matrices
modmat = make_matrices(~ s(tod, bs = "cp"), 
                       data = data.frame(tod = 1:24), 
                       knots = list(tod = c(0,24))) # wrapping points
dat$Z = modmat$Z # spline design matrix
dat$S = modmat$S # penalty matrix

# penalised negative log-likelihood function
pnll = function(par) {
  getAll(par, dat) # makes everything contained available without $
  Gamma = tpm(cbind(beta0, beta_spline), Z) # transition probabilities
  delta = stationary_p(Gamma, t = 1) # initial distribution
  mu = exp(logmu) # step mean
  sigma = exp(logsigma) # step sd
  # calculating all state-dependent densities
  allprobs = matrix(1, nrow = length(step), ncol = N)
  ind = which(!is.na(step)) # only for non-NA obs.
  for(j in 1:N) allprobs[ind,j] = dgamma2(step[ind],mu[j],sigma[j])
  -forward(delta, Gamma[,,tod], allprobs) +
      penalty(beta_spline, S, lambda) # this does all the penalization work
}

# model fitting
mod = qreml(pnll, par, dat, random = "beta_spline", silent = 2)
```
