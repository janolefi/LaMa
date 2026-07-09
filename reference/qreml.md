# Quasi restricted maximum likelihood (qREML) algorithm for models with penalised splines or simple i.i.d. random effects

This algorithm can be used very flexibly to fit statistical models that
involve **penalised splines** or simple **i.i.d. random effects**, i.e.
that have penalties of the form \$\$0.5 \sum\_{i} \lambda_i b_i^T S_i
b_i,\$\$ with smoothing parameters \\\lambda_i\\, coefficient vectors
\\b_i\\, and fixed penalty matrices \\S_i\\.

The **qREML** algorithm is typically much faster than REML or marginal
ML using the full Laplace approximation method, but may be slightly less
accurate regarding the estimation of the penalty strength parameters.

Under the hood, `qreml` uses the R package `RTMB` for automatic
differentiation in the inner optimisation. The user has to specify the
**penalised negative log-likelihood function** `pnll` structured as
dictated by `RTMB` and use the
[`penalty`](https://janolefi.github.io/LaMa/reference/penalty.md)
function to compute the quadratic-form penalty inside the likelihood.

## Usage

``` r
qreml(
  pnll,
  par,
  dat,
  random,
  map = NULL,
  silent = 1,
  psname = "lambda",
  alpha = 0.3,
  smoothing = 1,
  maxiter = 100,
  tol = 1e-04,
  method = "BFGS",
  control = list(),
  conv_crit = "relchange",
  spHess = FALSE,
  joint_unc = FALSE,
  saveall = FALSE
)
```

## Arguments

- pnll:

  penalised negative log-likelihood function that is structured as
  dictated by `RTMB` and uses the
  [`penalty`](https://janolefi.github.io/LaMa/reference/penalty.md)
  function from `LaMa` to compute the penalty

  Needs to be a function of the named list of initial parameters `par`
  only.

- par:

  named list of initial parameters

  The random effects/ spline coefficients can be vectors or matrices,
  the latter summarising several random effects of the same structure,
  each one being a row in the matrix.

- dat:

  initial data list that contains the data used in the likelihood
  function, hyperparameters, and the **initial penalty strength** vector

  If the initial penalty strength vector is **not** called `lambda`, the
  name it has in `dat` needs to be specified using the `psname` argument
  below. Its length needs to match the to the total number of random
  effects.

- random:

  vector of names of the random effects/ penalised parameters in `par`

  **Caution:** The ordering of `random` needs to match the order of the
  random effects passed to
  [`penalty`](https://janolefi.github.io/LaMa/reference/penalty.md)
  inside the likelihood function.

- map:

  optional map argument, containing factor vectors to indicate parameter
  sharing or fixing.

  Needs to be a named list for a subset of fixed effect parameters or
  penalty strength parameters. For example, if the model has four
  penalty strength parameters, `map[[psname]]` could be
  `factor(c(NA, 1, 1, 2))` to fix the first penalty strength parameter,
  estimate the second and third jointly, and estimate the fourth
  separately.

- silent:

  integer silencing level: 0 corresponds to full printing of inner and
  outer iterations, 1 to printing of outer iterations only, and 2 to no
  printing.

- psname:

  optional name given to the penalty strength parameter in `dat`.
  Defaults to `"lambda"`.

- alpha:

  optional hyperparamater for exponential smoothing of the penalty
  strengths.

  For larger values smoother convergence is to be expected but the
  algorithm may need more iterations.

- smoothing:

  optional scaling factor for the final penalty strength parameters

  Increasing this beyond one will lead to a smoother final model. Can be
  an integer or a vector of length equal to the length of the penalty
  strength parameter.

- maxiter:

  maximum number of iterations in the outer optimisation over the
  penalty strength parameters.

- tol:

  Convergence tolerance for the penalty strength parameters.

- method:

  optimisation method to be used by
  [`optim`](https://rdrr.io/r/stats/optim.html). Defaults to `"BFGS"`,
  but might be changed to `"L-BFGS-B"` for high-dimensional settings.

- control:

  list of control parameters for
  [`optim`](https://rdrr.io/r/stats/optim.html) to use in the inner
  optimisation. Here, `optim` uses the `BFGS` method which cannot be
  changed.

  We advise against changing the default values of `reltol` and `maxit`
  as this can decrease the accuracy of the Laplace approximation.

- conv_crit:

  character, convergence criterion for the penalty strength parameters.
  Can be `"relchange"` (default) or `"gradient"`.

- spHess:

  logical, if `TRUE`, sparse AD Hessian is used in each outer iteration.
  If your Hessian is large and sparse (many cross derivatives are 0),
  this will speed up the computations a lot. If your Hessian is dense,
  this will slow down the computations slightly and might require
  significantly more memory.

- joint_unc:

  logical, if `TRUE`, joint `RTMB` object is returned allowing for joint
  uncertainty quantification

- saveall:

  logical, if `TRUE`, then all model objects from each iteration are
  saved in the final model object.

## Value

model object of class 'qremlModel'. This is a list containing:

- ...:

  everything that is reported inside `pnll` using
  [`RTMB::REPORT()`](https://rdrr.io/pkg/RTMB/man/TMB-interface.html).
  When using `forward`, `tpm_g`, etc., this may involve automatically
  reported objects.

- obj:

  `RTMB` AD object containing the final conditional model fit

- psname:

  final penalty strength parameter vector

- all_psname:

  list of all penalty strength parameter vectors over the iterations

- par:

  named estimated parameter list in the same structure as the initial
  `par`. Note that the name `par` is not fixed but depends on the
  original name of your `par` list.

- relist_par:

  function to convert the estimated parameter vector to the estimated
  parameter list. This is useful for uncertainty quantification based on
  sampling from a multivariate normal distribution.

- par_vec:

  estimated parameter vector

- llk:

  unpenalised log-likelihood at the optimum

- n_fixpar:

  number of fixed, i.e. unpenalised, parameters

- edf:

  overall effective number of parameters

- all_edf:

  list of effective number of parameters for each smooth

- Hessian_condtional:

  final Hessian of the conditional penalised fit

- obj_joint:

  if `joint_unc = TRUE`, joint `RTMB` object for joint uncertainty
  quantification in model and penalty parameters.

## References

Koslik, J. O. (2024). Efficient smoothness selection for nonparametric
Markov-switching models via quasi restricted maximum likelihood. arXiv
preprint arXiv:2411.11498.

## See also

[`penalty`](https://janolefi.github.io/LaMa/reference/penalty.md) and
[`penalty2`](https://janolefi.github.io/LaMa/reference/penalty2.md) to
compute the penalty inside the likelihood function

## Examples

``` r
data = trex[1:1000,] # subset

# initial parameter list
par = list(logmu = log(c(0.3, 2.5)), # step mean
           logsigma = log(c(0.3, 1.5)), # step sd
           beta0 = c(-2,-2), # state process intercept
           betaspline = matrix(rep(0, 18), nrow = 2)) # state process spline coefs
          
# data object with initial penalty strength lambda
dat = list(step = data$step, # step length
           tod = data$tod, # time of day covariate
           N = 2, # number of states
           lambda = rep(10,2)) # initial penalty strength

# building model matrices
modmat = make_matrices(~ s(tod, bs = "cp"), 
                       data = data.frame(tod = 1:24), 
                       knots = list(tod = c(0,24))) # wrapping points
dat$Z = modmat$Z # spline design matrix
dat$S = modmat$S # penalty matrix

# penalised negative log-likelihood function
pnll = function(par) {
  getAll(par, dat) # makes everything contained available without $
  Gamma = tpm_g(Z, cbind(beta0, betaspline), ad = TRUE) # transition probabilities
  delta = stationary_p(Gamma, t = 1, ad = TRUE) # initial distribution
  mu = exp(logmu) # step mean
  sigma = exp(logsigma) # step sd
  # calculating all state-dependent densities
  allprobs = matrix(1, nrow = length(step), ncol = N)
  ind = which(!is.na(step)) # only for non-NA obs.
  for(j in 1:N) allprobs[ind,j] = dgamma2(step[ind],mu[j],sigma[j])
  -forward_g(delta, Gamma[,,tod], allprobs) +
      penalty(betaspline, S, lambda) # this does all the penalization work
}

# model fitting
mod = qreml(pnll, par, dat, random = "betaspline", silent = 2)
```
