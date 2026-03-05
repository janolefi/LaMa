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
[`penalty`](https://janoleko.github.io/reference/penalty.md) function to
compute the quadratic-form penalty inside the likelihood.

## Usage

``` r
qreml_old(
  pnll,
  par,
  dat,
  random,
  map = NULL,
  psname = "lambda",
  alpha = 0.25,
  smoothing = 1,
  maxiter = 100,
  tol = 1e-04,
  control = list(reltol = 1e-10, maxit = 1000),
  silent = 1,
  joint_unc = TRUE,
  saveall = FALSE
)
```

## Arguments

- pnll:

  penalised negative log-likelihood function that is structured as
  dictated by `RTMB` and uses the
  [`penalty`](https://janoleko.github.io/reference/penalty.md) function
  from `LaMa` to compute the penalty

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
  [`penalty`](https://janoleko.github.io/reference/penalty.md) inside
  the likelihood function.

- map:

  optional map argument, containing factor vectors to indicate parameter
  sharing or fixing.

  Needs to be a named list for a subset of fixed effect parameters or
  penalty strength parameters. For example, if the model has four
  penalty strength parameters, `map[[psname]]` could be
  `factor(c(NA, 1, 1, 2))` to fix the first penalty strength parameter,
  estimate the second and third jointly, and estimate the fourth
  separately.

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

- control:

  list of control parameters for
  [`optim`](https://rdrr.io/r/stats/optim.html) to use in the inner
  optimisation. Here, `optim` uses the `BFGS` method which cannot be
  changed.

  We advise against changing the default values of `reltol` and `maxit`
  as this can decrease the accuracy of the Laplace approximation.

- silent:

  integer silencing level: 0 corresponds to full printing of inner and
  outer iterations, 1 to printing of outer iterations only, and 2 to no
  printing.

- joint_unc:

  logical, if `TRUE`, joint `RTMB` object is returned allowing for joint
  uncertainty quantification

- saveall:

  logical, if `TRUE`, then all model objects from each iteration are
  saved in the final model object. \# @param epsilon vector of two
  values specifying the cycling detection parameters. If the relative
  change of the new penalty strength to the previous one is larger than
  `epsilon[1]` but the change to the one before is smaller than
  `epsilon[2]`, the algorithm will average the two last values to
  prevent cycling.

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

[`penalty`](https://janoleko.github.io/reference/penalty.md) to compute
the penalty inside the likelihood function

## Examples

``` r
# no example
```
