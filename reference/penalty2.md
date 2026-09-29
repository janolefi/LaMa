# Computes generalised quadratic-form penalties

This function computes a quadratic penalty of the form \$\$0.5 \sum\_{i}
\lambda_i b^T S_i b,\$\$ with smoothing parameters \\\lambda_i\\,
coefficient vector \\b\\, and fixed penalty matrices \\S_i\\. This
generalises the
[`penalty`](https://janolefi.github.io/LaMa/reference/penalty.md) by
allowing subsets of the coefficient vector \\b\\ to be penalised
multiple times with different smoothing parameters, which is necessary
for **tensor products**, **functional random effects** or **adaptive
smoothing**.

It is intended to be used inside the **penalised negative log-likelihood
function** when fitting models with penalised splines or simple random
effects via **quasi restricted maximum likelihood** (qREML) with the
[`qreml`](https://janolefi.github.io/LaMa/reference/qreml.md) function.
For [`qreml`](https://janolefi.github.io/LaMa/reference/qreml.md) to
work, the likelihood function needs to be compatible with the `RTMB` R
package to enable automatic differentiation.

## Usage

``` r
penalty2(re_coef, S, lambda)
```

## Arguments

- re_coef:

  list of coefficient vectors/ matrices

  Each list entry corresponds to a different smooth/ random effect with
  its own associated penalty matrix or penalty-matrix list in `S`. When
  several smooths/ random effects of the same kind are present, it is
  convenient to pass them as a matrix, where each row corresponds to one
  smooth/ random effect. This way all rows can use the same penalty
  matrix.

- S:

  list of fixed penalty matrices matching the structure of `re_coef`.

  This means if `re_coef` is of length 3, then `S` needs to be a list of
  length 3. Each entry needs to be either a penalty matrix, matching the
  dimension of the corresponding entry in `re_coef`, or a list with
  multiple penalty matrices for tensor products.

- lambda:

  penalty strength parameter vector that has a length corresponding to
  the provided `re_coef` and `S`.

  Specifically, for entries with one penalty matrix,
  `nrow(re_coef[[i]])` parameters are needed. For entries with `k`
  penalty matrices, `k * nrow(re_coef[[i]])` parameters are needed.

  E.g. if `re_coef[[1]]` is a vector and `re_coef[[2]]` a matrix with 4
  rows, `S[[1]]` is a list of length 2 and `S[[2]]` is a matrix, then
  `lambda` needs to be of length 1 \* 2 + 4 = 6.

## Value

returns the penalty value and reports to
[`qreml`](https://janolefi.github.io/LaMa/reference/qreml.md).

## Details

**Caution:** The formatting of `re_coef` needs to match the structure of
the parameter list in your penalised negative log-likelihood function,
i.e. you cannot have two random effect vectors of different names
(different list elements in the parameter list), combine them into a
matrix inside your likelihood and pass the matrix to `penalty`. If these
are seperate random effects, each with its own name, they need to be
passed as a list to `penalty`. Moreover, the ordering of `re_coef` needs
to match the character vector `random` specified in
[`qreml`](https://janolefi.github.io/LaMa/reference/qreml.md).

## See also

[`qreml`](https://janolefi.github.io/LaMa/reference/qreml.md) for the
**qREML** algorithm

## Examples

``` r
# Example with a single random effect
re = rep(0, 5)
S = diag(5)
lambda = 1
penalty(re, S, lambda)
#> [1] 0

# Example with two random effects, 
# where one element contains two random effects of similar structure
re = list(matrix(0, 2, 5), rep(0, 4))
S = list(diag(5), diag(4))
lambda = c(1,1,2) # length = total number of random effects
penalty(re, S, lambda)
#> [1] 0

# Full model-fitting example
# \donttest{
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
  Gamma = tpm_g(Z, cbind(beta0, betaspline)) # transition probabilities
  delta = stationary_p(Gamma, t = 1) # initial distribution
  mu = exp(logmu) # step mean
  sigma = exp(logsigma) # step sd
  # calculating all state-dependent densities
  allprobs = matrix(1, nrow = length(step), ncol = N)
  ind = which(!is.na(step)) # only for non-NA obs.
  for(j in 1:N) allprobs[ind,j] = dgamma2(step[ind],mu[j],sigma[j])
  -forward_g(delta, Gamma[,,tod], allprobs) +
      penalty(betaspline, S, lambda) # this does all the penalisation work
}

# model fitting
mod = qreml(pnll, par, dat, random = "betaspline")
#> Creating AD function
#> Initialising with lambda: 10 10
#> outer 1 - lambda: 3.636 2.858 
#> outer 2 - lambda: 0.786 0.954 
#> outer 3 - lambda: 0.26 0.387 
#> outer 4 - lambda: 0.521 0.267 
#> outer 5 - lambda: 0.468 0.185 
#> outer 6 - lambda: 0.49 0.103 
#> outer 7 - lambda: 0.486 0.079 
#> outer 8 - lambda: 0.486 0.08 
#> outer 9 - lambda: 0.486 0.08 
#> outer 10 - lambda: 0.487 0.08 
#> outer 11 - lambda: 0.487 0.08 
#> Converged
#> Final model fit with lambda: 0.486 0.079
# }
```
