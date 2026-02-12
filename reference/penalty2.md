# Computes generalised quadratic-form penalties

This function computes a quadratic penalty of the form \$\$0.5 \sum\_{i}
\lambda_i b^T S_i b,\$\$ with smoothing parameters \\\lambda_i\\,
coefficient vector \\b\\, and fixed penalty matrices \\S_i\\. This
generalises the
[`penalty`](https://janoleko.github.io/reference/penalty.md) by allowing
subsets of the coefficient vector \\b\\ to be penalised multiple times
with different smoothing parameters, which is necessary for **tensor
products**, **functional random effects** or **adaptive smoothing**.

It is intended to be used inside the **penalised negative log-likelihood
function** when fitting models with penalised splines or simple random
effects via **quasi restricted maximum likelihood** (qREML) with the
[`qreml`](https://janoleko.github.io/reference/qreml.md) function. For
[`qreml`](https://janoleko.github.io/reference/qreml.md) to work, the
likelihood function needs to be compatible with the `RTMB` R package to
enable automatic differentiation.

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
[`qreml`](https://janoleko.github.io/reference/qreml.md).

## Details

**Caution:** The formatting of `re_coef` needs to match the structure of
the parameter list in your penalised negative log-likelihood function,
i.e. you cannot have two random effect vectors of different names
(different list elements in the parameter list), combine them into a
matrix inside your likelihood and pass the matrix to `penalty`. If these
are seperate random effects, each with its own name, they need to be
passed as a list to `penalty`. Moreover, the ordering of `re_coef` needs
to match the character vector `random` specified in
[`qreml`](https://janoleko.github.io/reference/qreml.md).

## See also

[`qreml`](https://janoleko.github.io/reference/qreml.md) for the
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
data = trex[1:1000,] # subset

# initial parameter list
par = list(logmu = log(c(0.3, 1)), # step mean
           logsigma = log(c(0.2, 0.7)), # step sd
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
      penalty(betaspline, S, lambda) # this does all the penalization work
}

# model fitting
mod = qreml(pnll, par, dat, random = "betaspline")
#> Creating AD function
#> Performance tip: Consider running `TapeConfig(matmul = 'plain')` before `MakeADFun()` to speed up the forward algorithm.
#> Initialising with lambda: 10 10 
#> outer 1 - lambda: 5.545 5.001 
#> outer 2 - lambda: 3.289 3.001 
#> outer 3 - lambda: 2.081 2.091 
#> outer 4 - lambda: 1.406 1.599 
#> outer 5 - lambda: 1.018 1.293 
#> outer 6 - lambda: 0.79 1.08 
#> outer 7 - lambda: 0.654 0.92 
#> outer 8 - lambda: 0.573 0.794 
#> outer 9 - lambda: 0.525 0.69 
#> outer 10 - lambda: 0.497 0.602 
#> outer 11 - lambda: 0.48 0.526 
#> outer 12 - lambda: 0.471 0.46 
#> outer 13 - lambda: 0.467 0.403 
#> outer 14 - lambda: 0.465 0.353 
#> outer 15 - lambda: 0.465 0.31 
#> outer 16 - lambda: 0.466 0.272 
#> outer 17 - lambda: 0.467 0.24 
#> outer 18 - lambda: 0.469 0.214 
#> outer 19 - lambda: 0.472 0.191 
#> outer 20 - lambda: 0.474 0.172 
#> outer 21 - lambda: 0.476 0.157 
#> outer 22 - lambda: 0.478 0.144 
#> outer 23 - lambda: 0.48 0.134 
#> outer 24 - lambda: 0.481 0.126 
#> outer 25 - lambda: 0.483 0.119 
#> outer 26 - lambda: 0.484 0.114 
#> outer 27 - lambda: 0.485 0.109 
#> outer 28 - lambda: 0.486 0.106 
#> outer 29 - lambda: 0.487 0.103 
#> outer 30 - lambda: 0.488 0.101 
#> outer 31 - lambda: 0.488 0.099 
#> outer 32 - lambda: 0.489 0.097 
#> outer 33 - lambda: 0.489 0.096 
#> outer 34 - lambda: 0.489 0.095 
#> outer 35 - lambda: 0.49 0.094 
#> outer 36 - lambda: 0.49 0.094 
#> outer 37 - lambda: 0.49 0.093 
#> outer 38 - lambda: 0.49 0.093 
#> Converged
#> Final model fit with lambda: 0.49 0.093 
```
