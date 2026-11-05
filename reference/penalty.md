# Computes penalty based on quadratic form

This function computes quadratic penalties of the form \$\$0.5 \sum\_{i}
\lambda_i b_i^T S_i b_i,\$\$ with smoothing parameters \\\lambda_i\\,
coefficient vectors \\b_i\\, and fixed penalty matrices \\S_i\\.

It is intended to be used inside the **penalised negative log-likelihood
function** when fitting models with penalised splines or simple random
effects via **quasi restricted maximum likelihood** (qREML) with the
[`qreml`](https://janolefi.github.io/LaMa/reference/qreml.md) function.
For [`qreml`](https://janolefi.github.io/LaMa/reference/qreml.md) to
work, the likelihood function needs to be compatible with the `RTMB` R
package to enable automatic differentiation.

## Usage

``` r
penalty(re_coef, S, sp, lambda)
```

## Arguments

- re_coef:

  coefficient vector/ matrix or list of coefficient vectors/ matrices

  Each list entry corresponds to a different smooth/ random effect with
  its own associated penalty matrix in `S`. When several smooths/ random
  effects of the same kind are present, it is convenient to pass them as
  a matrix, where each row corresponds to one smooth/ random effect.
  This way all rows can use the same penalty matrix.

- S:

  fixed penalty matrix or list of penalty matrices matching the
  structure of `re_coef` and also the dimension of the individuals
  smooths/ random effects

- sp:

  penalty strength parameter vector that has a length corresponding to
  the **total number** of random effects/ spline coefficients in
  `re_coef`

  E.g. if `re_coef` contains one vector and one matrix with 4 rows, then
  `sp` needs to be of length 5.

- lambda:

  former name of `sp`, still works

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

## References

Koslik, J. O. (2024). Efficient smoothness selection for nonparametric
Markov-switching models via quasi restricted maximum likelihood. arXiv
preprint arXiv:2411.11498.

## See also

[`qreml`](https://janolefi.github.io/LaMa/reference/qreml.md) for the
**qREML** algorithm

## Examples

``` r
# Example with a single random effect
re = rep(0, 5)
S = diag(5)
sp = 1
penalty(re, S, sp)
#> [1] 0

# Example with two random effects, 
# where one element contains two random effects of similar structure
re = list(matrix(0, 2, 5), rep(0, 4))
S = list(diag(5), diag(4))
sp = c(1,1,2) # length = total number of random effects
penalty(re, S, sp)
#> [1] 0

# Full model-fitting example
# \donttest{
data = trex[1:1000,] # subset

# initial parameter list
par = list(logmu = log(c(0.3, 2.5)), # step mean
           logsigma = log(c(0.3, 1.5)), # step sd
           beta0 = c(-2,-2), # state process intercept
           betaspline = matrix(rep(0, 18), nrow = 2)) # state process spline coefs
          
# data object with initial penalty strength sp
dat = list(step = data$step, # step length
           tod = data$tod, # time of day covariate
           N = 2, # number of states
           sp = rep(10,2)) # initial penalty strength

# building model matrices
modmat = modmat(~ s(tod, bs = "cp"), 
                       data = data.frame(tod = 1:24), 
                       knots = list(tod = c(0,24))) # wrapping points
dat$Z = modmat$Z # spline design matrix
dat$S = modmat$S # penalty matrix

# penalised negative log-likelihood function
pnll = function(par) {
  getAll(par, dat) # makes everything contained available without $
  Gamma = tpm_g(Z, cbind(beta0, betaspline)) # transition probabilities
  delta = stationary_p(Gamma, t = tod[1] %% 24 + 1) # initial distribution
  mu = exp(logmu) # step mean
  sigma = exp(logsigma) # step sd
  # calculating all state-dependent densities
  allprobs = matrix(1, nrow = length(step), ncol = N)
  ind = which(!is.na(step)) # only for non-NA obs.
  for(j in 1:N) allprobs[ind,j] = dgamma2(step[ind],mu[j],sigma[j])
  -forward_g(delta, Gamma[,,tod], allprobs) +
      penalty(betaspline, S, sp) # this does all the penalization work
}

# model fitting
mod = qreml(pnll, par, dat, random = "betaspline")
#> Creating AD function
#> Initialising with sp: 10 10
#> outer 1 - sp: 3.588 2.862 
#> outer 2 - sp: 1.651 1.653 
#> outer 3 - sp: 0.941 1.186 
#> outer 4 - sp: 0.656 0.921 
#> outer 5 - sp: 0.537 0.742 
#> outer 6 - sp: 0.488 0.606 
#> outer 7 - sp: 0.468 0.499 
#> outer 8 - sp: 0.461 0.411 
#> outer 9 - sp: 0.46 0.339 
#> outer 10 - sp: 0.462 0.28 
#> outer 11 - sp: 0.465 0.234 
#> outer 12 - sp: 0.469 0.197 
#> outer 13 - sp: 0.472 0.169 
#> outer 14 - sp: 0.476 0.148 
#> outer 15 - sp: 0.478 0.133 
#> outer 16 - sp: 0.481 0.121 
#> outer 17 - sp: 0.483 0.113 
#> outer 18 - sp: 0.485 0.107 
#> outer 19 - sp: 0.486 0.103 
#> outer 20 - sp: 0.487 0.1 
#> outer 21 - sp: 0.488 0.097 
#> outer 22 - sp: 0.488 0.096 
#> outer 23 - sp: 0.488 0.095 
#> outer 24 - sp: 0.489 0.094 
#> outer 25 - sp: 0.489 0.093 
#> outer 26 - sp: 0.489 0.093 
#> Converged
#> Final sp: 0.489 0.093
# }
```
