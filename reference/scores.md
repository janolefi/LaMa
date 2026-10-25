# Per-observation score contributions

Computes the score contribution of each observation, \$\$u_t =
\frac{\partial}{\partial \theta} \log p(x_t \mid x_1, \dots,
x\_{t-1}),\$\$ for a model whose negative log-likelihood function calls
[`forward`](https://janolefi.github.io/LaMa/reference/forward.md). Their
column sums are the gradient of the log-likelihood. They are the
building block for robust (sandwich) standard errors, e.g.
cluster-robust ones via `rowsum(U, ID)`.

## Usage

``` r
scores(obj, par = NULL)
```

## Arguments

- obj:

  optimised `RTMB` object from `MakeADFun`, or a model object from
  [`lama`](https://janolefi.github.io/LaMa/reference/lama.md) or
  [`qreml`](https://janolefi.github.io/LaMa/reference/qreml.md)

- par:

  optional parameter vector like `obj$par`, e.g. `opt$par`. Defaults to
  the estimates (for
  [`qreml`](https://janolefi.github.io/LaMa/reference/qreml.md) models
  the final fit).

## Value

matrix of score contributions of dimension `c(nObs, length(obj$par))`

## Details

Each log-likelihood contribution in
[`forward`](https://janolefi.github.io/LaMa/reference/forward.md) gets a
weight \\w_t\\. The gradient of the weighted log-likelihood is linear in
the weights, so its derivative w.r.t. \\w_t\\ is \\u_t\\. This needs one
tape of the gradient and one reverse sweep per parameter. The weights
are handed to
[`forward`](https://janolefi.github.io/LaMa/reference/forward.md)
internally, so the likelihood function does not need to change.

The likelihood function is taken from `obj` and evaluated again, with
the data as they are when `scores()` is called. Parameters fixed via
`map` are respected, and the columns correspond to `obj$par`, like
`obj$gr()` and `obj$he()`.

If the likelihood function calls
[`forward`](https://janolefi.github.io/LaMa/reference/forward.md)
several times, e.g. once per track, the rows of the result follow the
order of these calls. Only
[`forward`](https://janolefi.github.io/LaMa/reference/forward.md) (and
its deprecated `_g` and `_p` versions) is supported, and models with
random effects are not.

## Examples

``` r
nll <- function(par) {
  getAll(par)
  Gamma <- tpm(eta)
  delta <- stationary(Gamma)
  allprobs <- cbind(dgamma2(step, mu[1], sigma[1]), 
                    dgamma2(step, mu[2], sigma[2]))
  -forward(delta, Gamma, allprobs)
}
step <- trex$step[1:200]
par <- list(eta = c(-2, -2), mu = c(0.3, 2.5), sigma = c(0.2, 1.5))
obj <- MakeADFun(nll, par, silent = TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)
mod <- lama(obj)
U <- scores(mod)
colSums(U) # ~ 0 at the MLE
#>           eta           eta            mu            mu         sigma 
#> -3.121369e-06  2.549282e-07 -4.984111e-05  1.896265e-05 -3.926205e-05 
#>         sigma 
#> -3.106319e-05 
```
