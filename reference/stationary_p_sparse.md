# Sparse version of [`stationary_p`](https://janolefi.github.io/LaMa/reference/stationary_p.md)

**Deprecated:** please use
[`stationary_p`](https://janolefi.github.io/LaMa/reference/stationary_p.md)
with an array of transition probability matrices instead, which gives
the same result.

## Usage

``` r
stationary_p_sparse(Gamma, t = NULL)
```

## Arguments

- Gamma:

  list of length L containing sparse transition probability matrices for
  one cycle.

- t:

  integer index of the time point in the cycle, for which to calculate
  the stationary distribution If t is not provided, the function
  calculates all stationary distributions for each time point in the
  cycle.

## Value

either the periodically stationary distribution at time t or all
periodically stationary distributions.

## Examples

``` r
Lambda = exp(cbind(2, 0.3 * sin(2 * pi * 1:24 / 24)))
dm = lapply(1:2, function(i) sapply(0:19, dpois, lambda = Lambda[, i]))
Gamma = tpm_phsmm(tpm_switch(), dm)
Delta = stationary_p_sparse(Gamma)
#> stationary_p_sparse() is deprecated, please use stationary_p() with an array of transition probability matrices instead.
```
