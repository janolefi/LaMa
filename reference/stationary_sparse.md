# Sparse version of [`stationary`](https://janolefi.github.io/LaMa/reference/stationary.md)

**Deprecated:** please use
[`stationary`](https://janolefi.github.io/LaMa/reference/stationary.md)
instead, which gives the same result.

## Usage

``` r
stationary_sparse(Gamma)
```

## Arguments

- Gamma:

  sparse transition probability matrix of dimension
  `c(nStates, nStates)`

## Value

stationary distribution of the Markov chain with the given transition
probability matrix

## Examples

``` r
Gamma = tpm_hsmm(tpm_switch(), list(dpois(0:19, 5), dpois(0:29, 11)))
delta = stationary_sparse(Gamma)
#> stationary_sparse() is deprecated, please use stationary() instead.
```
