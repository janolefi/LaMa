# Sparse version of [`stationary`](https://janolefi.github.io/reference/stationary.md)

This is function computes the stationary distribution of a Markov chain
with a given **sparse** transition probability matrix. Compatible with
automatic differentiation by `RTMB`

## Usage

``` r
stationary_sparse(Gamma)
```

## Arguments

- Gamma:

  sparse transition probability matrix of dimension c(N,N)

## Value

stationary distribution of the Markov chain with the given transition
probability matrix

## Examples

``` r
## HSMM example (here the approximating tpm is sparse)
# building the t.p.m. of the embedded Markov chain
omega = matrix(c(0,1,1,0), nrow = 2, byrow = TRUE)
# defining state aggregate sizes
sizes = c(20, 30)
# defining state dwell-time distributions
lambda = c(5, 11)
dm = list(dpois(1:sizes[1]-1, lambda[1]), dpois(1:sizes[2]-1, lambda[2]))
# calculating extended-state-space t.p.m.
Gamma = tpm_hsmm(omega, dm)
delta = stationary_sparse(Gamma)
```
