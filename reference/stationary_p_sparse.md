# Sparse version of [`stationary_p`](https://janoleko.github.io/reference/stationary_p.md)

This is function computes the periodically stationary distribution of a
Markov chain given a list of L **sparse** transition probability
matrices. Compatible with automatic differentiation by `RTMB`

## Usage

``` r
stationary_p_sparse(Gamma, t = NULL)
```

## Arguments

- Gamma:

  sist of length L containing sparse transition probability matrices for
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
## periodic HSMM example (here the approximating tpm is sparse)
N = 2 # number of states
L = 24 # cycle length
# time-varying mean dwell times
Z = cosinor(1:L, period = L) # trigonometric basis functions design matrix
beta = matrix(c(2, 2, 0.1, -0.1, -0.2, 0.2), nrow = 2)
Lambda = exp(cbind(1, Z) %*% t(beta))
sizes = c(20, 20) # approximating chain with 40 states
# state dwell-time distributions
dm = lapply(1:N, function(i) sapply(1:sizes[i]-1, dpois, lambda = Lambda[,i]))
omega = matrix(c(0,1,1,0), nrow = N, byrow = TRUE) # embedded t.p.m.

# calculating extended-state-space t.p.m.s
Gamma = tpm_phsmm(omega, dm)
# Periodically stationary distribution for specific time point
delta = stationary_p_sparse(Gamma, 4)

# All periodically stationary distributions
Delta = stationary_p_sparse(Gamma)
```
