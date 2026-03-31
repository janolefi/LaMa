# Compute the transition probability matrix of a thinned periodically inhomogeneous Markov chain.

If the transition probability matrix of an inhomogeneous Markov chain
varies only periodically (with period length \\L\\), it converges to a
so-called periodically stationary distribution. This happens, because
the thinned Markov chain, which has a full cycle as each time step, has
homogeneous transition probability matrix \$\$\Gamma_t = \Gamma^{(t)}
\Gamma^{(t+1)} \dots \Gamma^{(t+L-1)}\$\$ for all \\t = 1, \dots, L.\\
This function calculates the matrix above efficiently as a preliminery
step to calculating the periodically stationary distribution.

## Usage

``` r
tpm_thinned(Gamma, t)
```

## Arguments

- Gamma:

  array of transition probability matrices of dimension
  `c(nStates, nStates, L)`.

- t:

  integer index of the time point in the cycle, for which to calculate
  the thinned transition probility matrix

## Value

thinned transition probabilty matrix of dimension `c(nStates, nStates)`

## Examples

``` r
# setting parameters for trigonometric link
beta = matrix(c(-1, -2, 2, -1, 2, -4), nrow = 2, byrow = TRUE)
# calculating periodically varying t.p.m. array (of length 24 here)
Gamma = tpm_p(beta = beta)
# calculating t.p.m. of thinned Markov chain
tpm_thinned(Gamma, 4)
#>           [,1]      [,2]
#> [1,] 0.8926642 0.1073358
#> [2,] 0.8926642 0.1073358
```
