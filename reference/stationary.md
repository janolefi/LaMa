# Compute the stationary distribution of a homogeneous Markov chain

A homogeneous, finite state Markov chain that is irreducible and
aperiodic converges to a unique stationary distribution, here called
\\\delta\\. As it is stationary, this distribution satisfies \$\$\delta
\Gamma = \delta,\$\$ subject to \\\sum\_{j=1}^N \delta_j = 1\\, where
\\\Gamma\\ is the transition probability matrix. This function solves
the linear system of equations above.

## Usage

``` r
stationary(Gamma)
```

## Arguments

- Gamma:

  transition probability matrix of dimension `c(N,N)` or array of such
  matrices of dimension `c(N,N,nTracks)` if the stationary distribution
  should be computed for several matrices at once

## Value

either a single stationary distribution of the Markov chain (vector of
length `N`) or a matrix of stationary distributions of dimension
`c(nTracks,N)` with one stationary distribution in each row

## See also

[`tpm`](https://janolefi.github.io/LaMa/reference/tpm.md) to create a
transition probabilty matrix using the multinomial logistic link
(softmax)

Other stationary distribution functions:
[`stationary_cont()`](https://janolefi.github.io/LaMa/reference/stationary_cont.md),
[`stationary_p()`](https://janolefi.github.io/LaMa/reference/stationary_p.md)

## Examples

``` r
# single matrix
Gamma = tpm(c(rep(-2,3), rep(-3,3)))
delta = stationary(Gamma)
# multiple matrices
Gamma = array(Gamma, dim = c(3,3,10))
Delta = stationary(Gamma)
```
