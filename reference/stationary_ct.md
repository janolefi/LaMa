# Compute the stationary distribution of a continuous-time Markov chain

A well-behaved continuous-time Markov chain converges to a unique
stationary distribution, here called \\\pi\\. This distribution
satisfies \$\$\pi Q = 0,\$\$ subject to \\\sum\_{j=1}^N \pi_j = 1\\,
where \\Q\\ is the infinitesimal generator of the Markov chain. This
function solves the linear system of equations above for a given
generator matrix.

## Usage

``` r
stationary_ct(Q)

stationary_cont(Q)
```

## Arguments

- Q:

  infinitesimal generator matrix of dimension `c(N,N)` or array of such
  matrices of dimension `c(N,N,nTracks)` if the stationary distribution
  should be computed for several matrices at once

## Value

either a single stationary distribution of the continuous-time Markov
chain (vector of length `N`) or a matrix of stationary distributions of
dimension `c(nTracks,N)` with one stationary distribution in each row

## See also

[`generator`](https://janolefi.github.io/LaMa/reference/generator.md) to
create a generator matrix

Other stationary distribution functions:
[`stationary()`](https://janolefi.github.io/LaMa/reference/stationary.md),
[`stationary_p()`](https://janolefi.github.io/LaMa/reference/stationary_p.md)

## Examples

``` r
# single matrix
Q = generator(c(-2,-2))
Pi = stationary_ct(Q)
# multiple matrices
Q = array(Q, dim = c(2,2,10))
Pi = stationary_ct(Q)
```
