# Build the embedded transition probability matrix of an HSMM from unconstrained parameter vector

Hidden semi-Markov models are defined in terms of state durations and an
**embedded** transition probability matrix that contains the conditional
transition probabilities given that the **current state is left**. This
matrix necessarily has diagonal entries all equal to zero as
self-transitions are impossible.

This function builds such an embedded/ conditional transition
probability matrix from an unconstrained parameter vector. For each row
of the matrix, the inverse multinomial logistic link is applied.

For a matrix of dimension c(N,N), the number of free off-diagonal
elements is N\*(N-2), hence also the length of `param`. This means, for
2 states, the function needs to be called without any arguments, for
3-states with a vector of length 3, for 4 states with a vector of length
8, etc.

Compatible with automatic differentiation by `RTMB`

## Usage

``` r
tpm_emb(param = NULL)
```

## Arguments

- param:

  unconstrained parameter vector of length N\*(N-2) where N is the
  number of states of the Markov chain

  If the function is called without `param`, it will return the
  conditional transition probability matrix for a 2-state HSMM, which is
  fixed with 0 diagonal entries and off-diagonal entries equal to 1.

## Value

embedded/ conditional transition probability matrix of dimension c(N,N)

## See also

Other transition probability matrix functions:
[`generator()`](https://janolefi.github.io/LaMa/reference/generator.md),
[`tpm()`](https://janolefi.github.io/LaMa/reference/tpm.md),
[`tpm_cont()`](https://janolefi.github.io/LaMa/reference/tpm_cont.md),
[`tpm_emb_g()`](https://janolefi.github.io/LaMa/reference/tpm_emb_g.md),
[`tpm_g()`](https://janolefi.github.io/LaMa/reference/tpm_g.md),
[`tpm_g2()`](https://janolefi.github.io/LaMa/reference/tpm_g2.md),
[`tpm_p()`](https://janolefi.github.io/LaMa/reference/tpm_p.md)

## Examples

``` r
# 2 states: no free off-diagonal elements
omega = tpm_emb()

# 3 states: 3 free off-diagonal elements
param = rep(0, 3)
omega = tpm_emb(param)

# 4 states: 8 free off-diagonal elements
param = rep(0, 8)
omega = tpm_emb(param)
```
