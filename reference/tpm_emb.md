# Build the embedded transition probability matrix of an HSMM from unconstrained parameter vector

**Deprecated:** please use
[`tpm_switch`](https://janolefi.github.io/LaMa/reference/tpm_switch.md)
instead, which gives the same matrix.

## Usage

``` r
tpm_emb(param = NULL)
```

## Arguments

- param:

  unconstrained parameter vector of length `nStates * (nStates-2)`, or
  `NULL` for 2 states

## Value

embedded/ conditional transition probability matrix of dimension
`c(nStates, nStates)`

## Examples

``` r
omega = tpm_emb(rep(0, 3))
#> tpm_emb() is deprecated, please use tpm_switch() instead.
```
