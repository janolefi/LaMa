# Build the transition probability matrix of an HSMM-approximating HMM

Hidden semi-Markov models (HSMMs) are a flexible extension of HMMs. For
direct numerical maximum likelihood estimation, HSMMs can be represented
as HMMs on an enlarged state space (of size \\M\\) and with structured
transition probabilities. This function computes the transition matrix
of an HSMM.

## Usage

``` r
tpm_hsmm2(omega, dm, eps = 1e-10)
```

## Arguments

- omega:

  embedded transition probability matrix of dimension
  `c(nStates, nStates)`

- dm:

  state dwell-time distributions arranged in a list of length `nStates`.
  Each list element needs to be a vector of length N_i, where N_i is the
  state aggregate size.

- eps:

  rounding value: If an entry of the transition probability matrix is
  smaller, it is rounded to zero.

## Value

extended-state-space transition probability matrix of the approximating
HMM

## Examples

``` r
# building the t.p.m. of the embedded Markov chain
omega = matrix(c(0,1,1,0), nrow = 2, byrow = TRUE)
# defining state aggregate sizes
sizes = c(20, 30)
# defining state dwell-time distributions
lambda = c(5, 11)
dm = list(dpois(1:sizes[1]-1, lambda[1]), dpois(1:sizes[2]-1, lambda[2]))
# calculating extended-state-space t.p.m.
Gamma = tpm_hsmm(omega, dm)
```
