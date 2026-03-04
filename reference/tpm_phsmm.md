# Builds all transition probability matrices of an periodic-HSMM-approximating HMM

Hidden semi-Markov models (HSMMs) are a flexible extension of HMMs. For
direct numerical maximum likelhood estimation, HSMMs can be represented
as HMMs on an enlarged state space (of size \\M\\) and with structured
transition probabilities.

This function computes the transition matrices of a periodically
inhomogeneos HSMMs.

## Usage

``` r
tpm_phsmm(omega, dm, eps = 1e-10)
```

## Arguments

- omega:

  embedded transition probability matrix

  Either a matrix of dimension c(N,N) for homogeneous conditional
  transition probabilities (as computed by
  [`tpm_emb`](https://janoleko.github.io/reference/tpm_emb.md)), or an
  array of dimension c(N,N,L) for inhomogeneous conditional transition
  probabilities (as computed by
  [`tpm_emb_g`](https://janoleko.github.io/reference/tpm_emb_g.md)).

- dm:

  state dwell-time distributions arranged in a list of length N

  Each list element needs to be a matrix of dimension c(L, N_i), where
  each row t is the (approximate) probability mass function of state i
  at time t.

- eps:

  rounding value: If an entry of the transition probabily matrix is
  smaller, than it is rounded to zero. Usually, this should not be
  changed.

## Value

list of dimension length L, containing sparse extended-state-space
transition probability matrices of the approximating HMM for each time
point of the cycle.

## Examples

``` r
N = 2 # number of states
L = 24 # cycle length
# time-varying mean dwell times
Z = cosinor(1:L, period = L) # trigonometric basis functions design matrix
beta = matrix(c(2, 2, 0.1, -0.1, -0.2, 0.2), nrow = 2)
Lambda = exp(cbind(1, Z) %*% t(beta))
sizes = c(20, 20) # approximating chain with 40 states
# state dwell-time distributions
dm = lapply(1:N, function(i) sapply(1:sizes[i]-1, dpois, lambda = Lambda[,i]))

## homogeneous conditional transition probabilites
# diagonal elements are zero, rowsums are one
omega = matrix(c(0,1,1,0), nrow = N, byrow = TRUE)

# calculating extended-state-space t.p.m.s
Gamma = tpm_phsmm(omega, dm)

## inhomogeneous conditional transition probabilites
# omega can be an array
omega = array(omega, dim = c(N,N,L))

# calculating extended-state-space t.p.m.s
Gamma = tpm_phsmm(omega, dm)
```
