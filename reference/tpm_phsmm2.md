# Build all transition probability matrices of an periodic-HSMM-approximating HMM

Hidden semi-Markov models (HSMMs) are a flexible extension of HMMs. For
direct numerical maximum likelihood estimation, HSMMs can be represented
as HMMs on an enlarged state space (of size \\M\\) and with structured
transition probabilities. This function computes the transition matrices
of a periodically inhomogeneous HSMMs.

## Usage

``` r
tpm_phsmm2(omega, dm, eps = 1e-10)
```

## Arguments

- omega:

  embedded transition probability matrix

  Either a matrix of dimension `c(nStates, nStates)` for homogeneous
  conditional transition probabilities, or an array of dimension
  `c(nStates, nStates, L)` for inhomogeneous conditional transition
  probabilities.

- dm:

  state dwell-time distributions arranged in a list of length `nStates`

  Each list element needs to be a matrix of dimension `c(L, N_i)`, where
  each row t is the (approximate) probability mass function of state i
  at time t.

- eps:

  rounding value: If an entry of the transition probability matrix is
  smaller, it is rounded to zero.

## Value

array of dimension `c(M, M, L)`, containing the extended-state-space
transition probability matrices of the approximating HMM for each time
point of the cycle.

## Examples

``` r
N = 3
L = 24
# time-varying mean dwell times
Lambda = exp(matrix(rnorm(L*N, 2, 0.5), nrow = L))
sizes = c(25, 25, 25) # approximating chain with 75 states
# state dwell-time distributions
dm = list()
for(i in 1:3){
  dmi = matrix(nrow = L, ncol = sizes[i])
  for(t in 1:L){
    dmi[t,] = dpois(1:sizes[i]-1, Lambda[t,i])
  }
  dm[[i]] = dmi
}

## homogeneous conditional transition probabilities
# diagonal elements are zero, rowsums are one
omega = matrix(c(0,0.5,0.5,0.2,0,0.8,0.7,0.3,0), nrow = N, byrow = TRUE)

# calculating extended-state-space t.p.m.s
Gamma = tpm_phsmm(omega, dm)

## inhomogeneous conditional transition probabilities
# omega can be an array
omega = array(rep(omega,L), dim = c(N,N,L))
omega[1,,4] = c(0, 0.2, 0.8) # small change for inhomogeneity

# calculating extended-state-space t.p.m.s
Gamma = tpm_phsmm(omega, dm)
```
