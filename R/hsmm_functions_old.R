#' Build the transition probability matrix of an HSMM-approximating HMM
#'
#' @description
#' Hidden semi-Markov models (HSMMs) are a flexible extension of HMMs. 
#' For direct numerical maximum likelhood estimation, HSMMs can be represented as HMMs on an enlarged state space (of size \eqn{M}) and with structured transition probabilities.
#' This function computes the transition matrix of an HSMM.
#'
#' @param omega embedded transition probability matrix of dimension c(N,N)
#' @param dm state dwell-time distributions arranged in a list of length(N). Each list element needs to be a vector of length N_i, where N_i is the state aggregate size.
#' @param eps rounding value: If an entry of the transition probabily matrix is smaller, than it is rounded to zero.
#'
#' @return extended-state-space transition probability matrix of the approximating HMM
#' @export
#'
#' @examples
#' # building the t.p.m. of the embedded Markov chain
#' omega = matrix(c(0,1,1,0), nrow = 2, byrow = TRUE)
#' # defining state aggregate sizes
#' sizes = c(20, 30)
#' # defining state dwell-time distributions
#' lambda = c(5, 11)
#' dm = list(dpois(1:sizes[1]-1, lambda[1]), dpois(1:sizes[2]-1, lambda[2]))
#' # calculating extended-state-space t.p.m.
#' Gamma = tpm_hsmm(omega, dm)
tpm_hsmm2 = function(omega,dm,eps=1e-10){
  mv = vapply(dm, length, integer(1))
  m = length(mv)
  G = matrix(0,0,sum(mv))
  for (i in 1:m){
    mi = mv[[i]]
    F = cumsum(c(0,dm[[i]][-mi]))
    ci = ifelse(abs(1-F)>eps,dm[[i]]/(1-F),1)
    cim = ifelse(1-ci>0,1-ci,0)
    Gi = matrix(0,mi,0)
    for (j in 1:m){
      if(i==j) {
        if(mi==1){ 
          Gi = cbind(Gi,c(rep(0,mv[[j]]-1),cim))
        } else{ 
          Gi = cbind(Gi,rbind(cbind(rep(0,mi-1),diag(cim[-mi],mi-1,mi-1)),
                              c(rep(0,mi-1),cim[[mi]])))}
      } else   { if(mi==1)
      { Gi = cbind(Gi,matrix(c(omega[[i,j]]*ci,rep(0,mv[[j]]-1)),1))} else
      { Gi = cbind(Gi,cbind(omega[[i,j]]*ci,matrix(0,mv[[i]],mv[[j]]-1)))}
      }
    }
    G = rbind(G,Gi)
  }
  G 
}


#' Build all transition probability matrices of an periodic-HSMM-approximating HMM
#'
#' Hidden semi-Markov models (HSMMs) are a flexible extension of HMMs. For direct numerical maximum likelhood estimation, HSMMs can be represented as HMMs on an enlarged state space (of size \eqn{M}) and with structured transition probabilities.
#' This function computes the transition matrices of a periodically inhomogeneos HSMMs.
#'
#' @param omega embedded transition probability matrix
#' 
#' Either a matrix of dimension c(N,N) for homogeneous conditional transition probabilities, or an array of dimension c(N,N,L) for inhomogeneous conditional transition probabilities.
#' @param dm state dwell-time distributions arranged in a list of length(N)
#' 
#' Each list element needs to be a matrix of dimension c(L, N_i), where each row t is the (approximate) probability mass function of state i at time t.
#' @param eps rounding value: If an entry of the transition probabily matrix is smaller, than it is rounded to zero.
#'
#' @return array of dimension c(N,N,L), containing the extended-state-space transition probability matrices of the approximating HMM for each time point of the cycle.
#' @export
#'
#' @examples
#' N = 3
#' L = 24
#' # time-varying mean dwell times
#' Lambda = exp(matrix(rnorm(L*N, 2, 0.5), nrow = L))
#' sizes = c(25, 25, 25) # approximating chain with 75 states
#' # state dwell-time distributions
#' dm = list()
#' for(i in 1:3){
#'   dmi = matrix(nrow = L, ncol = sizes[i])
#'   for(t in 1:L){
#'     dmi[t,] = dpois(1:sizes[i]-1, Lambda[t,i])
#'   }
#'   dm[[i]] = dmi
#' }
#' 
#' ## homogeneous conditional transition probabilites
#' # diagonal elements are zero, rowsums are one
#' omega = matrix(c(0,0.5,0.5,0.2,0,0.8,0.7,0.3,0), nrow = N, byrow = TRUE)
#' 
#' # calculating extended-state-space t.p.m.s
#' Gamma = tpm_phsmm(omega, dm)
#' 
#' ## inhomogeneous conditional transition probabilites
#' # omega can be an array
#' omega = array(rep(omega,L), dim = c(N,N,L))
#' omega[1,,4] = c(0, 0.2, 0.8) # small change for inhomogeneity
#' 
#' # calculating extended-state-space t.p.m.s
#' Gamma = tpm_phsmm(omega, dm)
tpm_phsmm2 = function(omega, dm, eps = 1e-10){
  # dm list over states: entries matrices of dim c(L, N_i)
  L = nrow(dm[[1]]) # length of one cycle
  N = length(dm) # number of states
  # allowing for inhomogeneous cond. transition probs
  if(is.matrix(omega)){
    omega = array(rep(omega,L), dim = c(N,N,L))
  }
  mv = vapply(dm, ncol, integer(1)) # lengths of the pmf vectors
  M = sum(mv)
  G_all = array(dim = c(M,M,L))
  Fm = vector("list") # computing all cdfs
  for(i in 1:N){
    Fm[[i]] = cbind(0,t(apply(dm[[i]][,-mv[i]], 1, cumsum)))
  }
  for (t in 1:L){ # loop over time
    G = matrix(0,0,sum(mv)) # constructing an empty matrix with number of columns equal to final dimension (sum_i m_i)
    for (i in 1:N){ # loop over states
      Ni = mv[i]
      ci = numeric(Ni) # length of pmf for state aggregate i
      for (k in 0:(Ni-1)){
        l = ifelse(t==k, L, (t-k)%%L)
        l = ifelse(l==0, L, l)
        ci[k+1] = ifelse(abs(1-Fm[[i]][l,k+1]) > eps,
                         dm[[i]][l,k+1]/(1-Fm[[i]][l,k+1]), 1)
      }
      cim = ifelse(1-ci > 0, 1-ci, 0)
      Gi = matrix(0, Ni, 0)
      for (j in 1:N){
        if(i==j) {
          if(Ni==1){ 
            Gi = cbind(Gi, c(rep(0,mv[[j]]-1), cim))
          } else{ 
            Gi = cbind(Gi, rbind(cbind(rep(0,Ni-1),diag(cim[-Ni],Ni-1,Ni-1)),
                                 c(rep(0,Ni-1),cim[[Ni]])))}
        } else { if(Ni==1)
        { Gi = cbind(Gi, matrix(c(omega[[i,j,t]]*ci, rep(0,mv[[j]]-1)),1))} else
        { Gi = cbind(Gi, cbind(omega[[i,j,t]]*ci, matrix(0,mv[[i]],mv[[j]]-1)))}
        }
      }
      G = rbind(G,Gi)
    }
    G_all[,,t] = G
  }
  G_all
}


#' \href{https://www.taylorfrancis.com/books/mono/10.1201/b20790/hidden-markov-models-time-series-walter-zucchini-iain-macdonald-roland-langrock}{Forward algorithm} for hidden semi-Markov models with homogeneous transition probability matrix
#'
#' Hidden semi-Markov models (HSMMs) are a flexible extension of HMMs that can be approximated by HMMs on an enlarged state space (of size \eqn{M}) and with structured transition probabilities.
#'
#' @param delta initial or stationary distribution of length N, or matrix of dimension c(k,N) for k independent tracks, if \code{trackID} is provided
#' @param Gamma transition probability matrix of dimension c(M,M)
#' @param allprobs matrix of state-dependent probabilities/ density values of dimension c(n, N) which will automatically be converted to the appropriate dimension.
#' @param sizes state aggregate sizes that are used for the approximation of the semi-Markov chain.
#'
#' @return log-likelihood for given data and parameters
#' @export
#'
#' @examples
#' ## generating data from homogeneous 2-state HSMM
#' mu = c(0, 6)
#' lambda = c(6, 12)
#' omega = matrix(c(0,1,1,0), nrow = 2, byrow = TRUE)
#' # simulation
#' # for a 2-state HSMM the embedded chain always alternates between 1 and 2
#' s = rep(1:2, 100)
#' C = x = numeric(0)
#' for(t in 1:100){
#'   dt = rpois(1, lambda[s[t]])+1 # shifted Poisson
#'   C = c(C, rep(s[t], dt))
#'   x = c(x, rnorm(dt, mu[s[t]], 1.5)) # fixed sd 2 for both states
#' }
#' 
#' ## negative log likelihood function
#' mllk = function(theta.star, x, sizes){
#'   # parameter transformations for unconstraint optimization
#'   omega = matrix(c(0,1,1,0), nrow = 2, byrow = TRUE) # omega fixed (2-states)
#'   lambda = exp(theta.star[1:2]) # dwell time means
#'   dm = list(dpois(1:sizes[1]-1, lambda[1]), dpois(1:sizes[2]-1, lambda[2]))
#'   Gamma = tpm_hsmm2(omega, dm)
#'   delta = stationary(Gamma) # stationary
#'   mu = theta.star[3:4]
#'   sigma = exp(theta.star[5:6])
#'   # calculate all state-dependent probabilities
#'   allprobs = matrix(1, length(x), 2)
#'   for(j in 1:2){ allprobs[,j] = dnorm(x, mu[j], sigma[j]) }
#'   # return negative for minimization
#'   -forward_s(delta, Gamma, allprobs, sizes)
#' }
#' 
#' ## fitting an HSMM to the data
#' theta.star = c(log(5), log(10), 1, 4, log(2), log(2))
#' mod = nlm(mllk, theta.star, x = x, sizes = c(20, 30), stepmax = 5)
forward_s = function(delta, Gamma, allprobs, sizes){
  forward_cpp_s(allprobs, delta, Gamma, sizes)
}


#' \href{https://www.taylorfrancis.com/books/mono/10.1201/b20790/hidden-markov-models-time-series-walter-zucchini-iain-macdonald-roland-langrock}{Forward algorithm} for hidden semi-Markov models with periodically varying transition probability matrices
#'
#' Hidden semi-Markov models (HSMMs) are a flexible extension of HMMs that can be approximated by HMMs on an enlarged state space (of size \eqn{M}) and with structured transition probabilities.
#' Recently, this inference procedure has been generalised to allow either the dwell-time distributions or the conditional transition probabilities to depend on external covariates such as the time of day. This special case is implemented here.
#' This function allows for that, by expecting a transition probability matrix for each time point in a period, and an integer valued (\eqn{1, \dots, L}) time variable that maps the data index to the according time.
#'
#' @param delta initial or stationary distribution of length N, or matrix of dimension c(k,N) for k independent tracks, if \code{trackID} is provided
#' @param Gamma array of transition probability matrices of dimension c(M,M,L).
#' 
#' Here we use the definition \eqn{\Pr(S_t=j \mid S_{t-1}=i) = \gamma_{ij}^{(t)}} such that the transition probabilities between time point \eqn{t-1} and \eqn{t} are an element of \eqn{\Gamma^{(t)}}.
#' @param allprobs matrix of state-dependent probabilities/ density values of dimension c(n, N) which will automatically be converted to the appropriate dimension.
#' @param sizes state aggregate sizes that are used for the approximation of the semi-Markov chain.
#' @param tod (Integer valued) variable for cycle indexing in 1, ..., L, mapping the data index to a generalised time of day (length n).
#' For half-hourly data L = 48. It could, however, also be day of year for daily data and L = 365.
#'
#' @return log-likelihood for given data and parameters
#' @export
#'
#' @examples
#' ## generating data from homogeneous 2-state HSMM
#' mu = c(0, 6)
#' beta = matrix(c(log(4),log(6),-0.2,0.2,-0.1,0.4), nrow=2)
#' # time varying mean dwell time
#' Lambda = exp(cbind(1, cosinor(1:24))%*%t(beta))
#' omega = matrix(c(0,1,1,0), nrow = 2, byrow = TRUE)
#' # simulation
#' # for a 2-state HSMM the embedded chain always alternates between 1 and 2
#' s = rep(1:2, 100)
#' C = x = numeric(0)
#' tod = rep(1:24, 50) # time of day variable
#' time = 1
#' for(t in 1:100){
#'   dt = rpois(1, Lambda[tod[time], s[t]])+1 # dwell time depending on time of day
#'   time = time + dt
#'   C = c(C, rep(s[t], dt))
#'   x = c(x, rnorm(dt, mu[s[t]], 1.5)) # fixed sd 2 for both states
#' }
#' tod = tod[1:length(x)]
#' 
#' ## negative log likelihood function
#' mllk = function(theta.star, x, sizes, tod){
#'   # parameter transformations for unconstraint optimization
#'   omega = matrix(c(0,1,1,0), nrow = 2, byrow = TRUE) # omega fixed (2-states)
#'   mu = theta.star[1:2]
#'   sigma = exp(theta.star[3:4])
#'   beta = matrix(theta.star[5:10], nrow=2)
#'   # time varying mean dwell time
#'   Lambda = exp(cbind(1, cosinor(1:24))%*%t(beta))
#'   dm = list()
#'   for(j in 1:2){
#'     dm[[j]] = sapply(1:sizes[j]-1, dpois, lambda = Lambda[,j])
#'   }
#'   Gamma = tpm_phsmm2(omega, dm)
#'   delta = stationary_p(Gamma, tod[1])
#'   # calculate all state-dependent probabilities
#'   allprobs = matrix(1, length(x), 2)
#'   for(j in 1:2){ allprobs[,j] = dnorm(x, mu[j], sigma[j]) }
#'   # return negative for minimization
#'   -forward_sp(delta, Gamma, allprobs, sizes, tod)
#' }
#' 
#' ## fitting an HSMM to the data
#' theta.star = c(1, 4, log(2), log(2), # state-dependent parameters
#'                  log(4), log(6), rep(0,4)) # state process parameters dm
#' # mod = nlm(mllk, theta.star, x = x, sizes = c(10, 15), tod = tod, stepmax = 5)
forward_sp = function(delta, Gamma, allprobs, sizes, tod){
  if(min(tod)==1){
    tod = tod-1
  } 
  forward_cpp_sp(allprobs, delta, Gamma, sizes, tod)
}