// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "my_functions.h"

using namespace Rcpp;
using namespace arma;

// general forward algorithm, used by forward(), forward_g() and forward_p()
// Gamma.slice(gidx[t]) is the tpm for the transition into t
// trackStart: first index of each track (all indices 0-based)
// [[Rcpp::export]]
double forward_cpp(const arma::mat& allprobs, const arma::mat& Delta, const arma::cube& Gamma,
                   const IntegerVector& gidx, const IntegerVector& trackStart)
{
  int nObs = allprobs.n_rows;
  int K = trackStart.size();

  double l = 0;
  arma::rowvec foo;
  arma::rowvec phi;

  int k = 0; // track index

  for (int t = 0; t < nObs; t++){

    if(k < K && t == trackStart[k]) {
      // start of track k
      foo = Delta.row(k) % allprobs.row(t);
      k++;
    } else {
      foo = (phi * Gamma.slice(gidx[t])) % allprobs.row(t);
    }

    double sumfoo = sum(foo);
    l = l + log(sumfoo);
    phi = foo / sumfoo;
  }

  return l;
}

// [[Rcpp::export]]
double forward_cpp_s(const arma::mat& allprobs, const arma::rowvec& delta, const arma::mat& Gamma, const IntegerVector& agsizes)
{
  int N = allprobs.n_cols;
  int nObs = allprobs.n_rows;
  arma::rowvec foo(N);
  
  foo = delta % rep_times(allprobs.row(0), agsizes);
  double sumfoo = sum(foo);
  double l = log(sumfoo);
  arma::rowvec phi = foo/sumfoo;
  
  for (unsigned int i=1; i<nObs; i++){
    foo = (phi*Gamma) % rep_times(allprobs.row(i), agsizes);
    sumfoo = sum(foo);
    l = l + log(sumfoo);
    phi = foo/sumfoo;
  }
  return l;
}

// [[Rcpp::export]]
double forward_cpp_sp(const arma::mat& allprobs, const arma::rowvec& delta, const arma::cube& Gamma, 
                      const IntegerVector& agsizes, const std::vector<int> tod)
{
  int N = allprobs.n_cols;
  int nObs = allprobs.n_rows;
  arma::rowvec foo(N);
  
  foo = delta % rep_times(allprobs.row(0), agsizes);
  double sumfoo = sum(foo);
  double l = log(sumfoo);
  arma::rowvec phi = foo/sumfoo;
  
  for (unsigned int i=1; i<nObs; i++){
    foo = (phi*Gamma.slice(tod[i])) % rep_times(allprobs.row(i), agsizes);
    sumfoo = sum(foo);
    l = l + log(sumfoo);
    phi = foo/sumfoo;
  }
  return l;
}

// local state probabilities for all tracks, with the same inputs as forward_cpp()
// forecast: Pr(S_t = j | x_1, ..., x_t-1) instead of Pr(S_t = j | x_1, ..., x_T)
// rows are not normalised
// [[Rcpp::export]]
arma::mat stateprobs_cpp(const arma::mat& allprobs, const arma::mat& Delta, const arma::cube& Gamma,
                         const IntegerVector& gidx, const IntegerVector& trackStart, bool forecast)
{
  int nObs = allprobs.n_rows;
  int N = allprobs.n_cols;
  int K = trackStart.size();

  arma::mat phi(nObs, N); // filtered state probabilities
  arma::mat probs(nObs, N);
  arma::rowvec foo(N);
  arma::colvec beta(N);

  for (int k = 0; k < K; k++){
    int start = trackStart[k];
    int end = (k < K - 1) ? trackStart[k + 1] - 1 : nObs - 1;

    // forward pass
    foo = Delta.row(k);
    for (int t = start; t <= end; t++){
      if (t > start) foo = phi.row(t - 1) * Gamma.slice(gidx[t]);
      probs.row(t) = foo; // forecast
      foo = foo % allprobs.row(t);
      phi.row(t) = foo / sum(foo);
    }
    if (forecast) continue;

    // backward pass with scaled beta: smoothing probabilities are proportional to phi * beta
    beta.ones();
    probs.row(end) = phi.row(end);
    for (int t = end - 1; t >= start; t--){
      beta = Gamma.slice(gidx[t + 1]) * (allprobs.row(t + 1).t() % beta);
      beta = beta / sum(beta);
      probs.row(t) = phi.row(t) % beta.t();
    }
  }
  return probs;
}

// Viterbi algorithm for all tracks, with the same inputs as forward_cpp()
// [[Rcpp::export]]
arma::colvec viterbi_cpp(const arma::mat& allprobs, const arma::mat& Delta, const arma::cube& Gamma,
                         const IntegerVector& gidx, const IntegerVector& trackStart)
{
  int nObs = allprobs.n_rows;
  int N = allprobs.n_cols;
  int K = trackStart.size();

  arma::mat xi(nObs, N);
  arma::rowvec foo(N);
  arma::colvec iv(nObs);

  for (int k = 0; k < K; k++){
    int start = trackStart[k];
    int end = (k < K - 1) ? trackStart[k + 1] - 1 : nObs - 1;

    foo = Delta.row(k) % allprobs.row(start);
    xi.row(start) = foo / sum(foo);

    for (int t = start + 1; t <= end; t++){
      for (int j = 0; j < N; j++){
        foo(j) = arma::max(xi.row(t - 1).t() % Gamma.slice(gidx[t]).col(j)) * allprobs(t, j);
      }
      xi.row(t) = foo / sum(foo);
    }

    // backtracking
    iv(end) = arma::index_max(xi.row(end));
    for (int t = end - 1; t >= start; t--){
      iv(t) = arma::index_max(Gamma.slice(gidx[t + 1]).col(iv(t + 1)) % xi.row(t).t());
    }
  }
  return iv + 1;
}
