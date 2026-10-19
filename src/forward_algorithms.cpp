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
double forward_cpp(const arma::mat& allprobs, 
                   const arma::mat& Delta, 
                   const arma::cube& Gamma,
                   const IntegerVector& gidx, 
                   const IntegerVector& trackStart)
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
