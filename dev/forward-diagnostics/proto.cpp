// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// single non-AD kernel: Delta (nTracks x N), Gamma cube (N x N x nSlices),
// gidx (length nObs, 0-based slice of the transition INTO t), trackStart (0-based, length nTracks)
// [[Rcpp::export]]
double forward_cpp_proto(const arma::mat& allprobs, const arma::mat& Delta, const arma::cube& Gamma,
                         const IntegerVector& gidx, const IntegerVector& trackStart) {
  const int nObs = allprobs.n_rows;
  const int K = trackStart.size();
  double l = 0;
  arma::rowvec phi, foo;
  int k = 0;
  for (int t = 0; t < nObs; t++) {
    if (k < K && t == trackStart[k]) {
      foo = Delta.row(k) % allprobs.row(t);
      k++;
    } else {
      foo = (phi * Gamma.slice(gidx[t])) % allprobs.row(t);
    }
    double s = arma::accu(foo);
    l += std::log(s);
    phi = foo / s;
  }
  return l;
}
