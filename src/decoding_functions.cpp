// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// local state probabilities for all tracks, with the same inputs as forward_cpp()
// forecast: Pr(S_t = j | x_1, ..., x_t-1) instead of Pr(S_t = j | x_1, ..., x_T)
// rows are not normalised
// [[Rcpp::export]]
arma::mat stateprobs_cpp(const arma::mat& allprobs, 
                         const arma::mat& Delta, 
                         const arma::cube& Gamma,
                         const IntegerVector& gidx, 
                         const IntegerVector& trackStart, 
                         bool forecast)
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
arma::colvec viterbi_cpp(const arma::mat& allprobs, 
                         const arma::mat& Delta, 
                         const arma::cube& Gamma,
                         const IntegerVector& gidx, 
                         const IntegerVector& trackStart)
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
