// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// Forward algorithm and decoding for HSMMs, on the output of prepare_hsmm() (see R/hsmm.R)
// - stay.row(t): probabilities of all M states to stay for one more time point, for the step into t
// - omega.slice(oidx[t]): conditional transition probabilities for the step into t
// - Delta.row(k): initial distribution of the M states for track k
// - start, end: first and last observation used in each track (all indices 0-based)
// The transition matrices (M x M) are never built: each step uses their structure

// Sizes of the state aggregates, and where they are
struct Aggregates {
  int N, M;
  std::vector<int> I_low, I_high; // lowest and highest state of each aggregate
  std::vector<int> state;         // aggregate of each of the M states

  Aggregates(const IntegerVector& Ni) : N(Ni.size()), M(sum(Ni)) {
    for (int i = 0; i < N; i++){
      I_low.push_back(state.size());
      state.insert(state.end(), Ni[i], i);
      I_high.push_back(state.size() - 1);
    }
  }
};

// Forward step, phi * Gamma: where does the mass of each state go? As hsmm_step() in R
arma::rowvec step_forward(const arma::rowvec& phi, const arma::rowvec& stay,
                          const arma::mat& omega, const Aggregates& A)
{
  arma::rowvec left(A.N, fill::zeros); // mass that leaves each aggregate
  arma::rowvec phi_new(A.M, fill::zeros);

  for (int s = 0; s < A.M; s++){
    int i = A.state[s];
    double kept = phi(s) * stay(s);
    left(i) += phi(s) - kept;
    int up = std::min(s + 1, A.I_high[i]); // one state up, or wait in the highest state
    phi_new(up) += kept;
  }
  arma::rowvec arrive = left * omega;
  for (int j = 0; j < A.N; j++) phi_new(A.I_low[j]) += arrive(j);
  return phi_new;
}

// Backward step, Gamma * b: what is the value of each state, if b is the value after the step?
arma::colvec step_backward(const arma::colvec& b, const arma::rowvec& stay,
                           const arma::mat& omega, const Aggregates& A)
{
  arma::colvec b_low(A.N); // value of arriving in each aggregate
  for (int j = 0; j < A.N; j++) b_low(j) = b(A.I_low[j]);
  arma::colvec leave = omega * b_low; // value of leaving each aggregate

  arma::colvec b_new(A.M);
  for (int s = 0; s < A.M; s++){
    int i = A.state[s];
    int up = std::min(s + 1, A.I_high[i]);
    b_new(s) = stay(s) * b(up) + (1 - stay(s)) * leave(i);
  }
  return b_new;
}

// Viterbi step: as the forward step, with the largest mass instead of the sum of masses
// from: for each state, the state before the step that its largest mass comes from
arma::rowvec step_viterbi(const arma::rowvec& xi, const arma::rowvec& stay,
                          const arma::mat& omega, const Aggregates& A, arma::Row<int>& from)
{
  arma::rowvec left(A.N);          // largest mass that leaves each aggregate
  std::vector<int> left_from(A.N); // and its state
  left.fill(-1);
  arma::rowvec xi_new(A.M);
  xi_new.fill(-1);

  for (int s = 0; s < A.M; s++){
    int i = A.state[s];
    double kept = xi(s) * stay(s);
    double gone = xi(s) - kept;
    if (gone > left(i)) { left(i) = gone; left_from[i] = s; }
    int up = std::min(s + 1, A.I_high[i]);
    if (kept > xi_new(up)) { xi_new(up) = kept; from(up) = s; }
  }
  for (int j = 0; j < A.N; j++){
    int low = A.I_low[j];
    for (int i = 0; i < A.N; i++){
      double arrive = left(i) * omega(i, j);
      if (arrive > xi_new(low)) { xi_new(low) = arrive; from(low) = left_from[i]; }
    }
  }
  return xi_new;
}

// state-dependent probabilities of all M states at time t
arma::rowvec expand(const arma::mat& allprobs, int t, const Aggregates& A)
{
  arma::rowvec probs(A.M);
  for (int s = 0; s < A.M; s++) probs(s) = allprobs(t, A.state[s]);
  return probs;
}

// log-likelihood of all tracks
// [[Rcpp::export]]
double forward_hsmm_cpp(const arma::mat& allprobs, const arma::mat& Delta,
                        const arma::mat& stay, const arma::cube& omega, const IntegerVector& oidx,
                        const IntegerVector& Ni, const IntegerVector& start, const IntegerVector& end)
{
  Aggregates A(Ni);
  double l = 0;

  for (int k = 0; k < start.size(); k++){
    arma::rowvec phi = Delta.row(k);
    for (int t = start[k]; t <= end[k]; t++){
      if (t > start[k]) phi = step_forward(phi, stay.row(t), omega.slice(oidx[t]), A);
      phi = phi % expand(allprobs, t, A);
      double sumphi = sum(phi);
      l += log(sumphi);
      phi = phi / sumphi;
    }
  }
  return l;
}

// local state probabilities of the N states for all tracks, NA for unused observations
// forecast: Pr(S_t = j | x_1, ..., x_t-1) instead of Pr(S_t = j | x_1, ..., x_T)
// [[Rcpp::export]]
arma::mat stateprobs_hsmm_cpp(const arma::mat& allprobs, const arma::mat& Delta,
                              const arma::mat& stay, const arma::cube& omega, const IntegerVector& oidx,
                              const IntegerVector& Ni, const IntegerVector& start, const IntegerVector& end,
                              bool forecast)
{
  Aggregates A(Ni);
  arma::mat probs(allprobs.n_rows, A.N);
  probs.fill(NA_REAL);
  arma::mat filtered(A.M, allprobs.n_rows); // one column per observation

  // sums the M states in each aggregate, normalised
  auto save = [&](int t, const arma::vec& p) {
    probs.row(t).zeros();
    for (int s = 0; s < A.M; s++) probs(t, A.state[s]) += p(s);
    probs.row(t) /= sum(probs.row(t));
  };

  for (int k = 0; k < start.size(); k++){
    // forward pass
    arma::rowvec phi = Delta.row(k);
    for (int t = start[k]; t <= end[k]; t++){
      if (t > start[k]) phi = step_forward(phi, stay.row(t), omega.slice(oidx[t]), A);
      if (forecast) save(t, phi.t());
      phi = phi % expand(allprobs, t, A);
      phi = phi / sum(phi);
      filtered.col(t) = phi.t();
    }
    if (forecast) continue;

    // backward pass with scaled beta: smoothing probabilities are proportional to phi * beta
    arma::colvec beta(A.M, fill::ones);
    save(end[k], filtered.col(end[k]));
    for (int t = end[k] - 1; t >= start[k]; t--){
      beta = step_backward(expand(allprobs, t + 1, A).t() % beta, stay.row(t + 1), omega.slice(oidx[t + 1]), A);
      beta = beta / sum(beta);
      save(t, filtered.col(t) % beta);
    }
  }
  return probs;
}

// most probable path through the M states for all tracks: the state (1, ..., N) of each
// observation, NA for unused observations
// [[Rcpp::export]]
IntegerVector viterbi_hsmm_cpp(const arma::mat& allprobs, const arma::mat& Delta,
                               const arma::mat& stay, const arma::cube& omega, const IntegerVector& oidx,
                               const IntegerVector& Ni, const IntegerVector& start, const IntegerVector& end)
{
  Aggregates A(Ni);
  IntegerVector states(allprobs.n_rows, NA_INTEGER);
  arma::Mat<int> from(A.M, allprobs.n_rows); // one column per observation
  arma::Row<int> from_t(A.M, fill::zeros);

  for (int k = 0; k < start.size(); k++){
    arma::rowvec xi = Delta.row(k);
    for (int t = start[k]; t <= end[k]; t++){
      if (t > start[k]) {
        xi = step_viterbi(xi, stay.row(t), omega.slice(oidx[t]), A, from_t);
        from.col(t) = from_t.t();
      }
      xi = xi % expand(allprobs, t, A);
      xi = xi / sum(xi);
    }

    // backtracking
    int s = index_max(xi);
    for (int t = end[k]; t >= start[k]; t--){
      states[t] = A.state[s] + 1;
      if (t > start[k]) s = from(s, t);
    }
  }
  return states;
}
