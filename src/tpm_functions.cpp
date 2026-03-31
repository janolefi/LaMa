// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::cube tpm_g_cpp(const arma::mat& Z, const arma::mat& beta, 
                     const arma::uword N, const bool byrow)
{
  // get number of data points
  arma::uword n = Z.n_rows;
  
  arma::mat expEta(n, N*(N-1));
  expEta = arma::exp(Z * beta.t());
  
  // initialize matrix and cube
  arma::mat G(N,N);
  arma::cube Gamma(N,N,n);
  
  if(byrow == true) {
    // if byrow is true Gamma should be created rowwise: extra transpose step
    for(arma::uword i = 0; i < n; i++) {
      // Create a diagonal matrix G with N as diagonal elements
      G = arma::eye<arma::mat>(N, N);
      // Modify non-diagonal elements of G using values from expEta
      G(arma::find(G == 0)) = expEta.row(i).t();
      // transpose for byrow
      G = G.t();
      // Normalize each row of the matrix G by dividing it by the sum of its elements
      G = arma::normalise(G, 1, 1);
      Gamma.slice(i) = G;
    }
  } else {
    for(arma::uword i = 0; i < n; i++) {
      // Create a diagonal matrix G with N as diagonal elements
      G = arma::eye<arma::mat>(N, N);
      // Modify non-diagonal elements of G using values from expEta
      G(arma::find(G == 0)) = expEta.row(i).t();
      // Normalize each row of the matrix G by dividing it by the sum of its elements
      G = arma::normalise(G, 1, 1);
      Gamma.slice(i) = G;
    }
  }
  return Gamma;
}

// [[Rcpp::export]]
arma::cube tpm_g2_cpp(const arma::mat& Eta,
                      const arma::uword N, 
                      const bool byrow,
                      const arma::uvec& ref)
{
  // get number of data points
  arma::uword n = Eta.n_rows;
  
  arma::mat expEta = arma::exp(Eta);
  
  // initialize matrix and cube
  arma::cube Gamma(N,N,n);
  arma::mat G(N, N, arma::fill::zeros);  // Create G matrix
  
  if(byrow == true) {
    // if byrow is true Gamma should be created rowwise: extra transpose step
    for (arma::uword i = 0; i < n; ++i) {
      // Create a matrix with 1s in the ref column and 0s elsewhere
      G.zeros();  // Reset in-place without reallocating
      for (arma::uword j = 0; j < N; ++j) {
        G(j, ref[j] - 1) = 1.0;
      }
      // Modify non-diagonal elements of G using values from expEta
      G(arma::find(G == 0)) = expEta.row(i).t();
      // transpose for byrow
      G = G.t();
      // Normalize each row of the matrix G by dividing it by the sum of its elements
      G = arma::normalise(G, 1, 1);
      Gamma.slice(i) = G;
    }
  } else {
    for(arma::uword i = 0; i < n; i++) {
      // Create a matrix with 1s in the ref column and 0s elsewhere
      G.zeros();  // Reset in-place without reallocating
      for (arma::uword j = 0; j < N; ++j) {
        G(j, ref[j] - 1) = 1.0;
      }
      // Modify non-diagonal elements of G using values from expEta
      G(arma::find(G == 0)) = expEta.row(i).t();
      // Normalize each row of the matrix G by dividing it by the sum of its elements
      G = arma::normalise(G, 1, 1);
      Gamma.slice(i) = G;
    }
  }
  return Gamma;
}

// [[Rcpp::export]]
arma::cube tpm_g3_cpp(const arma::mat& Eta,
                      const arma::uword N,
                      const arma::uvec& ref,   // 1-based from R
                      const bool byrow = false)
{
  const arma::uword n = Eta.n_rows;
  // if (ref.n_elem != N) Rcpp::stop("ref must have length N");
  // const arma::uword n_trans = N * (N - 1);
  // if (Eta.n_cols != n_trans) Rcpp::stop("Eta must have N*(N-1) columns");
  
  arma::uvec ref0 = ref - 1; // convert to 0-based indexing
  if (arma::any(ref0 >= N)) Rcpp::stop("ref values must be in 1:N");
  
  arma::cube Gamma(N, N, n);
  const arma::mat expEta = arma::exp(Eta);
  
  for (arma::uword t = 0; t < n; ++t) {
    arma::mat G(N, N, arma::fill::ones); // temp slice
    arma::uword col_ind = 0;
    
    for (arma::uword i = 0; i < N; ++i) {
      for (arma::uword j = 0; j < N; ++j) {
        if (!byrow) {
          if (i != ref0[j]) G(j, i) = expEta(t, col_ind++);
        } else {
          if (j != ref0[i]) G(i, j) = expEta(t, col_ind++);
        }
      }
    }
    
    // normalize rows
    G = arma::normalise(G, 1, 1);
    Gamma.slice(t) = G;
  }
  
  return Gamma;
}



// [[Rcpp::export]]
arma::cube semigroup_cpp(const arma::mat& Q, const std::vector<double>& times) {
  int N = Q.n_cols;
  int n = times.size();
  
  arma::cube Gamma(N, N, n);
  
  for (unsigned int i = 0; i < n; i++) {
    Gamma.slice(i) = arma::expmat(Q * times.at(i));
  }
  
  return Gamma;
}

// [[Rcpp::export]]
arma::cube semigroup_g_cpp(const arma::cube& Q, const std::vector<double>& times) {
  int N = Q.n_cols;
  int n = times.size();
  
  arma::cube Gamma(N, N, n);
  
  for (unsigned int i = 0; i < n; i++) {
    Gamma.slice(i) = arma::expmat(Q.slice(i) * times.at(i));
  }
  
  return Gamma;
}

// [[Rcpp::export]]
arma::mat tpm_thinned_t_cpp(const arma::cube& Gamma, const int t)
{
  arma::uword L = Gamma.n_slices;
  arma::uword N = Gamma.n_rows;
  
  arma::mat GammaT(N,N); 
  // first half of the product form t to L
  GammaT = Gamma.slice(t-1);
  for(arma::uword i=t; i<L; i++) {
    GammaT = GammaT * Gamma.slice(i);
  }
  // second half of the product from 1 to t-1
  for(arma::uword i=0; i<t-1; i++) {
    GammaT = GammaT * Gamma.slice(i);
  }
  return GammaT;
}
