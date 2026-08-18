
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//' Pairwise squared distances between two point sets (rows)
//' @param A n×d
//' @param B m×d
//' @return n×m matrix of ||a_i - b_j||^2
// [[Rcpp::export]]
arma::mat pairwise_sqdist_cpp(const arma::mat& A, const arma::mat& B) {
  int n = A.n_rows, m = B.n_rows;
  arma::vec An = sum(square(A), 1);
  arma::vec Bn = sum(square(B), 1);
  arma::mat C = repmat(An, 1, m) + repmat(Bn.t(), n, 1) - 2.0 * (A * B.t());
  // Floating-point cancellation can yield tiny negatives for near-coincident
  // points; clamp to 0 to match the R fallback (dkge-transport.R).
  C.elem(arma::find(C < 0)).zeros();
  return C;
}
