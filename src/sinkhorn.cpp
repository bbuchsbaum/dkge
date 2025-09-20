#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

inline double log_sum_exp(const arma::rowvec& x) {
  double max_val = x.max();
  if (!std::isfinite(max_val)) {
    return max_val;
  }
  arma::rowvec shifted = x - max_val;
  return max_val + std::log(arma::sum(arma::exp(shifted)));
}

inline double log_sum_exp_col(const arma::colvec& x) {
  double max_val = x.max();
  if (!std::isfinite(max_val)) {
    return max_val;
  }
  arma::colvec shifted = x - max_val;
  return max_val + std::log(arma::sum(arma::exp(shifted)));
}

} // namespace

// [[Rcpp::export]]
Rcpp::List sinkhorn_plan_cpp(Rcpp::NumericMatrix C,
                             Rcpp::NumericVector mu,
                             Rcpp::NumericVector nu,
                             double epsilon,
                             int max_iter,
                             double tol,
                             Rcpp::Nullable<Rcpp::NumericVector> log_u_init = R_NilValue,
                             Rcpp::Nullable<Rcpp::NumericVector> log_v_init = R_NilValue,
                             bool keep_duals = true) {
  const int n = C.nrow();
  const int m = C.ncol();

  if (mu.size() != n) {
    Rcpp::stop("Length of mu must match number of rows in C");
  }
  if (nu.size() != m) {
    Rcpp::stop("Length of nu must match number of columns in C");
  }
  if (epsilon <= 0) {
    Rcpp::stop("epsilon must be positive");
  }
  if (max_iter <= 0) {
    Rcpp::stop("max_iter must be positive");
  }

  arma::mat logK = -arma::mat(C.begin(), n, m, false, true) / epsilon;
  arma::vec mu_vec(mu.begin(), n, false, true);
  arma::vec nu_vec(nu.begin(), m, false, true);

  arma::vec log_mu = arma::log(mu_vec);
  arma::vec log_nu = arma::log(nu_vec);

  arma::vec log_u(n, arma::fill::zeros);
  arma::vec log_v(m, arma::fill::zeros);

  if (log_u_init.isNotNull()) {
    arma::vec init = arma::vec(Rcpp::as<Rcpp::NumericVector>(log_u_init).begin(), n, false, true);
    if (init.n_elem == static_cast<unsigned>(n)) {
      log_u = init;
    }
  }
  if (log_v_init.isNotNull()) {
    arma::vec init = arma::vec(Rcpp::as<Rcpp::NumericVector>(log_v_init).begin(), m, false, true);
    if (init.n_elem == static_cast<unsigned>(m)) {
      log_v = init;
    }
  }

  arma::mat plan(n, m, arma::fill::zeros);
  arma::vec row_sums(n, arma::fill::zeros);
  arma::vec col_sums(m, arma::fill::zeros);

  const double* logK_ptr = logK.memptr();

  int iter_taken = max_iter;
  for (int it = 0; it < max_iter; ++it) {
    // Update log_u
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (int i = 0; i < n; ++i) {
      double max_val = -std::numeric_limits<double>::infinity();
      for (int j = 0; j < m; ++j) {
        double val = logK_ptr[i + j * n] + log_v(j);
        if (val > max_val) {
          max_val = val;
        }
      }
      double accum = 0.0;
      for (int j = 0; j < m; ++j) {
        double val = logK_ptr[i + j * n] + log_v(j);
        accum += std::exp(val - max_val);
      }
      log_u(i) = log_mu(i) - (max_val + std::log(accum));
    }

    // Update log_v
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (int j = 0; j < m; ++j) {
      const double* col_ptr = logK_ptr + j * n;
      double max_val = -std::numeric_limits<double>::infinity();
      for (int i = 0; i < n; ++i) {
        double val = col_ptr[i] + log_u(i);
        if (val > max_val) {
          max_val = val;
        }
      }
      double accum = 0.0;
      for (int i = 0; i < n; ++i) {
        double val = col_ptr[i] + log_u(i);
        accum += std::exp(val - max_val);
      }
      log_v(j) = log_nu(j) - (max_val + std::log(accum));
    }

    if ((it % 5 == 4) || (it == max_iter - 1)) {
      row_sums.zeros();
      col_sums.zeros();
      #ifdef _OPENMP
      #pragma omp parallel for schedule(static)
      #endif
      for (int j = 0; j < m; ++j) {
        const double* col_ptr = logK_ptr + j * n;
        double log_v_j = log_v(j);
        double col_sum_local = 0.0;
        for (int i = 0; i < n; ++i) {
          double val = std::exp(col_ptr[i] + log_u(i) + log_v_j);
          plan(i, j) = val;
          #ifdef _OPENMP
          #pragma omp atomic
          #endif
          row_sums(i) += val;
          col_sum_local += val;
        }
        col_sums(j) = col_sum_local;
      }

      double err1 = arma::abs(row_sums - mu_vec).max();
      double err2 = arma::abs(col_sums - nu_vec).max();
      if (std::max(err1, err2) < tol) {
        iter_taken = it + 1;
        Rcpp::List out = Rcpp::List::create(
          Rcpp::Named("plan") = plan,
          Rcpp::Named("iterations") = iter_taken
        );
        if (keep_duals) {
          out["log_u"] = log_u;
          out["log_v"] = log_v;
        }
        return out;
      }
    }
  }

  // If convergence criterion not met, return last plan
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int j = 0; j < m; ++j) {
    const double* col_ptr = logK_ptr + j * n;
    double log_v_j = log_v(j);
    for (int i = 0; i < n; ++i) {
      double val = std::exp(col_ptr[i] + log_u(i) + log_v_j);
      plan(i, j) = val;
    }
  }
  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("plan") = plan,
    Rcpp::Named("iterations") = max_iter
  );
  if (keep_duals) {
    out["log_u"] = log_u;
    out["log_v"] = log_v;
  }
  return out;
}
