/*
 * Efficient Rotation-Sign-Permutation Algorithm
 *
 * Author: Ricardo Rey-Sáez
 * email: ricardoreysaez95@gmail.com
 * Modification date: 19/02/2026

 * Solve rotation indeterminacy in factor loadings
 */

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <vector>
#include <limits>
#include <algorithm>
#include <cmath>

using namespace arma;

// Hungarian Algorithm for Linear Assignment Problem
// Adapted from here:
// https://github.com/kth-competitive-programming/kactl/blob/main/content/graph/WeightedMatching.h
std::vector<int> solve_hungarian(const arma::mat& cost_matrix) {

  // Small check
  if (cost_matrix.is_empty()) return {};

  int n_rows = cost_matrix.n_rows;
  int n_cols = cost_matrix.n_cols;
  // Internal logic uses 1-based indexing and a dummy worker at index 0
  int n = n_rows + 1;
  int m = n_cols + 1;
  // Auxiliar objects
  std::vector<double> u(n, 0), v(m, 0);
  std::vector<int> p(m, 0), ans(n - 1, 0);
  const double INF = std::numeric_limits<double>::infinity();

  for (int i = 1; i < n; ++i) {
    p[0] = i;
    int j0 = 0; // dummy worker
    std::vector<double> dist(m, INF);
    std::vector<int> pre(m, -1);
    std::vector<bool> done(m + 1, false);

    do { // Dijkstra phase
      done[j0] = true;
      int i0 = p[j0];
      int j1 = 0;
      double delta = INF;

      for (int j = 1; j < m; ++j) {
        if (!done[j]) {
          // Calculate cost (adjusting for 0-based arma matrix)
          double cur = cost_matrix(i0 - 1, j - 1) - u[i0] - v[j];

          if (cur < dist[j]) {
            dist[j] = cur;
            pre[j] = j0;
          }
          if (dist[j] < delta) {
            delta = dist[j];
            j1 = j;
          }
        }
      }

      for (int j = 0; j < m; ++j) {
        if (done[j]) {
          u[p[j]] += delta;
          v[j] -= delta;
        } else {
          dist[j] -= delta;
        }
      }
      j0 = j1;
    } while (p[j0] != 0);

    while (j0 != 0) { // Update alternating path
      int j1 = pre[j0];
      p[j0] = p[j1];
      j0 = j1;
    }
  }

  // Convert to 0-based indices for C++ output
  for (int j = 1; j < m; ++j) {
    if (p[j] != 0) {
      ans[p[j] - 1] = j - 1;
    }
  }
  return ans;
}

// Exact and Efficient Rotation-Sign-Permutation (RSP) post-processing algorithm
// [[Rcpp::export]]
Rcpp::List rsp_exact_eff(const arma::mat& Lambda_tilde,
                         const arma::mat& Lambda_star_init,
                         int J, int M,
                         int maxIter, double threshold) {

  const int n_draws = (int)Lambda_tilde.n_rows;
  const int n_pars  = J * M;
  const double adj_threshold = threshold * n_draws * n_pars;

  // Initialize Lambda_star from the reference matrix
  arma::mat Lambda_star = Lambda_star_init;

  // Output containers
  arma::mat s_vectors(n_draws, M, arma::fill::ones);
  arma::imat nu_vectors(n_draws, M, arma::fill::zeros);
  for (int m = 0; m < M; ++m) nu_vectors.col(m).fill(m + 1); // 1-based identity

  std::vector<double> objective_values(maxIter, 0.0);
  arma::mat Lambda_hat_vec(n_draws, Lambda_tilde.n_cols, arma::fill::zeros);

  bool criterion = true;
  int n_iter = 0;

  // Pre-compute buffer objects
  arma::mat Lambda_star_new(J, M);
  arma::mat Lambda_tilde_s(J, M);
  arma::mat Lambda_tilde_signed_s(J, M);
  arma::mat Lambda_hat_s(J, M);
  arma::vec sqnorm_Lambda_tilde(M);
  arma::vec sqnorm_Lambda_star(M);
  arma::mat innerprod(M, M);
  arma::mat C_s(M, M);
  arma::imat s_opt_pair(M, M);
  arma::rowvec s_t(M);

 // Main While loop
  while (criterion && n_iter < maxIter) {
    ++n_iter;

    Lambda_star_new.zeros();
    double objective = 0.0;

    // Column norms squared for posterior means
    sqnorm_Lambda_star = arma::sum(arma::square(Lambda_star), 0).t();

    for (int s = 0; s < n_draws; ++s) {

      // Lambda_tilde_s = reshape draw s
      Lambda_tilde_s = arma::reshape(Lambda_tilde.row(s).t(), J, M);

      // Column norms squared for Lambda_tilde_s
      sqnorm_Lambda_tilde = arma::sum(arma::square(Lambda_tilde_s), 0).t();

      // Cross-product between matrices
      // (i.e., trace-elements without sign adjustment
      innerprod = Lambda_tilde_s.t() * Lambda_star;

      // Fill cost matrix and best signs
      for (int m_orig = 0; m_orig < M; ++m_orig) {
        for (int m_hat = 0; m_hat < M; ++m_hat) {

          // Compute Squared Frobenius difference for all combinations
          const double base = sqnorm_Lambda_tilde(m_orig) + sqnorm_Lambda_star(m_hat);
          const double dist_posit = base - 2.0 * innerprod(m_orig, m_hat);
          const double dist_negat = base + 2.0 * innerprod(m_orig, m_hat);

          // If positive sign minimizes the Frobenius squared norm, select
          if (dist_posit <= dist_negat) {
            C_s(m_orig, m_hat) = dist_posit;
            s_opt_pair(m_orig, m_hat) =  1;
          // If not, choose negative sign
          } else {
            C_s(m_orig, m_hat) = dist_negat;
            s_opt_pair(m_orig, m_hat) = -1;
          }
        }
      }

      // Hungarian
      std::vector<int> perm = solve_hungarian(C_s);

      // Best sign vector and cost value
      double current_cost = 0.0;
      for (int m = 0; m < M; ++m) {
        const int assigned_hat = perm[m];
        s_t(m) = s_opt_pair(m, assigned_hat);
        current_cost += C_s(m, assigned_hat);
      }

      // Apply signs
      Lambda_tilde_signed_s = Lambda_tilde_s;
      Lambda_tilde_signed_s.each_row() %= s_t;

      // Apply permutation: draw col m goes to hat col perm[m]
      for (int m = 0; m < M; ++m) {
        Lambda_hat_s.col(perm[m]) = Lambda_tilde_signed_s.col(m);
      }

      // Accumulate
      objective += current_cost;
      Lambda_star_new += Lambda_hat_s;

      // Store states
      s_vectors.row(s) = s_t;
      for (int m = 0; m < M; ++m) nu_vectors(s, perm[m]) = m + 1; // 1-based
      Lambda_hat_vec.row(s) = arma::vectorise(Lambda_hat_s).t();
    }

    objective_values[n_iter - 1] = objective;

    if (n_iter > 1) {
      const double improv = objective_values[n_iter - 2] - objective_values[n_iter - 1];
      if (improv >= 0.0 && improv < adj_threshold) criterion = false;
    }

    Lambda_star = Lambda_star_new / n_draws;
  }

  objective_values.resize(n_iter);

  return Rcpp::List::create(
    Rcpp::Named("Lambda_hat_mcmc") = Lambda_hat_vec,
    Rcpp::Named("sign_vectors")    = s_vectors,
    Rcpp::Named("perm_vectors")    = nu_vectors,
    Rcpp::Named("Lambda_star")     = Lambda_star,
    Rcpp::Named("objective")       = objective_values
  );
}


