/*
 * Bayesian Factor Score computation via RcppArmadillo
 *
 * Author: Ricardo Rey-Sáez
 * email: ricardoreysaez95@gmail.com
 * Modification date: 19/02/2026
 *
 * Computes conditional posterior factor scores for each MCMC draw,
 * supporting both complete and incomplete (FIML) data patterns.
 */

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace arma;

// [[Rcpp::export]]
arma::cube compute_factor_scores_cpp(
    const arma::cube& post_draws,       // (n_iter x n_chains x n_vars)
    const arma::mat& data,              // (N x J) raw data
    const arma::uvec& lambda_idx,       // 0-based indices into var dimension
    const arma::uvec& sig_idx,          // 0-based indices for Sigma/Rho
    const arma::uvec& nu_idx,           // 0-based indices for nu (empty if no intercepts)
    int J, int M, int N,
    bool has_missing,
    const Rcpp::List& obs_patterns      // list of integer vectors (0-based observed indices per person)
) {

  const int n_iter   = post_draws.n_rows;
  const int n_chains = post_draws.n_cols;
  const int P        = N * M;

  // Output: (n_iter x n_chains x N*M) in column-major person order
  arma::cube eta_out(n_iter, n_chains, P, arma::fill::zeros);

  const arma::mat I_M = arma::eye<arma::mat>(M, M);
  const bool has_nu   = (nu_idx.n_elem > 0);

  for (int ch = 0; ch < n_chains; ++ch) {
    for (int s = 0; s < n_iter; ++s) {

      // ── Extract parameters for this draw ──

      // Lambda (J x M), column-major from Stan
      arma::vec lambda_vec(lambda_idx.n_elem);
      for (arma::uword k = 0; k < lambda_idx.n_elem; ++k) {
        lambda_vec(k) = post_draws(s, ch, lambda_idx(k));
      }
      arma::mat L_s = arma::reshape(lambda_vec, J, M);

      // Sigma or Rho (J x J)
      arma::vec sig_vec(sig_idx.n_elem);
      for (arma::uword k = 0; k < sig_idx.n_elem; ++k) {
        sig_vec(k) = post_draws(s, ch, sig_idx(k));
      }
      arma::mat Sigma_s = arma::reshape(sig_vec, J, J);

      // Intercepts
      arma::vec nu_s(J, arma::fill::zeros);
      if (has_nu) {
        for (arma::uword k = 0; k < nu_idx.n_elem; ++k) {
          nu_s(k) = post_draws(s, ch, nu_idx(k));
        }
      }

      // ── Uniquenesses: psi = diag(Sigma) - rowSums(Lambda^2) ──
      arma::vec h2_s = arma::sum(L_s % L_s, 1);  // row sums of squared loadings
      arma::vec psi_s = arma::diagvec(Sigma_s) - h2_s;
      arma::vec Psi_inv = 1.0 / psi_s;

      // ── Complete-data path: vectorized over all N persons ──
      if (!has_missing) {
        // L_s' * Psi^{-1} * L_s  (M x M)
        arma::mat L_weighted = L_s.each_col() % Psi_inv;  // (J x M) with rows scaled
        arma::mat LtPsiInvL  = L_s.t() * L_weighted;      // (M x M)

        // Sigma_eta = (I_M + L' Psi^{-1} L)^{-1}
        arma::mat Sigma_eta = arma::inv_sympd(I_M + LtPsiInvL);

        // Cholesky of Sigma_eta for sampling
        arma::mat chol_Sig;
        bool chol_ok = arma::chol(chol_Sig, Sigma_eta);
        if (!chol_ok) {
          chol_Sig = arma::chol(Sigma_eta);
        }

        // Center data: Y_centered = data - nu' (each row minus nu)
        arma::mat Y_centered = data.each_row() - nu_s.t();  // (N x J)

        // Conditional means for all persons: (N x M)
        // mu = Y_centered * diag(Psi_inv) * L_s * Sigma_eta
        arma::mat mu_all = Y_centered * L_weighted * Sigma_eta;

        // Sample: eta = mu + Z * chol(Sigma_eta)
        arma::mat Z = arma::randn<arma::mat>(N, M);
        arma::mat eta_NxM = mu_all + Z * chol_Sig;

        // Store in column-major order: all persons for factor 1, then factor 2, etc.
        for (int m = 0; m < M; ++m) {
          int idx_start = m * N;
          for (int i = 0; i < N; ++i) {
            eta_out(s, ch, idx_start + i) = eta_NxM(i, m);
          }
        }

        // ── Missing-data path: person by person ──
      } else {
        for (int i = 0; i < N; ++i) {
          Rcpp::IntegerVector obs_i_R = obs_patterns[i];
          int n_obs_i = obs_i_R.size();

          if (n_obs_i == 0) {
            // No observed data: sample from prior N(0, I)
            arma::vec eta_prior = arma::randn<arma::vec>(M);
            for (int m = 0; m < M; ++m) {
              eta_out(s, ch, m * N + i) = eta_prior(m);
            }
            continue;
          }

          // Convert to arma uvec (0-based)
          arma::uvec obs_i(n_obs_i);
          for (int k = 0; k < n_obs_i; ++k) {
            obs_i(k) = static_cast<arma::uword>(obs_i_R[k]);
          }

          // Subset Lambda, nu, Psi_inv for observed variables
          arma::mat L_i       = L_s.rows(obs_i);        // (n_obs x M)
          arma::vec nu_i      = nu_s.elem(obs_i);
          arma::vec psi_inv_i = Psi_inv.elem(obs_i);

          // Observed data for this person (convert row to col vec first)
          arma::vec y_full = data.row(i).t();
          arma::vec y_i = y_full.elem(obs_i);

          // Conditional posterior
          arma::mat L_w_i    = L_i.each_col() % psi_inv_i;
          arma::mat LtPsiInvL_i = L_i.t() * L_w_i;
          arma::mat Sigma_eta_i = arma::inv_sympd(I_M + LtPsiInvL_i);

          // Centered data
          arma::vec y_centered_i = y_i - nu_i;

          // Conditional mean
          arma::vec mu_i = Sigma_eta_i * (L_w_i.t() * y_centered_i);

          // Sample
          arma::mat chol_i;
          bool chol_ok = arma::chol(chol_i, Sigma_eta_i);
          if (!chol_ok) {
            chol_i = arma::chol(Sigma_eta_i + 1e-8 * I_M);
          }
          arma::vec eta_i = mu_i + chol_i.t() * arma::randn<arma::vec>(M);

          // Store
          for (int m = 0; m < M; ++m) {
            eta_out(s, ch, m * N + i) = eta_i(m);
          }
        }
      }
    }
  }

  return eta_out;
}
