// Unified Null (Independence) Model for Bayesian Fit Indices
//
// Flag:
//   model_type: 1 = raw (means + variances), 2 = cov (variances only)

functions {
  // Sufficient statistics log-likelihood with means (raw model)
  real sufficient_ll_raw(int N, int J, vector m_obs, vector nu,
                         matrix S_obs, matrix Sigma) {
    vector[J] diff = m_obs - nu;
    real logdet = log_determinant(Sigma);
    real tr_part = trace(mdivide_right_spd(S_obs, Sigma));
    real mean_quad = dot_product(diff, mdivide_left_spd(Sigma, diff));
    return -0.5 * N * (logdet + tr_part + mean_quad);
  }

  // Sufficient statistics log-likelihood without means (cov model)
  real sufficient_ll_cov(int N, matrix S, matrix Sigma) {
    return -0.5 * N * (log_determinant(Sigma) + trace(mdivide_right_spd(S, Sigma)));
  }
}

data {
  int<lower=0> N_complete;
  int<lower=0> N_incomplete;
  int<lower=0> J;

  // Model configuration flag
  int<lower=1, upper=2> model_type;  // 1=raw, 2=cov

  // Sufficient statistics
  vector[J] m_obs;
  cov_matrix[J] S_obs;

  // Incomplete data for FIML
  matrix[N_incomplete, J] Y_miss;
  int<lower=0> max_obs;
  array[N_incomplete] int<lower=0> n_obs;
  array[N_incomplete, max_obs] int<lower=0, upper=J> obs_idx;

  // Priors
  vector<lower=0>[2] pr_nu;       // Normal(mean, sd) for nu
  vector<lower=0>[3] pr_sigma;    // Student-t(df, loc, scale) for sigma
}

parameters {
  vector[model_type == 1 ? J : 0] nu;
  vector<lower=0>[J] sigma;
}

transformed parameters {
  cov_matrix[J] Sigma = diag_matrix(square(sigma));
}

model {
  // Priors
  if (model_type == 1) {
    target += normal_lpdf(nu | pr_nu[1], pr_nu[2]);
  }
  target += student_t_lpdf(sigma | pr_sigma[1], pr_sigma[2], pr_sigma[3]);
  target += J * -student_t_lccdf(0 | pr_sigma[1], pr_sigma[2], pr_sigma[3]);

  // Log-Likelihood: Complete observations
  if (N_complete > 0) {
    if (model_type == 1) {
      target += sufficient_ll_raw(N_complete, J, m_obs, nu, S_obs, Sigma);
    } else {
      target += sufficient_ll_cov(N_complete, S_obs, Sigma);
    }
  }

  // Log-Likelihood: Incomplete observations (FIML)
  // For null model with diagonal Sigma, use element-wise normal_lpdf
  if (N_incomplete > 0) {
    for (i in 1:N_incomplete) {
      int ni = n_obs[i];
      for (k in 1:ni) {
        int j = obs_idx[i, k];
        if (model_type == 1) {
          target += normal_lpdf(Y_miss[i, j] | nu[j], sigma[j]);
        } else {
          target += normal_lpdf(Y_miss[i, j] | 0, sigma[j]);
        }
      }
    }
  }
}
