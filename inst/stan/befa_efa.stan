// Unified Bayesian Exploratory Factor Analysis Model
// Flags:
//   model_type:  1 = raw (means + covariances), 2 = cov, 3 = cor
//   lambda_type: 1 = unit_vector (UV), 2 = unidimensional (M=1 Beta-Scaled), 3 = normal

functions {
  // Sufficient statistics log-likelihood WITH means (raw model)
  real sufficient_ll_raw(int N, int J, vector m_obs, vector nu,
                         matrix S_obs, matrix Sigma) {
    vector[J] diff = m_obs - nu;
    real logdet = log_determinant(Sigma);
    real tr_part = trace(mdivide_right_spd(S_obs, Sigma));
    real mean_quad = dot_product(diff, mdivide_left_spd(Sigma, diff));
    return -0.5 * N * (logdet + tr_part + mean_quad);
  }

  // Sufficient statistics log-likelihood WITHOUT means (cov/cor model)
  real sufficient_ll_cov(int N, matrix S, matrix Sigma) {
    return -0.5 * N * (log_determinant(Sigma) + trace(mdivide_right_spd(S, Sigma)));
  }

  // Scaled Beta density for unidimensional model
  real scaled_beta_lpdf(vector y, real alpha, real beta) {
    int K = num_elements(y);
    return sum((alpha - 1) * log1p(y) + (beta - 1) * log1m(y))
               - K * (alpha + beta - 1) * log(2)
               - K * lbeta(alpha, beta);
  }

  // FIML with precomputed indices — raw model (with means)
  real incomplete_obs_raw_lpdf(matrix Y_miss,
                               array[] int n_obs_vec,
                               array[,] int obs_idx,
                               int N_inc,
                               vector nu,
                               matrix Sigma) {
    real ll = 0;
    for (i in 1:N_inc) {
      int ni = n_obs_vec[i];
      array[ni] int idx_i = obs_idx[i, 1:ni];
      vector[ni] y_i = to_vector(Y_miss[i, idx_i]);
      vector[ni] nu_i = nu[idx_i];
      matrix[ni, ni] Sigma_i = Sigma[idx_i, idx_i];
      ll += multi_normal_lpdf(y_i | nu_i, Sigma_i);
    }
    return ll;
  }

  // FIML with precomputed indices — cov/cor model (zero-mean)
  real incomplete_obs_cov_lpdf(matrix Y_miss,
                               array[] int n_obs_vec,
                               array[,] int obs_idx,
                               int N_inc,
                               matrix Sigma) {
    real ll = 0;
    for (i in 1:N_inc) {
      int ni = n_obs_vec[i];
      array[ni] int idx_i = obs_idx[i, 1:ni];
      vector[ni] y_i = to_vector(Y_miss[i, idx_i]);
      matrix[ni, ni] Sigma_i = Sigma[idx_i, idx_i];
      ll += multi_normal_lpdf(y_i | rep_vector(0, ni), Sigma_i);
    }
    return ll;
  }
}

data {
  // Dimensions
  int<lower=0> N_complete;
  int<lower=0> N_incomplete;
  int<lower=0> J;
  int<lower=0> M;

  // Model configuration flags
  int<lower=1, upper=3> model_type;   // 1=raw, 2=cov, 3=cor
  int<lower=1, upper=3> lambda_type;  // 1=UV, 2=uni, 3=normal

  // Sufficient statistics
  vector[J] m_obs;          // sample means
  matrix[J, J] S_obs;       // sample covariance
  matrix[J, J] R_obs;       // sample correlation

  // Incomplete data for FIML
  matrix[N_incomplete, J] Y_miss;
  int<lower=0> max_obs;
  array[N_incomplete] int<lower=0> n_obs;
  array[N_incomplete, max_obs] int<lower=0, upper=J> obs_idx;

  // Prior hyperparameters
  vector[2] pr_nu;              // Normal(mean, sd) for nu
  vector<lower=0>[3] pr_sigma;  // Student-t(df, loc, scale)
  real<lower=0> pr_xi;          // Xi concentration for UV
  vector<lower=0>[2] pr_h2;     // Beta(a,b) for communalities
  vector[2] pr_Lambda;          // Beta(a,b) for uni OR Normal(m,s) for norm
  vector<lower=0>[2] pr_psi;    // InvGamma(a,b) for uniquenesses
}

parameters {
  vector[model_type == 1 ? J : 0] nu;
  vector<lower=0>[model_type <= 2 ? J : 0] sigma;

  // Unit Vector (UV) parameterization
  vector<lower=0, upper=1>[lambda_type == 1 ? J : 0] h2;
  matrix[lambda_type == 1 ? J : 0, M] Z;

  // Unidimensional (Beta-Scaled) parameterization
  matrix<lower=-1, upper=1>[lambda_type == 2 ? J : 0, 1] Lambda_uni;

  // Normal parameterization
  matrix[lambda_type == 3 ? J : 0, M] Lambda_norm;
  vector<lower=0>[lambda_type == 3 ? J : 0] psi;
}

transformed parameters {
  corr_matrix[J] Rho;
  matrix[J, J] Sigma;
  matrix[J, M] Lambda;
  vector[J] Psi;

  // UV-specific intermediates
  vector[lambda_type == 1 ? J : 0] Z_norm;

  if (lambda_type == 1) {
    // Unit Vector
    Z_norm = sqrt(diagonal(tcrossprod(Z)));
    {
      matrix[J, M] Z_UV = diag_pre_multiply(inv(Z_norm), Z);
      Lambda = diag_pre_multiply(sqrt(h2), Z_UV);
      Psi = 1 - h2;
      Rho = add_diag(tcrossprod(Lambda), Psi);
    }
  } else if (lambda_type == 2) {
    // Unidimensional (Beta-Scaled)
    Lambda = Lambda_uni;
    {
      matrix[J, J] h2_mat = tcrossprod(Lambda);
      Psi = 1 - diagonal(h2_mat);
      Rho = h2_mat + diag_matrix(Psi);
    }
  } else {
    // Normal
    Lambda = Lambda_norm;
    Psi = psi;
    {
      matrix[J, J] Sigma_raw = add_diag(tcrossprod(Lambda), psi);
      vector[J] inv_sd = inv_sqrt(diagonal(Sigma_raw));
      Rho = quad_form_diag(Sigma_raw, inv_sd);
    }
  }

  // Model-implied covariance matrix
  if (model_type <= 2) {
    Sigma = quad_form_diag(Rho, sigma);
  } else {
    Sigma = Rho;
  }
}

model {

  // Prior distribution per model type
  if (model_type == 1) {
    target += normal_lpdf(nu | pr_nu[1], pr_nu[2]);
  }

  if (model_type <= 2) {
    target += student_t_lpdf(sigma | pr_sigma[1], pr_sigma[2], pr_sigma[3]);
    target += J * -student_t_lccdf(0 | pr_sigma[1], pr_sigma[2], pr_sigma[3]);
  }

  // Prior distribution per factor loadings prior specification
  if (lambda_type == 1) {
    // UV: Unit-vector prior distribution
    target += beta_lpdf(h2 | pr_h2[1], pr_h2[2]);
    target += -0.5 * dot_self(Z_norm) + pr_xi * sum(log(Z_norm));
  } else if (lambda_type == 2) {
    // Uni: Scaled Beta prior distribution
    target += scaled_beta_lpdf(to_vector(Lambda_uni) | pr_Lambda[1], pr_Lambda[2]);
  } else {
    // Normal: Normal prior on loadings + InvGamma on uniquenesses
    target += normal_lpdf(to_vector(Lambda_norm) | pr_Lambda[1], pr_Lambda[2]);
    target += inv_gamma_lpdf(psi | pr_psi[1], pr_psi[2]);
  }

  // Log-Likelihood: Complete observations
  if (N_complete > 0) {
    if (model_type == 1) {
      target += sufficient_ll_raw(N_complete, J, m_obs, nu, S_obs, Sigma);
    } else if (model_type == 2) {
      target += sufficient_ll_cov(N_complete, S_obs, Sigma);
    } else {
      target += sufficient_ll_cov(N_complete, R_obs, Rho);
    }
  }

  // Log-Likelihood: Incomplete observations (FIML)
  if (N_incomplete > 0) {
    if (model_type == 1) {
      target += incomplete_obs_raw_lpdf(Y_miss | n_obs, obs_idx, N_incomplete, nu, Sigma);
    } else if (model_type == 2) {
      target += incomplete_obs_cov_lpdf(Y_miss | n_obs, obs_idx, N_incomplete, Sigma);
    } else {
      target += incomplete_obs_cov_lpdf(Y_miss | n_obs, obs_idx, N_incomplete, Rho);
    }
  }
}
