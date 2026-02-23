# Load BayesEFA
library(BayesEFA)

# Define a true loading matrix with a clear simple structure
L_true <- matrix(c(
  0.7, 0.0,
  0.7, 0.0,
  0.7, 0.0,
  0.0, 0.7,
  0.0, 0.7,
  0.0, 0.7
), ncol = 2, byrow = TRUE)

# Generate the corresponding correlation matrix
R_true <- tcrossprod(L_true)
diag(R_true) <- 1

# Simulate data for 300 observations
set.seed(2026)
X <- matrix(rnorm(300 * nrow(L_true)), ncol = nrow(L_true)) %*% chol(R_true)
colnames(X) <- paste0("Item_", 1:6)

# Fit the model without automatic post-processing
# We expect this to fail with current codebase
tryCatch(
  {
    fit_unaligned <- befa(
      data = X,
      n_factors = 2,
      compute_fit_indices = FALSE,
      compute_reliability = FALSE,
      rotate = "none",
      seed = 123,
      refresh = 0,
      show_message = FALSE,
      iter = 100,
      chains = 1
    )
    print("SUCCESS: Model ran without error.")
  },
  error = function(e) {
    print(paste("FAILURE:", e$message))
  }
)
