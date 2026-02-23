# Test script for get_befa_inits
source("R/utils_estimation.R")
source("R/utils_input.R") # For is_uni_model if needed, though get_befa_inits implements logic inline now.

# Mock data
stan_data <- list(J = 6, M = 2, m_obs = rep(0, 6), S_obs = diag(6))

# Helper to check dims
check_dims <- function(inits, param, expected_dims) {
  val <- inits[[1]][[param]]
  if (is.null(val)) stop(paste("Param", param, "is NULL"))

  dims <- if (is.vector(val)) length(val) else dim(val)

  if (length(dims) == 1 && length(expected_dims) == 1) {
    if (dims != expected_dims) stop(sprintf("%s: Expected %s, got %s", param, expected_dims, dims))
  } else {
    if (!all(dims == expected_dims)) stop(sprintf("%s: Expected (%s), got (%s)", param, paste(expected_dims, collapse = ","), paste(dims, collapse = ",")))
  }
}

cat("Testing NULL model...\n")
inits_null <- get_befa_inits("null", "unit_vector", stan_data, 1)
check_dims(inits_null, "nu", 0)
check_dims(inits_null, "sigma", 6)
check_dims(inits_null, "Z", c(0, 0)) # My current code does 0x0

cat("Testing RAW + UV model...\n")
inits_raw <- get_befa_inits("raw", "unit_vector", stan_data, 1)
check_dims(inits_raw, "nu", 6)
check_dims(inits_raw, "sigma", 6)
check_dims(inits_raw, "Z", c(6, 2))
check_dims(inits_raw, "Lambda_uni", c(0, 0))

cat("Testing COR + Normal model...\n")
inits_cor <- get_befa_inits("cor", "normal", stan_data, 1)
check_dims(inits_cor, "nu", 0)
check_dims(inits_cor, "sigma", 0)
check_dims(inits_cor, "Lambda_norm", c(6, 2))
check_dims(inits_cor, "Z", c(0, 0))

cat("All tests passed!\n")
