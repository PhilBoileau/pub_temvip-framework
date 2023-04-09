################################################################################
## Nonlinear Model with Treatment Effect Modification, Continuous Outcome
################################################################################

## Input:
##   n, an integer representing the number of observations to simulate
##   cov_mat, a 500x500 covariance matrix of the covariates
## Output:
##   A list object containing 500 covariates, a binary treatment indicator,
##   and the potential outcomes of n independently simulated observations in a
##   perfect RCT.
lm_tem_fun <- function(
  n = 125,
  cov_mat = diag(1, nrow = 500),
  simchef = TRUE
) {

  ## generate the biomarkers
  W <- MASS::mvrnorm(n = n, mu = rep(0, 500), Sigma = cov_mat)
  colnames(W) <- paste0("W", seq_len(500))

  ## simulate the binary treatment assignment
  prop_score <- plogis(0.25 * (W[, 1] - W[, 2] + W[, 3]))
  A <- rbinom(n, 1, prop_score)

  ## simulate the outcomes under treatment and controls

  ## define the sparse coefficients vector for the main effects, and controls
  ## and treatment interactions
  beta_main <- c(rep(2, 5), rep(0, 495))
  beta_0 <- c(rep(-2, 5), rep(0, 495))
  beta_1 <- c(rep(3, 5), rep(0, 495))

  ## generate the potential outcomes
  epsilon <- rnorm(n = n, mean = 0, sd = 0.5)
  W_t <- t(W)
  main_effects <- rowSums(1 + abs(t(W_t * beta_main)))
  Y_0 <- as.vector(main_effects + crossprod(W_t, beta_0) + epsilon)
  Y_1 <- as.vector(main_effects + crossprod(W_t, beta_1) + epsilon)

  ## define the observed outcome
  Y <- ifelse(A == 0, Y_0, Y_1)

  ## assembled into a tibble
  if (simchef) {
    sample_ls <- list(
      "Y" = Y,
      "A" = A,
      "W" = W
    )
  } else {
    sample_ls <- list(
      "Y" = Y,
      "Y_0" = Y_0,
      "Y_1" = Y_1,
      "prop_score" = prop_score,
      "A" = A,
      "W" = W
    )
  }

  return(sample_ls)
}

## compute the true parameter values using Monte Carlo
get_lm_tem_truth <- function(n = 1e5, cov_mat = diag(1, nrow = 500)) {

  ## generate the full data
  full_data_ls <- lm_tem_fun(n, cov_mat, simchef = FALSE)

  ## compute the risk-difference-scale TEM-VIP parameters
  pot_out_diff <- full_data_ls$Y_1 - full_data_ls$Y_0
  res_vec <- apply(
    full_data_ls$W, 2, function(mod) cov(pot_out_diff, mod) / var(mod)
  )

  ## return a tibble
  tibble(
    "modifier" = names(res_vec),
    "truth" = res_vec
  )
}

## save the true parameter values
library(here)
set.seed(11346512)
saveRDS(get_lm_tem_truth(), here("data/parameter-values/lm-tem-dgp.Rds"))
