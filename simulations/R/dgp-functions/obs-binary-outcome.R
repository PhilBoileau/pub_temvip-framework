################################################################################
## Linear Model with Treatment Effect Modification, Binary Outcome
################################################################################

## general toeplitz matrix
toeplitz_cov_mat <- function(p, rho, alpha) {
  times <- seq_len(p)
  H <- abs(outer(times, times, "-")) + diag(p)
  H <- H^-(1 + alpha) * rho
  covmat <- H + diag(p) * (1 - rho)
  return(covmat)
}

## Input:
##   n, an integer representing the number of observations to simulate
##   cov_mat, a 100x100 covariance matrix of the covariates
## Output:
##   A list object containing 100 covariates, a binary treatment indicator,
##   and the potential outcomes of n independently simulated observations in an
##   observational study.
binary_tem_fun <- function(
  n = 125,
  cov_mat = toeplitz_cov_mat(100, 0.1, 0.8),
  simchef = TRUE
) {

  ## generate the biomarkers
  W <- MASS::mvrnorm(n = n, mu = rep(0, 100), Sigma = cov_mat / 4)
  colnames(W) <- paste0("W", seq_len(100))

  ## simulate the binary treatment assignment
  prop_score <- plogis(0.25 * (W[, 1] + W[, 2] + W[, 3]))
  A <- rbinom(n, 1, prop_score)

  ## simulate the outcomes under treatment and controls

  ## define the sparse coefficients vector for the main effects, and controls
  ## and treatment interactions
  beta_0 <- c(rep(-0.5, 5), rep(0, 95))
  beta_1 <- c(rep(0.5, 5), rep(0, 95))

  ## generate the potential outcomes
  W_t <- t(W)
  q_Y_0 <- plogis(as.vector(-1 + crossprod(W_t, beta_0)))
  q_Y_1 <- plogis(as.vector(1 + crossprod(W_t, beta_1)))

  ## define the observed outcome
  Y_0 <- rbinom(n, 1, q_Y_0)
  Y_1 <- rbinom(n, 1, q_Y_1)
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
      "q_Y_0" = q_Y_0,
      "q_Y_1" = q_Y_1,
      "prop_score" = prop_score,
      "A" = A,
      "W" = W
    )
  }

  return(sample_ls)
}

## compute the true parameter values using Monte Carlo
get_binary_tem_truth <- function(
  n = 1e6,
  cov_mat = toeplitz_cov_mat(100, 0.1, 0.8)
) {

  ## generate the full data
  full_data_ls <- binary_tem_fun(n, cov_mat, simchef = FALSE)

  ## compute the risk-difference-scale TEM-VIP parameters
  log_pot_out_diff <- log(full_data_ls$q_Y_1) - log(full_data_ls$q_Y_0)
  res_vec <- apply(
    full_data_ls$W, 2, function(mod) cov(log_pot_out_diff, mod) / var(mod)
  )

  tibble(
    "modifier" = names(res_vec),
    "truth" = res_vec
  )
}

## ## save the true parameter values
## library(here)
## set.seed(917498)
## saveRDS(
##   get_binary_tem_truth(),
##   here("data/parameter-values/binary-tem-dgp.Rds")
## )
