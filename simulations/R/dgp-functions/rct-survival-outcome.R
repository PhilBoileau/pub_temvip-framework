################################################################################
## RCT with Treatment Effect Modification, Survival Outcome
################################################################################

library(data.table)

## define a block diagonal covariance matrix generator
block_diag_cov_mat <- function(p) {


  ## generate p %% 10 blocks
  sub_mat_ls <- lapply(
    seq_len(round(p /  10)),
    function(idx) {
      cov2cor(
        clusterGeneration::genPositiveDefMat(
        10, "onion", rangeVar = c(1, 1.1)
        )$Sigma
      )
    }
  )

  return(as.matrix(Matrix::bdiag(sub_mat_ls)))
}

## define the conditional survival hazard
cond_survival_hazard <- function(
  time, exposure, W_1, W_2, W_3, W_4, W_5, W_6, W_7, W_8, W_9, W_10
) {
  (time <= 9) * plogis(
    -(2 + exposure) + 5 * (2 * exposure - 1) *
      (W_1 + W_2 + W_3 + W_4 + W_5 + W_6 + W_7 + W_8 + W_9 + W_10)
  ) +
  (time == 10)
}

## define the conditional censoring hazard
cond_censoring_hazard <- function(exposure, W_1) {
  plogis(-(5 + exposure + W_1))
}

## Input:
##   n, an integer representing the number of observations to simulate
##   cov_mat, a 300x300 covariance matrix of the covariates
## Output:
##   A list object containing 300 covariates, a binary treatment indicator,
##   survival time and right-censoring indicator of data generated according an
##   RCT.
survival_tem_fun <- function(
  n = 125,
  cov_mat = thresholding_cov_mat(300)
) {

  ## define the conditional survival hazard
  cond_survival_hazard <- function(
    time, exposure, W_1, W_2, W_3, W_4, W_5, W_6, W_7, W_8, W_9, W_10
  ) {
    (time <= 9) * plogis(
      -(2 + exposure) + 5 * (2 * exposure - 1) *
        (W_1 + W_2 + W_3 + W_4 + W_5 + W_6 + W_7 + W_8 + W_9 + W_10)
    ) +
      (time == 10)
  }

  ## define the conditional censoring hazard
  cond_censoring_hazard <- function(exposure, W_1) {
    plogis(-(2 + exposure + W_1))
  }

  ## generate the biomarkers
  W <- MASS::mvrnorm(n = n, mu = rep(0, 300), Sigma = cov_mat / 2)
  colnames(W) <- paste0("W", seq_len(300))

  ## simulate the binary treatment assignment using 1:1 treatment/control arms
  A <- rbinom(n, 1, 0.5)

  ## simulate the outcomes under treatment and controls
  failure_time_sim <- function(A, W) {
    sapply(
      seq_len(n),
      function(obs) {
        failure_time <- NA
        for (t in 1:10) {
          prob <- cond_survival_hazard(
            t, A, W[obs, 1], W[obs, 2], W[obs, 3], W[obs, 4], W[obs, 5],
            W[obs, 6], W[obs, 7], W[obs, 8], W[obs, 9], W[obs, 10]
          )
          status <- rbinom(1, 1, prob)
          if (status == 1) {
            failure_time <- t
            break
          }
        }
        return(failure_time)
      }
    )
  }
  failure_time_1 <- failure_time_sim(1, W)
  failure_time_0 <- failure_time_sim(0, W)


  # generate the censoring events for t = 1 to 10
  censor_time_sim <- function(A, W) {
    sapply(
      seq_len(n),
      function(obs) {
        censor_time <- NA
        for (t in 1:10) {
          prob <- cond_censoring_hazard(A, W[obs, 1])
          status <- rbinom(1, 1, prob)
          if (status == 1) {
            censor_time <- t
            break
          }
        }
        if (is.na(censor_time)) censor_time <- 11
        return(censor_time)
      }
    )
  }
  censor_time_1 <- censor_time_sim(1, W)
  censor_time_0 <- censor_time_sim(0, W)

  # compile the failure and censoring times
  failure_time <- sapply(
    seq_len(n),
    function(obs) {
      if (A[obs] == 1) failure_time_1[obs] else failure_time_0[obs]
    }
  )
  censor_time <- sapply(
    seq_len(n),
    function(obs) {
      if (A[obs] == 1) censor_time_1[obs] else censor_time_0[obs]
    }
  )

  # assess the observed time-to-event and censoring indicator
  time <- sapply(
    seq_len(n),
    function(obs) {
      if (censor_time[obs] < failure_time[obs]) {
        censor_time[obs]
      } else {
        failure_time[obs]
      }
    }
  )
  failure <- sapply(
    seq_len(n),
    function(obs) if (time[obs] == censor_time[obs]) 0 else 1
  )
  ## assembled into a tibble
  sample_ls <- list(
    "T" = time,
    "F" = failure,
    "A" = A,
    "W" = W
  )

  return(sample_ls)
}

## compute the true parameter values using Monte Carlo
get_survival_tem_truth <- function(n = 1e5, cov_mat, t_0) {

  ## generate the full data
  full_data_ls <- survival_tem_fun(n, cov_mat)
  full_dt <- as.data.table(full_data_ls)

  ## compute the restricted mean survival times truncated at time t
  long_dt <- lapply(
    seq_len(nrow(full_dt)),
    function(obs) {
      obs_dt <- full_dt[obs]
      obs_dt <- obs_dt[rep(1:.N, t_0)]
      obs_dt$time <- seq_len(t_0)
      obs_dt
    }
  )
  long_dt <- rbindlist(long_dt, idcol = "id")

  # compute the true failure hazard at each timepoint under each condition
  exp_truth <- cond_survival_hazard(
    long_dt$time, 1,
    long_dt$W.W1, long_dt$W.W2, long_dt$W.W3, long_dt$W.W4, long_dt$W.W5,
    long_dt$W.W6, long_dt$W.W7, long_dt$W.W8, long_dt$W.W9, long_dt$W.W10
  )
  noexp_truth <- cond_survival_hazard(
    long_dt$time, 0,
    long_dt$W.W1, long_dt$W.W2, long_dt$W.W3, long_dt$W.W4, long_dt$W.W5,
    long_dt$W.W6, long_dt$W.W7, long_dt$W.W8, long_dt$W.W9, long_dt$W.W10
  )

  # compute the survival probability at each time
  long_dt$true_haz_exp <- exp_truth
  long_dt$true_haz_noexp <- noexp_truth
  long_dt[, surv_exp := cumprod(1 - true_haz_exp), by = "id"]
  long_dt[, surv_noexp := cumprod(1 - true_haz_noexp), by = "id"]

  # compute the restricted mean survival time differences
  long_dt[, rmst := cumsum(surv_exp - surv_noexp), by = "id"]

  # retain only the rmst differences at time_cutoff, and compute parameters
  res_dt <- long_dt[time == t_0]

  ## compute the risk-difference-scale TEM-VIP parameters
  rmst <- as.vector(res_dt$rmst)
  params <- sapply(
    seq_len(300) + 4,
    function(mod_idx) {
      mod <- res_dt[[mod_idx]]
      cov(mod, rmst) / var(mod)
    }
  )
  names(params) <- paste0("W", seq_len(300))

  return(list(
    params = dplyr::tibble(
      "modifier" = names(params),
      "truth" = params
    ),
    surv_exp = res_dt$surv_exp,
    surv_noexp = res_dt$surv_noexp
  ))
}

## save the true parameter values
library(here)
set.seed(61234)
cov_mat <- block_diag_cov_mat(300)
saveRDS(
  get_survival_tem_truth(1e5, cov_mat = cov_mat, t = 9),
  here("data/parameter-values/survival-tem-dgp.Rds")
)
