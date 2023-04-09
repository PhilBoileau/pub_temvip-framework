################################################################################
## Agumented modified covariates wrappers
################################################################################

library(personalized)
library(data.table)

## For observational continuous outcome simulations
aug_mod_cov_cont_rd_fun <- function(Y, A, W) {

  ## estimate the propensity score
  prop_func <- personalized::create.propensity.function(
    crossfit = TRUE,
    nfolds.crossfit = 5,
    cv.glmnet.args = list(type.measure = "auc", nfolds = 10)
  )

  ## specify the augmentation function
  aug_func <- personalized::create.augmentation.function(
    family = "gaussian",
    crossfit = TRUE,
    nfolds.crossfit = 5,
    cv.glmnet.args = list(type.measure = "mae", nfolds = 10)
  )

  ## estimate the conditional average treatment effect
  aug_mod_cov <- personalized::fit.subgroup(
    x = W,
    trt = A,
    y = Y,
    propensity.func = prop_func,
    augment.func = aug_func,
    loss = "sq_loss_lasso",
    nfolds = 10
  )

  ## massage the output to look like unihtee's
  data.table::data.table(
    "modifier" = rownames(aug_mod_cov$coefficients)[-c(1, 2)],
    "estimate" = aug_mod_cov$coefficients[-c(1, 2)],
    "se" = NA,
    "z" = NA,
    "p_value" = NA,
    "ci_lower" = NA,
    "ci_upper" = NA,
    "p_value_fdr" = NA
  )

}

## For observational continuous outcome simulations
aug_mod_cov_bin_rr_fun <- function(Y, A, W) {

  ## estimate the propensity score
  prop_func <- personalized::create.propensity.function(
    crossfit = TRUE,
    nfolds.crossfit = 5,
    cv.glmnet.args = list(type.measure = "auc", nfolds = 10)
  )

  ## specify the augmentation function
  aug_func <- personalized::create.augmentation.function(
    family = "binomial",
    crossfit = TRUE,
    nfolds.crossfit = 5,
    cv.glmnet.args = list(type.measure = "mae", nfolds = 10)
  )

  ## estimate the conditional average treatment effect
  aug_mod_cov <- personalized::fit.subgroup(
    x = W,
    trt = A,
    y = Y,
    propensity.func = prop_func,
    augment.func = aug_func,
    loss = "logistic_loss_lasso",
    nfolds = 10
  )

  ## massage the output to look like unihtee's
  data.table::data.table(
    "modifier" = rownames(aug_mod_cov$coefficients)[-c(1, 2)],
    "estimate" = aug_mod_cov$coefficients[-c(1, 2)],
    "se" = NA,
    "z" = NA,
    "p_value" = NA,
    "ci_lower" = NA,
    "ci_upper" = NA,
    "p_value_fdr" = NA
  )

}

## For RCT time-to-event outcome simulations
aug_mod_cov_tte_rd_fun <- function(T, F, A, W) {

  ## defin the propensity score
  prop_func <- function(x, trt) 0.5

  ## specify the augmentation function
  aug_func <- personalized::create.augmentation.function(
    family = "cox",
    crossfit = TRUE,
    nfolds.crossfit = 5,
    cv.glmnet.args = list(type.measure = "C", nfolds = 10)
  )

  ## estimate the conditional average treatment effect
  aug_mod_cov <- personalized::fit.subgroup(
    x = W,
    trt = A,
    y = survival::Surv(T, F),
    propensity.func = prop_func,
    augment.func = aug_func,
    loss = "cox_loss_lasso",
    nfolds = 10
  )

  ## massage the output to look like unihtee's
  data.table::data.table(
    "modifier" = rownames(aug_mod_cov$coefficients)[-1],
    "estimate" = aug_mod_cov$coefficients[-1],
    "se" = NA,
    "z" = NA,
    "p_value" = NA,
    "ci_lower" = NA,
    "ci_upper" = NA,
    "p_value_fdr" = NA
  )

}
