################################################################################
## unihtee wrappers
################################################################################

library(unihtee)
library(glmnet)
library(ranger)
library(xgboost)
library(earth)
library(data.table)
library(sl3)

## For observational continuous outcome simulations
unihtee_cont_rd_fun <- function(Y, A, W, estimator, cross_fit = FALSE) {

  ## orgnanize the data into a data.table
  dt <- data.table::data.table(Y, A, W)
  confounders <- paste0("W", seq_len(500))
  colnames(dt) <- c("Y", "A", confounders)

  ## define the super learner to use## define base learners
  ## define the interactions
  interactions <- lapply(confounders, function(w) c(w, "A"))
  lrnr_interactions <- sl3::Lrnr_define_interactions$new(interactions)
  lrnr_mean <- sl3::Lrnr_mean$new()
  lrnr_lasso_co <- sl3::make_learner(
    sl3::Pipeline, lrnr_interactions, sl3::Lrnr_glmnet$new()
  )
  lrnr_enet_co <- sl3::make_learner(
    sl3::Pipeline, lrnr_interactions, sl3::Lrnr_glmnet$new(alpha = 0.5)
  )
  lrnr_ridge_co <- sl3::make_learner(
    sl3::Pipeline, lrnr_interactions, sl3::Lrnr_glmnet$new(alpha = 0)
  )
  lrnr_earth <- sl3:::Lrnr_earth$new()
  lrnr_earth_interactions <- sl3::make_learner(
    sl3::Pipeline, lrnr_interactions, sl3::Lrnr_earth$new()
  )
  lrnr_lasso_ps <- sl3::Lrnr_glmnet$new()
  lrnr_enet_ps <- sl3::Lrnr_glmnet$new(alpha = 0.5)
  lrnr_ridge_ps <- sl3::Lrnr_glmnet$new(alpha = 0)
  lrnr_ranger <- sl3::Lrnr_ranger$new()

  ## define the super learners
  lrnr_sl_ps <- sl3::Lrnr_sl$new(
    learners = list(
      lrnr_mean, lrnr_lasso_ps, lrnr_enet_ps, lrnr_ridge_ps, lrnr_earth,
      lrnr_ranger
    ),
    metalearner = sl3::make_learner(
      sl3::Lrnr_solnp, sl3::metalearner_logistic_binomial,
      sl3::loss_squared_error
    )
  )
  lrnr_sl_co <- sl3::Lrnr_sl$new(
    learners = list(
      lrnr_mean, lrnr_lasso_co, lrnr_enet_co, lrnr_ridge_co,
      lrnr_earth_interactions, lrnr_ranger
    ),
    metalearner = sl3::make_learner(
      sl3::Lrnr_solnp, sl3::metalearner_linear, sl3::loss_squared_error
    )
  )

  ## apply unicate for RD TEM VIP to continuous outcome
  results <- unihtee::unihtee(
    data = dt,
    confounders = confounders,
    modifiers = confounders,
    exposure = "A",
    outcome = "Y",
    outcome_type = "continuous",
    effect = "absolute",
    cross_fit = cross_fit,
    estimator = estimator,
    prop_score_estimator = lrnr_sl_ps,
    cond_outcome_estimator = lrnr_sl_co
  )

  return(results)

}

## for observation binary outcome simulations
unihtee_bin_rr_fun <- function(Y, A, W, estimator, cross_fit = FALSE) {

  ## orgnanize the data into a data.table
  dt <- data.table::data.table(Y, A, W)
  confounders <- paste0("W", seq_len(100))
  colnames(dt) <- c("Y", "A", confounders)

  ## define the super learner to use## define base learners
  ## define the interactions
  interactions <- lapply(confounders, function(w) c(w, "A"))
  lrnr_interactions <- sl3::Lrnr_define_interactions$new(interactions)
  lrnr_lasso_co <- sl3::make_learner(
    sl3::Pipeline, lrnr_interactions, sl3::Lrnr_glmnet$new(family = "binomial")
  )
  lrnr_enet_co <- sl3::make_learner(
    sl3::Pipeline, lrnr_interactions,
    sl3::Lrnr_glmnet$new(family = "binomial", alpha = 0.5)
  )
  lrnr_ridge_co <- sl3::make_learner(
    sl3::Pipeline, lrnr_interactions,
    sl3::Lrnr_glmnet$new(family = "binomial", alpha = 0)
  )
  lrnr_lasso_ps <- sl3::Lrnr_glmnet$new(family = "binomial")
  lrnr_enet_ps <- sl3::Lrnr_glmnet$new(family = "binomial", alpha = 0.5)
  lrnr_ridge_ps <- sl3::Lrnr_glmnet$new(family = "binomial", alpha = 0)
  lrnr_earth_ps <- sl3:::Lrnr_earth$new(glm = list(family = "binomial"))
  lrnr_earth_co <- sl3::make_learner(
    sl3::Pipeline, lrnr_interactions,
    sl3::Lrnr_earth$new(
      glm = list(family = "binomial")
    )
  )
  lrnr_ranger <- sl3::Lrnr_ranger$new()

 ## define the super learners
  bin_metalearner <- sl3::make_learner(
    sl3::Lrnr_solnp, sl3::metalearner_logistic_binomial,
    sl3::loss_squared_error
  )
  lrnr_sl_ps <- sl3::Lrnr_sl$new(
    learners = list(
      lrnr_lasso_ps, lrnr_enet_ps, lrnr_ridge_ps, lrnr_earth_ps, lrnr_ranger
    ),
    metalearner = bin_metalearner
  )
  lrnr_sl_co <- sl3::Lrnr_sl$new(
    learners = list(
      lrnr_lasso_co, lrnr_enet_co, lrnr_ridge_co, lrnr_earth_co, lrnr_ranger
    ),
    metalearner = bin_metalearner
  )

  ## apply unicate for RR TEM VIP to binary outcome
  results <- unihtee::unihtee(
    data = dt,
    confounders = confounders,
    modifiers = confounders,
    exposure = "A",
    outcome = "Y",
    outcome_type = "binary",
    effect = "relative",
    cross_fit = cross_fit,
    estimator = estimator,
    prop_score_estimator = lrnr_sl_ps,
    cond_outcome_estimator = lrnr_sl_co
  )

  return(results)

}

## for observation survival outcome simulations
unihtee_tte_rd_fun <- function(
  T, F, A, W, estimator, cross_fit = FALSE
) {

  ## orgnanize the data into a data.table
  dt <- data.table::data.table(T, (1 - F), A, W)
  confounders <- paste0("W", seq_len(300))
  colnames(dt) <- c("T", "C", "A", confounders)
  dt$prop_scores <- 0.5

  ## apply unicate for RD TEM VIP to time-to-event outcome
  form_str <- paste(
    "~ A +",
    paste0("W", seq_len(300), collapse = " + "), "+",
    paste0("A:W", seq_len(300), collapse = " + ")
  )
  results <- unihtee::unihtee(
    data = dt,
    confounders = confounders,
    modifiers = confounders,
    exposure = "A",
    outcome = "T",
    censoring = "C",
    time_cutoff = 9,
    outcome_type = "time-to-event",
    effect = "absolute",
    estimator = estimator,
    cross_fit = cross_fit,
    prop_score_values = "prop_scores",
    failure_hazard_estimator = sl3::Lrnr_glmnet$new(
      formula = form_str, family = "binomial"
    ),
    censoring_hazard_estimator = sl3::Lrnr_glmnet$new(family = "binomial")
  )

  return(results)

}
