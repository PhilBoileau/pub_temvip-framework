################################################################################
# Observational Study, Continuous Outcome
################################################################################

# set the warnings to print as they occur
options(warn = 1)

# load required libraries
library(here)
library(data.table)
library(dplyr)
library(tidyr)
library(tibble)
library(simChef)
library(future)

# set up parallelization
plan(multisession, workers = 20L)

# define the data-generating process object
source(here("R/dgp-functions/obs-cont-outcome.R"))
dgp <- create_dgp(.dgp_fun = lm_tem_fun)

## define the method objects
source(here("R/method-functions/unihtee-wrapper.R"))
source(here("R/method-functions/mod-cov-wrapper.R"))
source(here("R/method-functions/aug-mod-cov-wrapper.R"))
tml_estimator <- create_method(
  .method_fun = unihtee_cont_rd_fun, estimator = "tmle"
)
onestep_estimator <- create_method(
  .method_fun = unihtee_cont_rd_fun, estimator = "onestep"
)
cf_tml_estimator <- create_method(
  .method_fun = unihtee_cont_rd_fun, estimator = "tmle", cross_fit = TRUE
)
cf_onestep_estimator <- create_method(
  .method_fun = unihtee_cont_rd_fun, estimator = "onestep", cross_fit = TRUE
)
mod_cov_estimator <- create_method(.method_fun = mod_cov_cont_rd_fun)
aug_mod_cov_estimator <- create_method(.method_fun = aug_mod_cov_cont_rd_fun)

## load the true parameter values
truth_df <- readRDS(here("data/parameter-values/lm-tem-dgp.Rds"))

## define the evaluator objects
source(here("R/evaluator-functions/bias.R"))
source(here("R/evaluator-functions/variance.R"))
source(here("R/evaluator-functions/coverage.R"))
source(here("R/evaluator-functions/fdr.R"))
source(here("R/evaluator-functions/tnr.R"))
source(here("R/evaluator-functions/tpr.R"))
bias_eval <- create_evaluator(.eval_fun = bias_fun, true_params_df = truth_df)
variance_eval <- create_evaluator(.eval_fun = variance_fun)
coverage_eval <- create_evaluator(
  .eval_fun = coverage_fun, true_params_df = truth_df
)
fdr_eval <- create_evaluator(.eval_fun = fdr_fun, true_params_df = truth_df)
tnr_eval <- create_evaluator(.eval_fun = tnr_fun, true_params_df = truth_df)
tpr_eval <- create_evaluator(.eval_fun = tpr_fun, true_params_df = truth_df)

## meal prep
experiment <- create_experiment(name = "observational-continuous") %>%
  add_dgp(dgp, name = "Observational Study, Continuous Outcome") %>%
  add_vary_across(
    .dgp = "Observational Study, Continuous Outcome",
    n = c(125, 250, 500, 1000, 2000)
  ) %>%
  add_method(tml_estimator, name = "TMLE" ) %>%
  add_method(onestep_estimator, name = "One Step") %>%
  ## add_method(cf_tml_estimator, name = "CF TMLE" ) %>%
  ## add_method(cf_onestep_estimator, name = "CF One Step") %>%
  add_method(mod_cov_estimator, name = "Mod. Cov.") %>%
  add_method(aug_mod_cov_estimator, name = "Aug. Mod. Cov.") %>%
  add_evaluator(bias_eval, name = "Empirical Bias") %>%
  add_evaluator(variance_eval, name = "Empirical Variance") %>%
  add_evaluator(coverage_eval, name = "Empirical Coverage") %>%
  add_evaluator(fdr_eval, name = "Empirical FDR") %>%
  add_evaluator(tnr_eval, name = "Empirical TNR") %>%
  add_evaluator(tpr_eval, name = "Empirical TPR")


## put it in the oven
set.seed(514)
results <- experiment$run(n_reps = 200, save = TRUE, verbose = 2)
