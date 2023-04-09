################################################################################
## Bias Evaluator
################################################################################

bias_fun <- function(fit_results, true_params_df) {
  group_vars <- c(".dgp_name", ".method_name", "n", "modifier")
  eval_out <- fit_results %>%
    dplyr::filter(!(.method_name %in% c("Mod. Cov.", "Aug. Mod. Cov."))) %>%
    tidyr::unnest(cols = c("modifier", "estimate")) %>%
    dplyr::left_join(true_params_df, by = "modifier") %>%
    dplyr::group_by(dplyr::across({{group_vars}})) %>%
    dplyr::summarize(bias = mean(estimate - truth), .groups = "drop")

  return(eval_out)
}
