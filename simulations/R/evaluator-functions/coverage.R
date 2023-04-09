################################################################################
## Coverage Evaluator
################################################################################

coverage_fun <- function(fit_results, true_params_df) {
  group_vars <- c(".dgp_name", ".method_name", "n", "modifier")
  eval_out <- fit_results %>%
    dplyr::filter(!(.method_name %in% c("Mod. Cov.", "Aug. Mod. Cov."))) %>%
    tidyr::unnest(cols = c("modifier", "ci_lower", "ci_upper")) %>%
    dplyr::left_join(true_params_df, by = "modifier") %>%
    dplyr::mutate(
      covered = dplyr::if_else(ci_lower < truth & truth < ci_upper, 1, 0)
    ) %>%
    dplyr::group_by(dplyr::across({{group_vars}})) %>%
    dplyr::summarize(coverage = mean(covered), .groups = "drop")

  return(eval_out)
}
