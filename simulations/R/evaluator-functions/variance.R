################################################################################
## Variance Evaluator
################################################################################

variance_fun <- function(fit_results) {
  group_vars <- c(".dgp_name", ".method_name", "n", "modifier")
  eval_out <- fit_results %>%
    dplyr::filter(!(.method_name %in% c("Mod. Cov.", "Aug. Mod. Cov."))) %>%
    tidyr::unnest(cols = c("modifier", "estimate")) %>%
    dplyr::group_by(dplyr::across({{group_vars}})) %>%
    dplyr::summarize(variance = var(estimate), .groups = "drop")

  return(eval_out)
}
