################################################################################
# True Negative Rate Function Evaluator
################################################################################

tnr_fun <- function(fit_results, true_params_df) {

  group_vars_tnp <- c(".rep", ".dgp_name", ".method_name", "n")
  group_vars_tnr <- c(".dgp_name", ".method_name", "n")
  eval_out <- fit_results %>%
    tidyr::unnest(cols = c("modifier", "estimate", "p_value_fdr")) %>%
    dplyr::left_join(true_params_df, by = "modifier") %>%
    dplyr::mutate(
      pred_tem = if_else(
        is.na(p_value_fdr),
        if_else(estimate == 0, 0, 1),     # for mod cov methods
        if_else(p_value_fdr < 0.05, 1, 0) # for unihtee methods
       ),
      true_tem = if_else(abs(truth) < 0.2, 0, 1)
    ) %>%
    dplyr::filter(true_tem == 0) %>%
    dplyr::group_by(dplyr::across({{group_vars_tnp}})) %>%
    dplyr::summarise(tnp = mean(1 - pred_tem), .groups = "drop") %>%
    dplyr::group_by(dplyr::across({{group_vars_tnr}})) %>%
    dplyr::summarize(tnr = mean(tnp), .groups = "drop")

  return(eval_out)
}
