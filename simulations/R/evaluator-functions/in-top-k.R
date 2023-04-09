################################################################################
## Top K Evaluator
################################################################################

top_k_fun <- function(fit_results, true_params_df, k) {

  # define the largest parameters
  top_modifiers <- true_params_df %>%
    dplyr::arrange(desc(abs(truth))) %>%
    dplyr::slice_head(n = k) %>%
    dplyr::pull(modifier)

  init_group_vars <- c(".dgp_name", ".method_name", "n", ".rep")
  group_vars <- c(".dgp_name", ".method_name", "n")
  eval_out <- fit_results %>%
    tidyr::unnest(cols = c("modifier", "estimate")) %>%
    dplyr::select(dplyr::all_of(c(init_group_vars, "modifier", "estimate"))) %>%
    dplyr::left_join(true_params_df, by = "modifier") %>%
    dplyr::group_by(dplyr::across({{init_group_vars}})) %>%
    dplyr::arrange(dplyr::desc(abs(estimate))) %>%
    dplyr::select(-"estimate") %>%
    dplyr::slice_head(n = k) %>%
    dplyr::mutate(in_top_k = (modifier %in% top_modifiers)) %>%
    dplyr::group_by(dplyr::across({{group_vars}})) %>%
    dplyr::summarise(prop_top_k = mean(in_top_k), .groups = "drop")

  return(eval_out)
}
