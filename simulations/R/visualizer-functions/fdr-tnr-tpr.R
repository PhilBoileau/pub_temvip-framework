################################################################################
# FDR, TNR, TPR plotting function
################################################################################

summary_plot_fun <- function(eval_results) {

  # combine the evaluator results into a single tibble
  emp_fdr_tbl <- eval_results$`Empirical FDR` %>%
    dplyr::mutate(
      metric = "FDR",
      value = fdr,
      yintercept = 0.05
    ) %>%
    dplyr::select(-fdr)
  emp_tpr_tbl <- eval_results$`Empirical TPR` %>%
    dplyr::mutate(
      metric = "TPR",
      value = tpr,
      yintercept = 1
    ) %>%
    dplyr::select(-tpr)
  emp_tnr_tbl <- eval_results$`Empirical TNR` %>%
    dplyr::mutate(
      metric = "TNR",
      value = tnr,
      yintercept = 1
    ) %>%
    dplyr::select(-tnr)
  comb_tbl <- dplyr::bind_rows(emp_fdr_tbl, emp_tpr_tbl, emp_tnr_tbl)

  # order the methods
  comb_tbl$.method_name <- factor(
    comb_tbl$.method_name,
    levels = c(
      "Mod. Cov.", "Aug. Mod. Cov.",
      "One Step", "CF One Step",
      "TMLE", "CF TMLE"
    )
  )

  plt <- ggplot2::ggplot(comb_tbl) +
    ggplot2::aes(x = n, y = value, colour = .method_name) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::geom_line(alpha = 0.7) +
    ggplot2::geom_hline(aes(yintercept = yintercept), colour = "red",
                        linetype = 2, alpha = 0.5) +
    ggplot2::facet_grid(rows = ggplot2::vars(metric)) +
    ggplot2::scale_color_discrete(name = "Method") +
    ggplot2::xlab("Sample Size") +
    ggplot2::ylab("Value") +
    ggplot2::theme_bw()

  return(plt)
}
