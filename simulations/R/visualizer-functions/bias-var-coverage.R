################################################################################
# Bias, variance, coverage plotting function
################################################################################

bias_var_coverage_fun <- function(eval_results, tem_indices) {

  ## combine the evaluator results into a single tibble
  emp_bias_tbl <- eval_results$`Empirical Bias` %>%
    dplyr::mutate(
      metric = "Bias",
      value = bias,
      yintercept = 0
    ) %>%
    dplyr::select(-bias)
  emp_variance_tbl <- eval_results$`Empirical Variance` %>%
    dplyr::mutate(
      metric = "Variance",
      value = variance,
      yintercept = 0
    ) %>%
    dplyr::select(-variance)
  comb_tbl <- dplyr::bind_rows(emp_bias_tbl, emp_variance_tbl)

  # order the methods
  comb_tbl$.method_name <- factor(
    comb_tbl$.method_name,
    levels = c("One Step", "CF One Step", "TMLE", "CF TMLE")
  )

  # order metric
  comb_tbl$metric <- factor(
    comb_tbl$metric,
    levels = c("Bias", "Variance")
  )

  # order the modifiers by index
  comb_tbl$modifier <- as.numeric(stringr::str_sub(comb_tbl$modifier, 2))

  ## create the plots
  tem_plt <- comb_plt %>%
    dplyr::filter(modifier %in% tem_indices) %>%
    ggplot2::ggplot() +
    ggplot2::aes(
      x = modifier, y = value, colour = .method_name, shape = .method_name
    ) +
    ggplot2::geom_point(alpha = 0.3, size = 0.5) +
    ggplot2::geom_hline(aes(yintercept = yintercept), colour = "black",
                        linetype = 2, alpha = 0.5) +
    ggplot2::facet_grid(cols = ggplot2::vars(n),
                        rows = ggplot2::vars(metric),
                        scales = "free_y") +
    ggplot2::scale_color_discrete(name = "Method") +
    ggplot2::scale_shape_discrete(name = "Method") +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        override.aes = list(size = 3, alpha = 1)
      )
    ) +
    ggplot2::xlab("Modifier Index") +
    ggplot2::ylab("Value") +
    ggplot2::ggtitle(("Treatment Effect Modifiers"))
    ggplot2::theme_bw()

  nontem_plt <- comb_plt %>%
    dplyr::filter(!(modifier %in% tem_indices)) %>%
    ggplot2::ggplot() +
    ggplot2::aes(
      x = modifier, y = value, colour = .method_name, shape = .method_name
    ) +
    ggplot2::geom_point(alpha = 0.3, size = 0.5) +
    ggplot2::geom_hline(aes(yintercept = yintercept), colour = "black",
                        linetype = 2, alpha = 0.5) +
    ggplot2::facet_grid(cols = ggplot2::vars(n),
                        rows = ggplot2::vars(metric),
                        scales = "free_y") +
    ggplot2::scale_color_discrete(name = "Method") +
    ggplot2::scale_shape_discrete(name = "Method") +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        override.aes = list(size = 3, alpha = 1)
      )
    ) +
    ggplot2::xlab("Modifier Index") +
    ggplot2::ylab("Value") +
    ggplot2::ggtitle(("Non-Modifiers"))
    ggplot2::theme_bw()

  return(plt)
}
