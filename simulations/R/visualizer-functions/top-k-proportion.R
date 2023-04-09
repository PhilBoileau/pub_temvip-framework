################################################################################
# Top k modifier selection plotting function
################################################################################

top_selection_prop_plot_fun <- function(eval_results) {

  results <- eval_results$`Proportion Correct Top 10`
  # order the methods
  results$.method_name <- factor(
    results$.method_name,
    levels = c(
      "Mod. Cov.", "Aug. Mod. Cov.",
      "One Step", "CF One Step",
      "TMLE", "CF TMLE"
    )
  )

  plt <- ggplot2::ggplot(results) +
    ggplot2::aes(x = n, y = prop_top_k, colour = .method_name) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::geom_line(alpha = 0.7) +
    ggplot2::scale_color_discrete(name = "Method") +
    ggplot2::xlab("Sample Size") +
    ggplot2::ylab("Average Proportion Top 10 Modiers Selected") +
    ggplot2::theme_bw()

  return(plt)
}
