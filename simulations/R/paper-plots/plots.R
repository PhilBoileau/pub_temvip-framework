################################################################################
## Generate Figures for Paper
################################################################################

## load required libraries
library(here)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)

## establish common orderings
method_names <- c(
  "One Step", "TMLE", "Mod. Cov.", "Aug. Mod. Cov."
)
method_labels <- c(
  "One Step", "TML", "Mod. Cov.", "Aug. Mod. Cov."
)
bias_var_labels <- c("Absolute Bias", "Variance")
classification_labels <- c("FDR", "TNR", "TPR")

## bias/variance plotting function
bias_var_plot_fun <- function(eval_results, title, tem_indices) {

  ## combine the evaluator results into a single tibble
  emp_bias_tbl <- eval_results$`Empirical Bias` %>%
    mutate(
      metric = "Absolute Bias",
      value = abs(bias),
      yintercept = 0
    ) %>%
    select(-bias)
  emp_variance_tbl <- eval_results$`Empirical Variance` %>%
    mutate(
      metric = "Variance",
      value = variance,
      yintercept = 0
    ) %>%
    select(-variance)
  comb_tbl <- bind_rows(emp_bias_tbl, emp_variance_tbl) %>%
    filter(
      .method_name %in% c("One Step", "TMLE", "Mod. Cov.", "Aug. Mod. Cov.")
    )


  ## order the methods
  comb_tbl$.method_name <- factor(
    comb_tbl$.method_name,
    levels = method_names
  )
  levels(comb_tbl$.method_name) <- method_labels

  ## order metric
  comb_tbl$metric <- factor(
    comb_tbl$metric,
    levels = bias_var_labels
  )

  ## order and rename the sample sizes
  comb_tbl <- comb_tbl %>%
    mutate(
      n = paste0("n = ", n),
      n = factor(
        n, levels = c("n = 125", "n = 250", "n = 500", "n = 1000", "n = 2000"))
    )

  ## order the modifiers by index
  comb_tbl$modifier <- as.numeric(stringr::str_sub(comb_tbl$modifier, 2))

  ## create the plot
  tem_plt <- comb_tbl %>%
    filter(modifier %in% tem_indices) %>%
    ggplot() +
    aes(x = modifier, y = value, colour = .method_name, shape = .method_name) +
    geom_point(size = 1, alpha = 0.7) +
    geom_hline(
      aes(yintercept = yintercept), colour = "black", linetype = 2, alpha = 0.5
    ) +
    facet_grid(
      cols = ggplot2::vars(n), rows = ggplot2::vars(metric), scales = "free_y",

    ) +
    scale_color_brewer(palette = "Dark2", name = "Estimator") +
    scale_shape_discrete(name = "Estimator") +
    guides(
      colour = ggplot2::guide_legend(override.aes = list(size = 3, alpha = 1))
    ) +
    xlab("Confounder Index") +
    ylab(title) +
    ggtitle("Treatment Effect Modifiers") +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
    theme_bw()

  nontem_plt <- comb_tbl %>%
    filter(!(modifier %in% tem_indices)) %>%
    ggplot() +
    aes(x = modifier, y = value, colour = .method_name, shape = .method_name) +
    geom_point(alpha = 0.5, size = 0.5) +
    geom_hline(
      aes(yintercept = yintercept), colour = "black", linetype = 2, alpha = 0.5
    ) +
    facet_grid(
      cols = ggplot2::vars(n), rows = ggplot2::vars(metric), scales = "free_y",

    ) +
    scale_color_brewer(palette = "Dark2", name = "Estimator") +
    scale_shape_discrete(name = "Estimator") +
    guides(
      colour = ggplot2::guide_legend(override.aes = list(size = 3, alpha = 1))
    ) +
    xlab("Confounder Index") +
    ggtitle("Non-Modifiers") +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
    theme_bw() +
    theme(axis.title.y = element_blank())


    return(list(tem_plt, nontem_plt))
}

## TEM classification results plotting function
classification_plot_fun <- function(eval_results, title) {

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
    dplyr::select(-tnr) %>%
    filter(
      .method_name %in% c("One Step", "TMLE", "Mod. Cov.", "Aug. Mod. Cov.")
    )

  comb_tbl <- dplyr::bind_rows(emp_fdr_tbl, emp_tpr_tbl, emp_tnr_tbl)


  ## order the methods
  comb_tbl$.method_name <- factor(
    comb_tbl$.method_name,
    levels = method_names
  )
  levels(comb_tbl$.method_name) <- method_labels

  ## order metric
  comb_tbl$metric <- factor(
    comb_tbl$metric,
    levels = classification_labels
  )

  ## order and rename the sample sizes
  comb_tbl$n <- factor(
    comb_tbl$n,
    levels = c("125", "250", "500", "1000", "2000")
  )

  plt <- ggplot2::ggplot(comb_tbl) +
    aes(x = n, y = value, colour = .method_name, group = .method_name) +
    geom_point(alpha = 0.7) +
    geom_line(alpha = 0.7) +
    geom_hline(
      aes(yintercept = yintercept), colour = "black", linetype = 2, alpha = 0.5
    ) +
    facet_grid(rows = ggplot2::vars(metric)) +
    scale_color_brewer(palette = "Dark2", name = "Estimator") +
    xlab("Sample Size") +
    ylab("Empirical Value") +
    ggtitle(title) +
    theme_bw()

  return(plt)
}

## build the continuous outcome result plots ###################################

## load the simulation results
eval_results <- readRDS(
  here(paste0(
    "results/observational-continuous/Observational Study, Continuous Outcome/",
    "Varying n/eval_results.rds"
  ))
)

## bias and variance plot
continuous_bias_variance_plot <- bias_var_plot_fun(
  eval_results = eval_results,
  title = "Continuous Outcome",
  tem_indices = seq_len(5)
)

## classification results plot
continuous_classification_plot <- classification_plot_fun(
  eval_results = eval_results,
  title = "Continuous Outcome"
)

## build the binary outcome result plots #######################################

## load the simulation results
eval_results <- readRDS(
  here(paste0(
    "results/observational-binary/Observational Study, Binary Outcome/",
    "Varying n/eval_results.rds"
  ))
)

## bias and variance plot
binary_bias_variance_plot <- bias_var_plot_fun(
  eval_results = eval_results,
  title = "Binary Outcome",
  tem_indices = seq_len(5)
)

## classification results plot
binary_classification_plot <- classification_plot_fun(
  eval_results = eval_results,
  title = "Binary Outcome"
)

## build the time-to-event result plots ########################################

## load the simulation results
eval_results <- readRDS(
  here(paste0("results/RCT TTE/RCT, TTE Outcome/Varying n/eval_results.rds"))
)

## bias and variance plot
tte_bias_variance_plot <- bias_var_plot_fun(
  eval_results = eval_results,
  title = "Time-to-Event Outcome",
  tem_indices = seq_len(10)
)

## classification results plot
tte_classification_plot <- classification_plot_fun(
  eval_results = eval_results,
  title = "Time-to-Event Outcome"
)

## assembler bias-variance plots ###############################################

bias_variance_plot <- ggarrange(
  plotlist = list(
    continuous_bias_variance_plot[[1]], continuous_bias_variance_plot[[2]],
    binary_bias_variance_plot[[1]], binary_bias_variance_plot[[2]],
    tte_bias_variance_plot[[1]], tte_bias_variance_plot[[2]]
  ),
  ncol = 2,
  nrow = 3,
  labels = c("A", "", "B", "", "C", ""),
  label.x = 0,
  common.legend = TRUE,
  legend = "bottom"
)

ggexport(
  bias_variance_plot,
  filename = here("results/paper-plots/bias-variance-plot.jpeg"),
  width = 1600,
  height = 1200,
  res = 125
)

## assembler classification results plots ######################################

classification_results_plot <- ggarrange(
  plotlist = list(continuous_classification_plot, binary_classification_plot,
                  tte_classification_plot),
  ncol = 1,
  labels = "AUTO",
  label.x = 0,
  common.legend = TRUE,
  legend = "bottom",
  align = "hv"
)

ggexport(
  classification_results_plot,
  filename = here("results/paper-plots/classification-plot.jpeg"),
  width = 1200,
  height = 1200,
  res = 120
)
