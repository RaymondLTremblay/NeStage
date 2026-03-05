# =============================================================================
# NeStage ggplot2 plots
# =============================================================================

nestage_plot_effective_sizes_gg <- function(df, title = "Effective sizes by model") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for this plot.")
  }
  stopifnot(all(c("model","NyN","NeN") %in% names(df)))
  library(ggplot2)

  df_long <- rbind(
    data.frame(model = df$model, metric = "Ny/N", value = df$NyN),
    data.frame(model = df$model, metric = "Ne/N", value = df$NeN)
  )

  ggplot(df_long, aes(x = value, y = reorder(model, value), color = metric)) +
    geom_point(size = 2) +
    facet_wrap(~ metric, scales = "free_x") +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")
}

nestage_plot_compare_gg <- function(cmp, title = "Model differences (b vs a)") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for this plot.")
  }
  library(ggplot2)
  stopifnot(all(c("compare","d_NyN","d_NeN") %in% names(cmp)))

  cmp_long <- rbind(
    data.frame(compare = cmp$compare, metric = "ΔNy/N", value = cmp$d_NyN),
    data.frame(compare = cmp$compare, metric = "ΔNe/N", value = cmp$d_NeN)
  )

  ggplot(cmp_long, aes(x = value, y = reorder(compare, value), color = metric)) +
    geom_point(size = 2) +
    facet_wrap(~ metric, scales = "free_x") +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")
}
