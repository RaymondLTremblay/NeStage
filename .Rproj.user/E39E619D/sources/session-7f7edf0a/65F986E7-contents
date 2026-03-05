## scripts/02_plots_gg.R

nestage_rename_populations <- function(df, map, order_levels = TRUE) {
  stopifnot("population" %in% names(df))
  idx <- match(df$population, names(map))
  df$population <- ifelse(is.na(idx), df$population, unname(map[idx]))
  if (order_levels) df$population <- factor(df$population, levels = unique(unname(map)))
  df
}

nestage_plot_ratios_gg <- function(df) {
  req <- c("population","model","NyN","NeN")
  if (!all(req %in% names(df))) stop("df must contain: ", paste(req, collapse=", "))

  dd <- rbind(
    data.frame(population=df$population, model=df$model, metric="Ny/N", value=df$NyN),
    data.frame(population=df$population, model=df$model, metric="Ne/N", value=df$NeN)
  )
  dd <- dd[is.finite(dd$value), , drop = FALSE]
  if (nrow(dd) == 0L) stop("No finite Ny/N or Ne/N to plot (estimators still NA?).")

  dd$population <- factor(dd$population, levels = unique(dd$population))

  ggplot(dd, aes(x = population, y = value, color = model, shape = model)) +
    geom_point(size = 2.2, position = position_dodge(width = 0.6), na.rm = TRUE) +
    facet_wrap(~ metric, scales = "free_y", nrow = 1) +
    labs(title = "Ny/N and Ne/N across models", x = "Population", y = NULL,
         color = "Model", shape = "Model") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold"))
}

nestage_plot_generation_time_gg <- function(df) {
  req <- c("population","model","L")
  if (!all(req %in% names(df))) stop("df must contain: ", paste(req, collapse=", "))

  dL <- df[grepl("^Sex-only", df$model) & is.finite(df$L), , drop = FALSE]
  if (nrow(dL) == 0L) stop("No finite L for Sex-only models (is L computed?).")

  dL$population <- factor(dL$population, levels = unique(dL$population))

  ggplot(dL, aes(x = population, y = L)) +
    geom_point(size = 2.4, na.rm = TRUE) +
    labs(title = "Generation time (L) — Sex-only", x = "Population", y = "L (years)") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold"))
}
