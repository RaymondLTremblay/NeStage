# R/nestage-plots.R

#' @importFrom graphics par axis plot
NULL

#' Plot Ny/N and Ne/N by population and model (points only)
#' @param df Tidy table from nestage_run_all() or combined tables.
#' @export
nestage_plot_ratios <- function(df) {
  op <- par(mfrow = c(1, 2))
  on.exit(par(op), add = TRUE)
  x <- seq_len(nrow(df))
  # Ny/N
  plot(x, df$NyN, pch = 16, xaxt = "n", xlab = "Population x Model", ylab = "Ny/N",
       main = "Ny/N across models")
  axis(1, at = x, labels = paste(df$population, df$model, sep = " | "), las = 2, cex.axis = 0.6)
  # Ne/N
  plot(x, df$NeN, pch = 16, xaxt = "n", xlab = "Population x Model", ylab = "Ne/N",
       main = "Ne/N across models")
  axis(1, at = x, labels = paste(df$population, df$model, sep = " | "), las = 2, cex.axis = 0.6)
}

#' Plot generation time (L) for Sex-only models (points only)
#' @param df Tidy table from nestage_run_all() or combined tables.
#' @export
nestage_plot_generation_time <- function(df) {
  dL <- df[grepl("^Sex-only", df$model), , drop = FALSE]
  x <- seq_len(nrow(dL))
  plot(x, dL$L, pch = 16, xaxt = "n", xlab = "Population", ylab = "L (years)",
       main = "Generation time (L) - Sex-only")
  axis(1, at = x, labels = dL$population, las = 2, cex.axis = 0.7)
}
