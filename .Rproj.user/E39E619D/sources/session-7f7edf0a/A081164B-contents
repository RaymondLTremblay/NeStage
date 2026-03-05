# =============================================================================
# NeStage base R plots
# =============================================================================

# Dot plot: Ny/N and Ne/N by model (points only, no bars)
nestage_plot_effective_sizes_base <- function(df, main = "Effective sizes by model") {
  stopifnot(all(c("model","NyN","NeN") %in% names(df)))
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(mfrow = c(1,2), mar = c(6,4,3,1))

  # Ny/N
  ord <- order(df$NyN, decreasing = TRUE)
  dotchart(df$NyN[ord], labels = df$model[ord], main = paste(main, "Ny/N"), xlab = "Ny/N")

  # Ne/N
  ord2 <- order(df$NeN, decreasing = TRUE)
  dotchart(df$NeN[ord2], labels = df$model[ord2], main = paste(main, "Ne/N"), xlab = "Ne/N")
  invisible(TRUE)
}

# Difference plot (delta) between models (if you computed nestage_compare)
nestage_plot_compare_base <- function(cmp, main = "Model differences (b vs a)") {
  stopifnot(all(c("compare","d_NyN","d_NeN") %in% names(cmp)))
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(mfrow = c(1,2), mar = c(6,4,3,1))
  dotchart(cmp$d_NyN, labels = cmp$compare, main = paste(main, "ΔNy/N"), xlab = "ΔNy/N")
  dotchart(cmp$d_NeN, labels = cmp$compare, main = paste(main, "ΔNe/N"), xlab = "ΔNe/N")
  invisible(TRUE)
}
