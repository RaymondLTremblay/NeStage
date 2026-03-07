# =============================================================================
# Ne_compadre_sexual_analysis.R
# -----------------------------------------------------------------------------
# Analysis script: compute Ne/N for sexually reproducing plant species
# from the COMPADRE Plant Matrix Database, summarise by family, and plot.
#
# Workflow:
#   Step 1  Install packages (once)
#   Step 2  Load packages + source function file
#   Step 3  Filter COMPADRE             [~1 min, requires internet]
#   Step 4  Save / reload filtered db  [avoid re-fetching every session]
#   Step 5  Run Ne analysis             [seconds]
#   Step 6  Inspect results
#   Step 7  Plots
#   Step 8  Save results to CSV
#
# Author:  Raymond L. Tremblay  <raymond.tremblay@upr.edu>
# Date:    2026-03-06
# =============================================================================


# =============================================================================
# STEP 1 -- Install packages (run once, then comment out)
# =============================================================================
# install.packages("Rcompadre")
# install.packages("ggplot2")
# install.packages("dplyr")
# install.packages("remotes")
# remotes::install_github("RaymondLTremblay/NeStage")


# =============================================================================
# STEP 2 -- Load packages and source the function file
# =============================================================================
library(Rcompadre)
library(NeStage)
library(ggplot2)
library(dplyr)

# Source the COMPADRE-NeStage bridge functions.
# Adjust the path to wherever you saved Ne_compadre_sexual.R
source("Ne_compadre_sexual.R")


# =============================================================================
# STEP 3 -- Filter COMPADRE
# (requires internet; prints a step-by-step count and family table)
# =============================================================================
db_sexual <- Ne_compadre_sexual_filter(verbose = TRUE)

# Quick check
nrow(db_sexual)              # how many species passed all filters
table(db_sexual$Family)      # species per family
table(db_sexual$Order)       # species per order


# =============================================================================
# STEP 4 -- Save / reload the filtered database
# Save once so you can reload without re-fetching COMPADRE every session.
# =============================================================================

## -- save (run after Step 3) -------------------------------------------------
saveRDS(db_sexual, "db_sexual_compadre.rds")

## -- reload (use instead of Steps 3-4 in future sessions) -------------------
# db_sexual <- readRDS("db_sexual_compadre.rds")


# =============================================================================
# STEP 5 -- Run Ne analysis
# D is derived from the stable stage distribution of A = matU + matF.
# Poisson reproductive variance assumed (Vk/k = 1) -- most conservative.
# =============================================================================
out <- Ne_compadre_sexual(
  db        = db_sexual,
  Ne_target = 5000,       # Ne conservation threshold
  Vk_over_k = NULL,       # NULL = Poisson; set e.g. 2 for over-dispersed
  verbose   = TRUE        # one line per species; set FALSE to suppress
)

# Pull out the summary table
summary_df <- out$summary


# =============================================================================
# STEP 6 -- Inspect results
# =============================================================================

# --- 6a: Full table -----------------------------------------------------------
print(summary_df, digits = 4)

# --- 6b: Overview statistics --------------------------------------------------
cat("\n==============================\n")
cat("  COMPADRE Ne/N OVERVIEW\n")
cat("==============================\n")
cat(sprintf("  Species analysed : %d\n",   nrow(summary_df)))
cat(sprintf("  Families covered : %d\n",   length(unique(summary_df$Family))))
cat(sprintf("  Ne/N range       : %.4f -- %.4f\n",
            min(summary_df$NeN), max(summary_df$NeN)))
cat(sprintf("  Median Ne/N      : %.4f\n", median(summary_df$NeN)))
cat(sprintf("  Min_N range      : %d -- %d\n",
            min(summary_df$Min_N), max(summary_df$Min_N)))
cat("==============================\n\n")

# --- 6c: Species below Frankham (1995) empirical median (Ne/N = 0.10) --------
cat("Species with Ne/N < 0.10  [below Frankham 1995 empirical median]:\n")
low_Ne <- summary_df |>
  filter(NeN < 0.10) |>
  arrange(NeN) |>
  select(Species, Family, NeN, Min_N, L, lambda)
print(low_Ne)

# --- 6d: Species needing the largest census populations ----------------------
cat("\nTop 10 species by minimum census size needed for Ne >= 5000:\n")
summary_df |>
  arrange(desc(Min_N)) |>
  slice_head(n = 10) |>
  select(Species, Family, NeN, Min_N, L, MatrixDimension) |>
  print()

# --- 6e: Access a full Ne_sexual_Y2000 result for one species ----------------
# Replace the species name below with any name in out$results
first_sp <- names(out$results)[1]
cat("\nFull result for:", first_sp, "\n")
print(out$results[[first_sp]])

# --- 6f: Family-level summary (mean Ne/N, n species) -------------------------
fam_summary <- summary_df |>
  group_by(Family) |>
  summarise(
    n_species   = n(),
    mean_NeN    = round(mean(NeN),    3),
    median_NeN  = round(median(NeN),  3),
    min_NeN     = round(min(NeN),     3),
    max_NeN     = round(max(NeN),     3),
    mean_L      = round(mean(L),      2),
    .groups     = "drop"
  ) |>
  arrange(median_NeN)

cat("\nFamily-level Ne/N summary (sorted by median Ne/N):\n")
print(fam_summary, n = Inf)


# =============================================================================
# STEP 7 -- Plots
# =============================================================================

# --- Plot 1: Distribution of Ne/N across all species -------------------------
p1 <- ggplot(summary_df, aes(x = NeN)) +
  geom_histogram(bins = 30, fill = "#2166ac", colour = "white", alpha = 0.85) +
  geom_vline(xintercept = 0.10, linetype = "dashed",
             colour = "#d6604d", linewidth = 0.9) +
  annotate("text", x = 0.105, y = Inf,
           label = "Frankham (1995)\nmedian = 0.10",
           hjust = 0, vjust = 1.6, size = 3.2, colour = "#d6604d") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    title    = expression("Distribution of " * N[e]/N * " across plant species (COMPADRE)"),
    subtitle = "Sexual species -- mean matrices -- wild, unmanipulated populations",
    x        = expression(N[e] / N),
    y        = "Number of species",
    caption  = "Yonezawa et al. (2000) Eq. 6 (d = 0). D from stable stage distribution."
  ) +
  theme_classic(base_size = 12)

print(p1)
ggsave("Ne_compadre_distribution.png", p1, width = 7, height = 4.5, dpi = 180)


# --- Plot 2: Ne/N by family (families with >= 2 species) ---------------------
fam_min2 <- names(which(table(summary_df$Family) >= 2))
df_fam   <- summary_df |> filter(Family %in% fam_min2)

p2 <- ggplot(df_fam,
             aes(x = reorder(Family, NeN, FUN = median), y = NeN)) +
  geom_boxplot(fill = "#4dac26", alpha = 0.55,
               outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.15, size = 1.8, alpha = 0.75, colour = "#1b7837") +
  geom_hline(yintercept = 0.10, linetype = "dashed",
             colour = "#d6604d", linewidth = 0.8) +
  annotate("text", x = 0.7, y = 0.103,
           label = "Frankham 0.10", hjust = 0, size = 2.9, colour = "#d6604d") +
  coord_flip() +
  labs(
    title    = expression("Ne/N by plant family (COMPADRE)"),
    subtitle = "Families with \u2265 2 species -- sexual, mean matrices",
    x        = NULL,
    y        = expression(N[e] / N),
    caption  = "Red dashed line = Frankham (1995) empirical median."
  ) +
  theme_classic(base_size = 12)

print(p2)
ggsave("Ne_compadre_by_family.png", p2, width = 7, height = 6, dpi = 180)


# --- Plot 3: Ne/N vs generation time L, coloured by family -------------------
p3 <- ggplot(summary_df,
             aes(x = L, y = NeN, colour = Family, label = Species)) +
  geom_smooth(method  = "lm", se = TRUE,
              colour  = "grey50", fill = "grey85",
              linetype = "dashed", linewidth = 0.7,
              inherit.aes = FALSE,
              mapping = aes(x = L, y = NeN)) +
  geom_point(size = 2.5, alpha = 0.85) +
  geom_hline(yintercept = 0.10, linetype = "dotted",
             colour = "#d6604d", linewidth = 0.6) +
  labs(
    title    = expression(N[e]/N * " vs generation time " * italic(L)),
    subtitle = "Sexually reproducing plant species from COMPADRE",
    x        = "Generation time L (years)",
    y        = expression(N[e] / N),
    colour   = "Family",
    caption  = "Red dotted line = Frankham (1995) median."
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position    = "right",
        legend.text        = element_text(size = 7),
        legend.key.size    = unit(0.35, "cm"),
        legend.title       = element_text(size = 8))

print(p3)
ggsave("Ne_compadre_NeN_vs_L.png", p3, width = 8.5, height = 5.5, dpi = 180)


# --- Plot 4: Min_N (census size for Ne >= 5000) by family --------------------
p4 <- ggplot(summary_df,
             aes(x = reorder(Family, Min_N, FUN = median), y = Min_N)) +
  geom_point(colour = "#762a83", size = 2.2, alpha = 0.8) +
  geom_hline(yintercept = 5000, linetype = "dashed",
             colour = "grey40", linewidth = 0.7) +
  scale_y_log10(labels = scales::comma) +
  coord_flip() +
  labs(
    title    = "Minimum census size (N) needed for Ne \u2265 5000",
    subtitle = "One point per species, sorted by family median",
    x        = NULL,
    y        = "Minimum N (log scale)",
    caption  = "Grey dashed line = Ne target of 5000."
  ) +
  theme_classic(base_size = 12)

print(p4)
ggsave("Ne_compadre_minN_by_family.png", p4, width = 7, height = 6, dpi = 180)


# =============================================================================
# STEP 8 -- Save results to CSV
# =============================================================================
write.csv(summary_df, "Ne_compadre_sexual_results.csv", row.names = FALSE)
write.csv(fam_summary, "Ne_compadre_family_summary.csv", row.names = FALSE)

cat("\nFiles saved:\n")
cat("  Ne_compadre_sexual_results.csv  -- one row per species\n")
cat("  Ne_compadre_family_summary.csv  -- family-level statistics\n")
cat("  Ne_compadre_distribution.png\n")
cat("  Ne_compadre_by_family.png\n")
cat("  Ne_compadre_NeN_vs_L.png\n")
cat("  Ne_compadre_minN_by_family.png\n")
cat("  db_sexual_compadre.rds          -- filtered COMPADRE db\n\n")


# =============================================================================
# Session info (for reproducibility)
# =============================================================================
sessionInfo()
