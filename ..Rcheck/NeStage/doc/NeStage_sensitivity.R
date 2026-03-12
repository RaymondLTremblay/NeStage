## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse  = TRUE,
  comment   = "#>",
  fig.width = 7,
  fig.height = 4.5,
  warning   = FALSE,
  message   = FALSE
)

## ----load-nestage, include=FALSE----------------------------------------------
library(NeStage)

## ----packages, include=FALSE--------------------------------------------------
library(ggplot2)
library(dplyr)
library(purrr)
library(tidyr)
library(knitr)

## ----setup-populations--------------------------------------------------------
# ---- L. rupestris population 1 (Ne/N = 0.572) ----
MatU_Lrup1 <- matrix(c(
  0,     0,     0,     0,     0,     0,
  0.738, 0.738, 0,     0,     0,     0,
  0,     0,     0.515, 0,     0.076, 0.013,
  0,     0.038, 0,     0.777, 0,     0,
  0,     0.002, 0.368, 0.011, 0.730, 0.171,
  0,     0,     0.037, 0,     0.169, 0.790
), nrow = 6, byrow = TRUE,
dimnames = list(paste0("s", 1:6), paste0("s", 1:6)))

F_Lrup1 <- c(0, 0, 0, 0, 0.120, 0.414)   # monthly fecundity
D_Lrup1 <- c(22, 22, 36, 36, 48.5, 48.5)
D_Lrup1 <- D_Lrup1 / sum(D_Lrup1)

# ---- L. rubripetala population 1 (Ne/N = 0.930) ----
MatU_Lrub1 <- matrix(c(
  0,     0,     0,     0,     0,
  0.667, 0.667, 0,     0,     0,
  0,     0,     0.825, 0,     0.227,
  0,     0.037, 0,     0.812, 0,
  0,     0,     0.137, 0.010, 0.739
), nrow = 5, byrow = TRUE,
dimnames = list(paste0("s", 1:5), paste0("s", 1:5)))

F_Lrub1 <- c(0, 0, 0, 0, 0.043)
N_Lrub1 <- c(9/2, 9/2, 44/2, 44/2, 101)
D_Lrub1 <- N_Lrub1 / sum(N_Lrub1)

# ---- L. rubripetala population 5 (Ne/N = 0.250, lowest in Lrub) ----
MatU_Lrub5 <- matrix(c(
  0,     0,     0,     0,     0,
  0.667, 0.667, 0,     0,     0,
  0,     0,     0.726, 0,     0.250,
  0,     0.037, 0,     0.812, 0,
  0,     0,     0.237, 0.010, 0.739
), nrow = 5, byrow = TRUE,
dimnames = list(paste0("s", 1:5), paste0("s", 1:5)))

F_Lrub5 <- c(0, 0, 0, 0, 0.067)
N_Lrub5 <- c(52/2, 52/2, 4/2, 4/2, 46)
D_Lrub5 <- N_Lrub5 / sum(N_Lrub5)

# Baseline Ne estimates for reference
Ne_base_Lrup1 <- Ne_sexual_Y2000(
  T_mat = MatU_Lrup1, F_vec = F_Lrup1, D = D_Lrup1,
  Ne_target = 50, census_N = 40, population = "L. rupestris pop 1"
)
Ne_base_Lrub1 <- Ne_sexual_Y2000(
  T_mat = MatU_Lrub1, F_vec = F_Lrub1, D = D_Lrub1,
  Ne_target = 50, census_N = 40, population = "L. rubripetala pop 1"
)
Ne_base_Lrub5 <- Ne_sexual_Y2000(
  T_mat = MatU_Lrub5, F_vec = F_Lrub5, D = D_Lrub5,
  Ne_target = 50, census_N = 40, population = "L. rubripetala pop 5"
)

# Baseline summary
data.frame(
  Population = c("L. rupestris pop 1", "L. rubripetala pop 1",
                 "L. rubripetala pop 5"),
  Stages     = c(6, 5, 5),
  NeN        = round(c(Ne_base_Lrup1$NeN, Ne_base_Lrub1$NeN,
                        Ne_base_Lrub5$NeN), 3),
  Ne_at_40   = c(Ne_base_Lrup1$Ne_at_census, Ne_base_Lrub1$Ne_at_census,
                  Ne_base_Lrub5$Ne_at_census),
  Min_N_50   = c(Ne_base_Lrup1$Min_N, Ne_base_Lrub1$Min_N,
                  Ne_base_Lrub5$Min_N)
) |>
  knitr::kable(
    col.names = c("Population", "Stages", "Ne/N",
                  "Ne at N=40", "Min N (Ne≥50)"),
    caption   = "Baseline Ne/N estimates for the three study populations
                 (monthly matrices, Ne_target = 50, census_N = 40)."
  )

## ----vk-sensitivity-lrup1, fig.cap="Sensitivity of Ne/N to sexual reproductive variance in stage 5 of *L. rupestris* pop 1. Dashed vertical line = Poisson default."----
sens_Vk_Lrup1 <- Ne_sensitivity_Vk(
  model_fn    = Ne_sexual_Y2000,
  T_mat       = MatU_Lrup1,
  F_vec       = F_Lrup1,
  D           = D_Lrup1,
  stage_index = 5,              # stage 5: single-shoot reproductive adults
  Vk_range    = seq(0.5, 8, by = 0.5),
  Ne_target   = 50,
  population  = "L. rupestris pop 1"
)

print(sens_Vk_Lrup1)

## ----vk-both-stages, fig.cap="Sensitivity of Ne/N to sexual reproductive variance: stage 5 (left) vs stage 6 (right) in *L. rupestris* pop 1.", fig.width=10----
# Also sweep stage 6 (two-or-more-shoot adults) for comparison
sens_Vk_Lrup1_s6 <- Ne_sensitivity_Vk(
  model_fn    = Ne_sexual_Y2000,
  T_mat       = MatU_Lrup1,
  F_vec       = F_Lrup1,
  D           = D_Lrup1,
  stage_index = 6,
  Vk_range    = seq(0.5, 8, by = 0.5),
  Ne_target   = 50,
  population  = "L. rupestris pop 1 — stage 6"
)

# Combined data for comparison plot
df_vk_compare <- bind_rows(
  sens_Vk_Lrup1$data    %>% mutate(Stage = "Stage 5 (1-shoot repro)"),
  sens_Vk_Lrup1_s6$data %>% mutate(Stage = "Stage 6 (2+-shoot repro)")
)

ggplot(df_vk_compare, aes(x = Vk_over_k, y = NeN, color = Stage)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("#2166ac", "#d6604d")) +
  labs(
    title    = expression(paste("Ne/N vs. sexual reproductive variance  ", V[k]/bar(k))),
    subtitle = "L. rupestris pop 1 — comparing focal stage",
    x        = expression(V[k]/bar(k)),
    y        = expression(N[e]/N),
    color    = NULL,
    caption  = "Dashed line = Poisson default (Vk/k̄ = 1). Monthly matrices."
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.4)
  )

## ----L-sensitivity, fig.cap="Sensitivity of Ne/N to generation time L across all three study populations.", fig.width=10, fig.height=5----
L_range_use <- seq(5, 80, by = 5)

# Run for all three populations
make_L_sens <- function(MatU, F_vec, D, Ne_target, pop_label) {
  Ne_sensitivity_L(
    model_fn   = Ne_sexual_Y2000,
    T_mat      = MatU,
    F_vec      = F_vec,
    D          = D,
    L_range    = L_range_use,
    Ne_target  = Ne_target,
    population = pop_label
  )
}

sens_L_Lrup1 <- make_L_sens(MatU_Lrup1, F_Lrup1, D_Lrup1, 50,
                              "L. rupestris pop 1")
sens_L_Lrub1 <- make_L_sens(MatU_Lrub1, F_Lrub1, D_Lrub1, 50,
                              "L. rubripetala pop 1")
sens_L_Lrub5 <- make_L_sens(MatU_Lrub5, F_Lrub5, D_Lrub5, 50,
                              "L. rubripetala pop 5")

# Observed L values from monthly matrices
L_obs <- c(Ne_base_Lrup1$L, Ne_base_Lrub1$L, Ne_base_Lrub5$L)

df_L <- bind_rows(
  sens_L_Lrup1$data %>% mutate(Population = "L. rupestris pop 1",
                                 L_obs = L_obs[1]),
  sens_L_Lrub1$data %>% mutate(Population = "L. rubripetala pop 1",
                                 L_obs = L_obs[2]),
  sens_L_Lrub5$data %>% mutate(Population = "L. rubripetala pop 5",
                                 L_obs = L_obs[3])
)

ggplot(df_L, aes(x = L, y = NeN, color = Population)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.8) +
  geom_vline(aes(xintercept = L_obs, color = Population),
             linetype = "dashed", linewidth = 0.7, alpha = 0.7) +
  scale_color_manual(values = c("#2166ac", "#1b7837", "#d6604d")) +
  labs(
    title    = expression(paste("Sensitivity of  ", N[e]/N, "  to generation time  ", L)),
    subtitle = "Dashed vertical lines = observed L from monthly matrices",
    x        = "Generation time L (months)",
    y        = expression(N[e]/N),
    color    = NULL,
    caption  = "All populations: Ne_sexual_Y2000, Poisson defaults, Ne_target = 50."
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title      = element_text(face = "bold"),
    legend.position = "bottom",
    legend.text     = element_text(face = "italic"),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.4)
  )

## ----L-print------------------------------------------------------------------
# Print elasticity summary for the most sensitive population
print(sens_L_Lrub5)

## ----d-sensitivity, fig.cap="Effect of clonal fraction on Ne/N: flat under Poisson variance (left), declining under super-Poisson clonal variance (right).", fig.width=10----
# Use L. rubripetala pop 1 matrix structure with a mixed reproduction scenario
# Add a clonal component: imagining a hypothetical scenario where stage 5
# (reproductive adults) can reproduce clonally via offshoots

# Poisson case — d has no effect
sens_d_poisson <- Ne_sensitivity_d(
  T_mat       = MatU_Lrub1,
  F_vec       = F_Lrub1,
  D           = D_Lrub1,
  d_fixed     = c(0, 0, 0, 0, 0.5),   # base: stage 5 is 50% clonal
  stage_index = 5,
  d_range     = seq(0, 1, by = 0.1),
  Vc_over_c   = rep(1, 5),             # Poisson clonal variance
  Ne_target   = 50,
  population  = "L. rubripetala pop 1 (Poisson)"
)

# Super-Poisson clonal variance — d starts to matter
sens_d_superpoisson <- Ne_sensitivity_d(
  T_mat       = MatU_Lrub1,
  F_vec       = F_Lrub1,
  D           = D_Lrub1,
  d_fixed     = c(0, 0, 0, 0, 0.5),
  stage_index = 5,
  d_range     = seq(0, 1, by = 0.1),
  Vc_over_c   = c(1, 1, 1, 1, 4),     # stage 5 clonal output highly variable
  Ne_target   = 50,
  population  = "L. rubripetala pop 1 (Vc/c̄ = 4)"
)

df_d <- bind_rows(
  sens_d_poisson$data      %>% mutate(Scenario = "Poisson (Vc/c̄ = 1)"),
  sens_d_superpoisson$data %>% mutate(Scenario = "Super-Poisson (Vc/c̄ = 4)")
)

ggplot(df_d, aes(x = d_focal, y = NeN, color = Scenario)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  geom_vline(xintercept = 0.5, linetype = "dashed",
             color = "grey50", linewidth = 0.7) +
  annotate("text", x = 0.51, y = max(df_d$NeN) * 0.98,
           label = "Base d = 0.5", hjust = 0, size = 3.2, color = "grey40") +
  scale_color_manual(values = c("#1b7837", "#d6604d")) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
  labs(
    title    = expression(paste("Ne/N vs. clonal fraction  ", d[i], "  in stage 5")),
    subtitle = "L. rubripetala pop 1 — hypothetical mixed reproduction scenario",
    x        = expression(paste("Clonal fraction  ", d[i],
                                "  (0 = fully sexual,  1 = fully clonal)")),
    y        = expression(N[e]/N),
    color    = NULL,
    caption  = "Ne_mixed_Y2000. Flat curve confirms Poisson result of Yonezawa et al. (2000)."
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title      = element_text(face = "bold"),
    legend.position = "bottom",
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.4)
  )

## ----d-print-superpoisson-----------------------------------------------------
print(sens_d_superpoisson)

## ----elasticity-comparison----------------------------------------------------
# Extract local elasticity from each sensitivity object
# Local elasticity approximated by central finite difference at Vk/k_bar = 1
# E = (dNeN/dtheta) * (theta_ref / NeN_ref)

local_elasticity <- function(sens_obj, param_col, ref_val) {
  df  <- sens_obj$data
  idx <- which.min(abs(df[[param_col]] - ref_val))
  # Use neighbouring points for central difference
  if (idx == 1)        idx <- 2
  if (idx == nrow(df)) idx <- nrow(df) - 1
  dNeN   <- df$NeN[idx + 1] - df$NeN[idx - 1]
  dtheta <- df[[param_col]][idx + 1] - df[[param_col]][idx - 1]
  NeN_ref <- df$NeN[idx]
  (dNeN / dtheta) * (ref_val / NeN_ref)
}

# L. rupestris pop 1: sweep Vk stages 5 and 6, and L
E_Vk5_Lrup1 <- local_elasticity(sens_Vk_Lrup1,    "Vk_over_k", 1)
E_Vk6_Lrup1 <- local_elasticity(sens_Vk_Lrup1_s6, "Vk_over_k", 1)
E_L_Lrup1   <- local_elasticity(sens_L_Lrup1,      "L",
                                  Ne_base_Lrup1$L)

# L. rubripetala pop 1
sens_Vk_Lrub1 <- Ne_sensitivity_Vk(
  model_fn    = Ne_sexual_Y2000,
  T_mat       = MatU_Lrub1,
  F_vec       = F_Lrub1,
  D           = D_Lrub1,
  stage_index = 5,
  Vk_range    = seq(0.5, 8, by = 0.5),
  Ne_target   = 50,
  population  = "L. rubripetala pop 1"
)
E_Vk5_Lrub1 <- local_elasticity(sens_Vk_Lrub1, "Vk_over_k", 1)
E_L_Lrub1   <- local_elasticity(sens_L_Lrub1,   "L", Ne_base_Lrub1$L)

# L. rubripetala pop 5
sens_Vk_Lrub5 <- Ne_sensitivity_Vk(
  model_fn    = Ne_sexual_Y2000,
  T_mat       = MatU_Lrub5,
  F_vec       = F_Lrub5,
  D           = D_Lrub5,
  stage_index = 5,
  Vk_range    = seq(0.5, 8, by = 0.5),
  Ne_target   = 50,
  population  = "L. rubripetala pop 5"
)
E_Vk5_Lrub5 <- local_elasticity(sens_Vk_Lrub5, "Vk_over_k", 1)
E_L_Lrub5   <- local_elasticity(sens_L_Lrub5,   "L", Ne_base_Lrub5$L)

# Assemble table
elast_tbl <- data.frame(
  Population = c("*L. rupestris* pop 1",
                 "*L. rubripetala* pop 1",
                 "*L. rubripetala* pop 5"),
  Base_NeN   = round(c(Ne_base_Lrup1$NeN,
                         Ne_base_Lrub1$NeN,
                         Ne_base_Lrub5$NeN), 3),
  E_Vk_repro = round(c(E_Vk5_Lrup1, E_Vk5_Lrub1, E_Vk5_Lrub5), 3),
  E_L        = round(c(E_L_Lrup1,   E_L_Lrub1,   E_L_Lrub5),   3)
)

knitr::kable(elast_tbl,
  col.names = c("Population", "Base Ne/N",
                "E(Vk/k̄) repro stage", "E(L)"),
  escape  = FALSE,
  caption = "Local elasticities at the Poisson reference point (Vk/k̄ = 1)
             and at the observed L. E(Vk/k̄) < 0: increasing reproductive
             variance reduces Ne/N. E(L) < 0: longer generation time
             reduces Ne/N. Larger absolute values indicate higher sensitivity."
)

## ----elasticity-plot, fig.cap="Absolute local elasticity of Ne/N with respect to two parameters across three populations. Larger bars indicate parameters where management intervention would have the greatest proportional effect on Ne/N.", fig.height=4----
elast_long <- elast_tbl %>%
  select(Population, E_Vk_repro, E_L) %>%
  pivot_longer(cols = c(E_Vk_repro, E_L),
               names_to = "Parameter",
               values_to = "Elasticity") %>%
  mutate(
    Abs_E     = abs(Elasticity),
    Parameter = recode(Parameter,
                       E_Vk_repro = "Sexual repro variance\n(Vk/k̄, repro stage)",
                       E_L        = "Generation time (L)")
  )

ggplot(elast_long, aes(x = Population, y = Abs_E, fill = Parameter)) +
  geom_col(position = "dodge", width = 0.6) +
  scale_fill_manual(values = c("#2166ac", "#d6604d")) +
  scale_x_discrete(labels = function(x) gsub("\\*", "", x)) +
  labs(
    title    = "|Elasticity| of Ne/N: which parameter drives genetic drift most?",
    subtitle = "Higher = Ne/N more sensitive to this parameter",
    x        = NULL,
    y        = "|Local elasticity|",
    fill     = "Parameter",
    caption  = "Local elasticity at Poisson reference (Vk/k̄ = 1) and observed L.
                All populations: Ne_sexual_Y2000, monthly matrices."
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title      = element_text(face = "bold"),
    axis.text.x     = element_text(face = "italic", size = 9),
    legend.position = "right",
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.4)
  )

## ----final-comparison, fig.cap="Ne/N as a function of Vk/k̄ in the reproductive stage for all three study populations. Populations with higher baseline Ne/N are also less sensitive to reproductive variance.", fig.height=5----
df_final <- bind_rows(
  sens_Vk_Lrup1$data %>% mutate(Population = "*L. rupestris* pop 1"),
  sens_Vk_Lrub1$data %>% mutate(Population = "*L. rubripetala* pop 1"),
  sens_Vk_Lrub5$data %>% mutate(Population = "*L. rubripetala* pop 5")
)

ggplot(df_final, aes(x = Vk_over_k, y = NeN, color = Population)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_vline(xintercept = 1, linetype = "dashed",
             color = "grey50", linewidth = 0.7) +
  geom_hline(yintercept = 50/40, linetype = "dotted",
             color = "grey30", linewidth = 0.6) +
  annotate("text", x = 7.5, y = 50/40 + 0.03,
           label = "Ne=50 at N=40", hjust = 1, size = 3, color = "grey30") +
  scale_color_manual(
    values = c("#2166ac", "#1b7837", "#d6604d"),
    labels = c(
      expression(italic("L. rubripetala")~"pop 1"),
      expression(italic("L. rubripetala")~"pop 5"),
      expression(italic("L. rupestris")~"pop 1")
    )
  ) +
  labs(
    title    = expression(paste("Ne/N vs. sexual reproductive variance  ",
                                V[k]/bar(k), "  — three populations")),
    subtitle = "Dotted horizontal line = Ne = 50 threshold at N = 40",
    x        = expression(V[k]/bar(k)~"(reproductive stage)"),
    y        = expression(N[e]/N),
    color    = NULL,
    caption  = "Ne_sexual_Y2000, monthly matrices, Ne_target = 50, census_N = 40."
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title      = element_text(face = "bold"),
    legend.position = "bottom",
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.4)
  )

## ----session-info-------------------------------------------------------------
sessionInfo()

