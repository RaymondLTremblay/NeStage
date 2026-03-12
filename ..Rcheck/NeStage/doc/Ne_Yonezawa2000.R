## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  echo    = TRUE,
  eval    = TRUE,
  message = FALSE,
  warning = FALSE,
  comment = "#>"
)

## ----load-nestage, include=FALSE----------------------------------------------
library(NeStage)

## ----helpers, message=FALSE, warning=FALSE------------------------------------
if (!requireNamespace("gt",      quietly = TRUE)) install.packages("gt")
if (!requireNamespace("popbio",  quietly = TRUE)) install.packages("popbio")
library(gt)
library(popbio)

# Helper: generation time — Yonezawa (2000) T^x iteration method
gen_time_yonezawa <- function(T_mat, F_mat, xmax = 500) {
  s    <- nrow(T_mat)
  Fvec <- as.numeric(F_mat[1, ])
  A    <- T_mat + F_mat
  ev   <- eigen(A)
  w    <- abs(Re(ev$vectors[, 1])); w <- w / sum(w)
  Tx   <- diag(s); num <- 0; den <- 0
  for (x in seq_len(xmax)) {
    Tx  <- T_mat %*% Tx
    l_x <- 0; m_x <- 0
    for (i in seq_len(s)) {
      l_x <- l_x + w[i] * sum(Tx[, i])
      m_x <- m_x + w[i] * sum(Fvec * Tx[, i])
    }
    num <- num + x * m_x * l_x
    den <- den + m_x * l_x
  }
  if (den == 0) stop("Generation time undefined: no reproductive output over xmax years.")
  num / den
}

# Helper: stage-weighted mean of squared survivals
over_u2_from_T <- function(T_mat, D) { u <- colSums(T_mat); sum(D * u^2) }

# Clonal-dominant formulas (Eq. 10 and Eq. 11) — for pedagogical Steps 4–5
Ny_over_N        <- function(over_u2)    1 / (1 - over_u2)
Ne_over_N_approx <- function(over_u2, L) 1 / ((1 - over_u2) * L)

## ----data-entry, message=FALSE, warning=FALSE---------------------------------
stages <- c("Stage 1 (one-leaf)", "Stage 2 (multileaf NF)", "Stage 3 (multileaf FL)")

# --- Population Miz ---
T_Miz <- matrix(c(
  0.789, 0.121, 0.054,
  0.007, 0.621, 0.335,
  0.001, 0.258, 0.611
), nrow = 3, byrow = TRUE,
   dimnames = list(paste0("to_", 1:3), paste0("from_", 1:3)))

D_obs_Miz <- c(0.935, 0.038, 0.027)
D_exp_Miz <- c(0.921, 0.046, 0.033)
F_vec_Miz <- c(0.055, 1.328, 2.398)
F_Miz     <- matrix(0, 3, 3, dimnames = dimnames(T_Miz))
F_Miz[1,] <- F_vec_Miz

# --- Population Nan ---
T_Nan <- matrix(c(
  0.748, 0.137, 0.138,
  0.006, 0.669, 0.374,
  0.001, 0.194, 0.488
), nrow = 3, byrow = TRUE,
   dimnames = list(paste0("to_", 1:3), paste0("from_", 1:3)))

D_obs_Nan <- c(0.958, 0.027, 0.015)
D_exp_Nan <- c(0.951, 0.034, 0.015)
F_vec_Nan <- c(0.138, 2.773, 5.016)
F_Nan     <- matrix(0, 3, 3, dimnames = dimnames(T_Nan))
F_Nan[1,] <- F_vec_Nan

## ----table2-gt, message=FALSE, warning=FALSE----------------------------------
tbl2 <- data.frame(
  Variable = c(
    "u11","u21","u31","u12","u22","u32","u13","u23","u33",
    "D1 (obs)","D2 (obs)","D3 (obs)",
    "D1 (exp)","D2 (exp)","D3 (exp)",
    "u·1","u·2","u·3",
    "F1","F2","F3",
    "r1","r2","r3","Sum r_i"
  ),
  Description = c(
    "Stage 1→1","Stage 1→2","Stage 1→3",
    "Stage 2→1","Stage 2→2","Stage 2→3",
    "Stage 3→1","Stage 3→2","Stage 3→3",
    "Frac. stage 1 (obs)","Frac. stage 2 (obs)","Frac. stage 3 (obs)",
    "Frac. stage 1 (exp)","Frac. stage 2 (exp)","Frac. stage 3 (exp)",
    "Survival stage 1","Survival stage 2","Survival stage 3",
    "Newborns/plant stage 1","Newborns/plant stage 2","Newborns/plant stage 3",
    "Newborn contrib. r1","Newborn contrib. r2","Newborn contrib. r3",
    "Recruitment rate"
  ),
  Miz = c(
    T_Miz[1,1],T_Miz[2,1],T_Miz[3,1],
    T_Miz[1,2],T_Miz[2,2],T_Miz[3,2],
    T_Miz[1,3],T_Miz[2,3],T_Miz[3,3],
    D_obs_Miz, D_exp_Miz, colSums(T_Miz), F_vec_Miz,
    F_vec_Miz * D_obs_Miz, sum(F_vec_Miz * D_obs_Miz)
  ),
  Nan = c(
    T_Nan[1,1],T_Nan[2,1],T_Nan[3,1],
    T_Nan[1,2],T_Nan[2,2],T_Nan[3,2],
    T_Nan[1,3],T_Nan[2,3],T_Nan[3,3],
    D_obs_Nan, D_exp_Nan, colSums(T_Nan), F_vec_Nan,
    F_vec_Nan * D_obs_Nan, sum(F_vec_Nan * D_obs_Nan)
  )
)

gt(tbl2) |>
  tab_header(
    title    = md("**Table 2** — Demographic and reproductive variables"),
    subtitle = md("*Fritillaria camtschatcensis*, Yonezawa et al. (2000)")
  ) |>
  tab_row_group(label = md("**Transition matrix** [u_ji]"),    rows = 1:9)  |>
  tab_row_group(label = md("**Stage fractions** D_i"),         rows = 10:15)|>
  tab_row_group(label = md("**Survival rates** u·i"),          rows = 16:18)|>
  tab_row_group(label = md("**Newborns per plant** F_i"),      rows = 19:21)|>
  tab_row_group(label = md("**Newborn contributions** r_i"),   rows = 22:25)|>
  cols_label(Variable = "Symbol", Description = "Definition",
             Miz = "Miz", Nan = "Nan") |>
  fmt_number(columns = c(Miz, Nan), decimals = 3) |>
  tab_style(style = cell_text(weight = "bold"),
            locations = cells_row_groups()) |>
  tab_options(table.width = pct(85), table.font.size = px(13))

## ----eigen-and-L, message=FALSE, warning=FALSE--------------------------------
A_Miz <- T_Miz + F_Miz
A_Nan <- T_Nan + F_Nan

eg_Miz <- popbio::eigen.analysis(A_Miz)
eg_Nan <- popbio::eigen.analysis(A_Nan)

# L: authoritative values from Table 4 of Yonezawa et al. (2000)
# used in all replication calculations below.
L_Miz <- 13.399
L_Nan <-  8.353

# Exploratory: compare to computational approximations
L_Miz_Y <- gen_time_yonezawa(T_Miz, F_Miz, xmax = 500)
L_Nan_Y <- gen_time_yonezawa(T_Nan, F_Nan, xmax = 500)
L_Miz_P <- as.numeric(popbio::generation.time(A_Miz))
L_Nan_P <- as.numeric(popbio::generation.time(A_Nan))

## ----eigen-gt, message=FALSE, warning=FALSE-----------------------------------
data.frame(
  Population = c("Miz", "Nan"),
  lambda     = round(c(eg_Miz$lambda1, eg_Nan$lambda1), 4),
  L_paper    = c(L_Miz, L_Nan),
  L_Yonezawa = round(c(L_Miz_Y, L_Nan_Y), 3),
  L_popbio   = round(c(L_Miz_P, L_Nan_P), 3)
) |>
  gt() |>
  tab_header(
    title    = md("**Eigenvalues and generation times**"),
    subtitle = md("*L* (paper) = Table 4 values used for replication; others shown for comparison only")
  ) |>
  cols_label(
    Population = "Population", lambda = md("*λ*"),
    L_paper    = md("*L* (paper)"),
    L_Yonezawa = md("*L* (T^x iteration)"),
    L_popbio   = md("*L* (popbio / Caswell)")
  ) |>
  tab_spanner(label = "Generation time L — exploratory comparison",
              columns = c(L_Yonezawa, L_popbio)) |>
  tab_style(style = cell_fill(color = "#e8f5e9"),
            locations = cells_body(columns = L_paper)) |>
  tab_style(style = cell_fill(color = "#f0f7fb"),
            locations = cells_column_spanners()) |>
  tab_options(table.width = pct(80))

## ----compute-table4-observed, message=FALSE, warning=FALSE--------------------
ou2_Miz_obs <- over_u2_from_T(T_Miz, D_obs_Miz)
ou2_Nan_obs <- over_u2_from_T(T_Nan, D_obs_Nan)

NyN_Miz_obs  <- Ny_over_N(ou2_Miz_obs)
NyN_Nan_obs  <- Ny_over_N(ou2_Nan_obs)
NeN_Miz_obs  <- Ne_over_N_approx(ou2_Miz_obs, L_Miz)
NeN_Nan_obs  <- Ne_over_N_approx(ou2_Nan_obs, L_Nan)
MinN_Miz_obs <- ceiling(5000 / NeN_Miz_obs)
MinN_Nan_obs <- ceiling(5000 / NeN_Nan_obs)

## ----compute-parenthetical, message=FALSE, warning=FALSE----------------------
ou2_Miz_exp <- over_u2_from_T(T_Miz, D_exp_Miz)
ou2_Nan_exp <- over_u2_from_T(T_Nan, D_exp_Nan)

NyN_Miz_exp  <- Ny_over_N(ou2_Miz_exp)
NyN_Nan_exp  <- Ny_over_N(ou2_Nan_exp)
NeN_Miz_exp  <- Ne_over_N_approx(ou2_Miz_exp, L_Miz)
NeN_Nan_exp  <- Ne_over_N_approx(ou2_Nan_exp, L_Nan)
MinN_Miz_exp <- ceiling(5000 / NeN_Miz_exp)
MinN_Nan_exp <- ceiling(5000 / NeN_Nan_exp)

## ----replication-check, message=FALSE, warning=FALSE--------------------------
pop  <- c("Miz", "Nan")

# Paper values (Table 4, Yonezawa et al. 2000)
pL   <- c(13.399,  8.353)
pNyO <- c( 2.932,  2.428);  pNyE <- c(2.977, 2.444)
pNeO <- c( 0.219,  0.291);  pNeE <- c(0.222, 0.293)
pMnO <- c(22832,  17183);   pMnE <- c(22523, 17065)

cL   <- round(c(L_Miz, L_Nan), 3)
cNyO <- round(c(NyN_Miz_obs, NyN_Nan_obs), 3)
cNyE <- round(c(NyN_Miz_exp, NyN_Nan_exp), 3)
cNeO <- round(c(NeN_Miz_obs, NeN_Nan_obs), 3)
cNeE <- round(c(NeN_Miz_exp, NeN_Nan_exp), 3)
cMnO <- round(c(MinN_Miz_obs, MinN_Nan_obs))
cMnE <- round(c(MinN_Miz_exp, MinN_Nan_exp))

all_pass <- abs(cL - pL) <= 0.005 &
            abs(cNyO - pNyO) <= 0.002 & abs(cNyE - pNyE) <= 0.003 &
            abs(cNeO - pNeO) <= 0.002 & abs(cNeE - pNeE) <= 0.002 &
            abs(cMnO - pMnO) <= 25    & abs(cMnE - pMnE) <= 25

pass_rows <- which(all_pass)
fail_rows <- which(!all_pass)

# ---- Table A: L and Ny/N ----
data.frame(
  Pop  = pop,
  pL   = pL,   cL   = cL,
  pNyO = pNyO, cNyO = cNyO,
  pNyE = pNyE, cNyE = cNyE
) |>
  gt() |>
  tab_header(
    title    = md("**Replication check (Part 1 of 2): *L* and *N*_y/*N***"),
    subtitle = md("*Fritillaria camtschatcensis* — computed vs. Table 4, Yonezawa et al. (2000). Green = match within rounding.")
  ) |>
  tab_spanner(label = md("Generation time *L*"),         columns = c(pL, cL))   |>
  tab_spanner(label = md("*N*_y/*N*  (observed *D*)"),   columns = c(pNyO, cNyO)) |>
  tab_spanner(label = md("*N*_y/*N*  (expected *D*)"),   columns = c(pNyE, cNyE)) |>
  cols_label(Pop = "Pop",
             pL = "Paper", cL = "Computed",
             pNyO = "Paper", cNyO = "Computed",
             pNyE = "Paper", cNyE = "Computed") |>
  fmt_number(columns = c(pL, cL, pNyO, cNyO, pNyE, cNyE), decimals = 3) |>
  tab_style(style = cell_fill(color = "#e8f5e9"),
            locations = cells_body(rows = pass_rows)) |>
  tab_style(style = cell_fill(color = "#fff3e0"),
            locations = cells_body(rows = fail_rows)) |>
  tab_style(style = cell_fill(color = "#f0f7fb"),
            locations = cells_column_spanners()) |>
  tab_options(table.width = pct(72), table.font.size = px(13))

# ---- Table B: Ne/N and Min N ----
data.frame(
  Pop  = pop,
  pNeO = pNeO, cNeO = cNeO,
  pNeE = pNeE, cNeE = cNeE,
  pMnO = pMnO, cMnO = cMnO,
  pMnE = pMnE, cMnE = cMnE
) |>
  gt() |>
  tab_header(
    title    = md("**Replication check (Part 2 of 2): *N*_e/*N* and minimum *N***"),
    subtitle = md("*Fritillaria camtschatcensis* — computed vs. Table 4, Yonezawa et al. (2000). Green = match within rounding.")
  ) |>
  tab_spanner(label = md("*N*_e/*N*  (observed *D*)"),           columns = c(pNeO, cNeO)) |>
  tab_spanner(label = md("*N*_e/*N*  (expected *D*)"),           columns = c(pNeE, cNeE)) |>
  tab_spanner(label = md("Min *N* for *N*_e = 5000 (obs *D*)"), columns = c(pMnO, cMnO)) |>
  tab_spanner(label = md("Min *N* for *N*_e = 5000 (exp *D*)"), columns = c(pMnE, cMnE)) |>
  cols_label(Pop = "Pop",
             pNeO = "Paper", cNeO = "Computed",
             pNeE = "Paper", cNeE = "Computed",
             pMnO = "Paper", cMnO = "Computed",
             pMnE = "Paper", cMnE = "Computed") |>
  fmt_number(columns = c(pNeO, cNeO, pNeE, cNeE), decimals = 3) |>
  fmt_integer(columns = c(pMnO, cMnO, pMnE, cMnE)) |>
  tab_style(style = cell_fill(color = "#e8f5e9"),
            locations = cells_body(rows = pass_rows)) |>
  tab_style(style = cell_fill(color = "#fff3e0"),
            locations = cells_body(rows = fail_rows)) |>
  tab_style(style = cell_fill(color = "#f0f7fb"),
            locations = cells_column_spanners()) |>
  tab_options(table.width = pct(85), table.font.size = px(13))

## ----clonal-inputs-explained, eval=FALSE--------------------------------------
# Ne_clonal_Y2000(
#   T_mat      = T_Miz,    # MatU: survival-transition matrix (s x s)
#                           # Column sums must be <= 1
#   F_vec      = F_vec_Miz, # Per-capita clonal offspring production per stage
#                           # (row 1 of MatF, or a numeric vector of length s)
#   D          = D_obs_Miz, # Stage frequency vector, sums to 1
#   L          = L_Miz,    # Generation time (years). If NULL, computed
#                           # internally via the Yonezawa T^x iteration.
#                           # Supply the paper value here for exact replication.
#   Ne_target  = 5000,     # Ne viability threshold (Lande 1995).
#                           # For Fritillaria (widespread alpine species) the
#                           # long-term evolutionary threshold is appropriate.
#                           # For small endemic species use 50 (Franklin 1980).
#   census_N   = 200,      # Your actual or expected census population size.
#                           # Reports Ne_at_census = NeN * census_N directly.
#   population = "Miz"     # Label for printed output
# )

## ----clonal-observed, message=FALSE, warning=FALSE----------------------------
# --- Observed stage fractions ---
Ne_Miz_obs <- Ne_clonal_Y2000(
  T_mat      = T_Miz,
  F_vec      = F_vec_Miz,
  D          = D_obs_Miz,
  L          = L_Miz,
  Ne_target  = 5000,
  census_N   = 200,
  population = "Miz (observed D)"
)

Ne_Nan_obs <- Ne_clonal_Y2000(
  T_mat      = T_Nan,
  F_vec      = F_vec_Nan,
  D          = D_obs_Nan,
  L          = L_Nan,
  Ne_target  = 5000,
  census_N   = 200,
  population = "Nan (observed D)"
)

# --- Expected (equilibrium) stage fractions ---
Ne_Miz_exp <- Ne_clonal_Y2000(
  T_mat      = T_Miz,
  F_vec      = F_vec_Miz,
  D          = D_exp_Miz,
  L          = L_Miz,
  Ne_target  = 5000,
  census_N   = 200,
  population = "Miz (expected D)"
)

Ne_Nan_exp <- Ne_clonal_Y2000(
  T_mat      = T_Nan,
  F_vec      = F_vec_Nan,
  D          = D_exp_Nan,
  L          = L_Nan,
  Ne_target  = 5000,
  census_N   = 200,
  population = "Nan (expected D)"
)

## ----print-Miz-obs------------------------------------------------------------
print(Ne_Miz_obs)

## ----nestage-vs-paper, message=FALSE, warning=FALSE---------------------------
data.frame(
  Population = c("Miz", "Nan", "Miz", "Nan"),
  D_type     = c("Observed", "Observed", "Expected", "Expected"),
  L_used     = c(L_Miz, L_Nan, L_Miz, L_Nan),
  NyN_NeStage = round(c(Ne_Miz_obs$NyN, Ne_Nan_obs$NyN,
                         Ne_Miz_exp$NyN, Ne_Nan_exp$NyN), 3),
  NyN_paper   = c(2.932, 2.428, 2.977, 2.444),
  NeN_NeStage = round(c(Ne_Miz_obs$NeN, Ne_Nan_obs$NeN,
                         Ne_Miz_exp$NeN, Ne_Nan_exp$NeN), 3),
  NeN_paper   = c(0.219, 0.291, 0.222, 0.293),
  MinN_NeStage = c(Ne_Miz_obs$Min_N, Ne_Nan_obs$Min_N,
                   Ne_Miz_exp$Min_N, Ne_Nan_exp$Min_N),
  MinN_paper   = c(22832, 17183, 22523, 17065)
) |>
  gt(groupname_col = "D_type") |>
  tab_header(
    title    = md("**`Ne_clonal_Y2000()` vs. Table 4 — Yonezawa et al. (2000)**"),
    subtitle = md("Parenthetical values in the paper correspond to rows using *Expected* stage fractions")
  ) |>
  cols_label(
    Population   = "Population",
    L_used       = md("*L* (years)"),
    NyN_NeStage  = md("*N*_y/*N* (NeStage)"),
    NyN_paper    = md("*N*_y/*N* (paper)"),
    NeN_NeStage  = md("*N*_e/*N* (NeStage)"),
    NeN_paper    = md("*N*_e/*N* (paper)"),
    MinN_NeStage = md("Min *N* (NeStage)"),
    MinN_paper   = md("Min *N* (paper)")
  ) |>
  tab_spanner(label = md("*N*_y/*N*"),  columns = c(NyN_NeStage, NyN_paper))  |>
  tab_spanner(label = md("*N*_e/*N*"),  columns = c(NeN_NeStage, NeN_paper))  |>
  tab_spanner(label = md("Min *N* for *N*_e ≥ 5000"),
              columns = c(MinN_NeStage, MinN_paper)) |>
  fmt_number(columns = c(NyN_NeStage, NyN_paper, NeN_NeStage, NeN_paper),
             decimals = 3) |>
  fmt_integer(columns = c(MinN_NeStage, MinN_paper)) |>
  tab_style(
    style     = list(cell_fill(color = "#f5f5f5"), cell_text(weight = "bold")),
    locations = cells_row_groups()
  ) |>
  tab_footnote(
    footnote  = "Min N = minimum census size required for Ne >= 5000 (Lande 1995).",
    locations = cells_column_spanners(spanners = "Min *N* for *N*_e ≥ 5000")
  ) |>
  tab_options(table.width = pct(90))

## ----mixed-consistency, message=FALSE, warning=FALSE--------------------------
# Run the general model with clonal-dominant Poisson defaults
Ne_Miz_gen <- Ne_mixed_Y2000(
  T_mat      = T_Miz,
  F_vec      = F_vec_Miz,
  D          = D_obs_Miz,
  d          = rep(1, 3),        # fully clonal
  Vc_over_c  = rep(1, 3),        # Poisson clonal variance
  Vk_over_k  = rep(1, 3),        # Poisson sexual variance (irrelevant when d=1)
  a          = 0,
  L          = L_Miz,
  Ne_target  = 5000,
  census_N   = 200,
  population = "Miz (general model, clonal defaults)"
)

Ne_Nan_gen <- Ne_mixed_Y2000(
  T_mat      = T_Nan,
  F_vec      = F_vec_Nan,
  D          = D_obs_Nan,
  d          = rep(1, 3),
  Vc_over_c  = rep(1, 3),
  Vk_over_k  = rep(1, 3),
  a          = 0,
  L          = L_Nan,
  Ne_target  = 5000,
  census_N   = 200,
  population = "Nan (general model, clonal defaults)"
)

data.frame(
  Population  = c("Miz", "Nan"),
  NeN_clonal  = round(c(Ne_Miz_obs$NeN, Ne_Nan_obs$NeN), 6),
  NeN_general = round(c(Ne_Miz_gen$NeN, Ne_Nan_gen$NeN), 6),
  Difference  = formatC(
    c(abs(Ne_Miz_gen$NeN - Ne_Miz_obs$NeN),
      abs(Ne_Nan_gen$NeN - Ne_Nan_obs$NeN)),
    format = "e", digits = 2)
) |>
  gt() |>
  tab_header(
    title    = md("**Internal consistency check**"),
    subtitle = md("`Ne_mixed_Y2000()` with clonal defaults must equal `Ne_clonal_Y2000()`")
  ) |>
  cols_label(
    Population  = "Population",
    NeN_clonal  = md("*N*_e/*N* (`Ne_clonal_Y2000`)"),
    NeN_general = md("*N*_e/*N* (`Ne_mixed_Y2000`, clonal defaults)"),
    Difference  = "Difference"
  ) |>
  tab_footnote("Difference should be < 1e-10 under clonal-dominant Poisson defaults.") |>
  tab_options(table.width = pct(70))

## ----mixed-sensitivity, message=FALSE, warning=FALSE--------------------------
# Vary d from 0 (fully sexual) to 1 (fully clonal) under Poisson defaults
d_vals <- c(0, 0.25, 0.5, 0.75, 1.0)

sensitivity_d <- purrr::map_dfr(d_vals, function(d_val) {
  r <- Ne_mixed_Y2000(
    T_mat     = T_Miz,
    F_vec     = F_vec_Miz,
    D         = D_obs_Miz,
    d         = rep(d_val, 3),
    Vc_over_c = rep(1, 3),
    Vk_over_k = rep(1, 3),
    a         = 0,
    L         = L_Miz,
    Ne_target = 5000,
    population = paste0("Miz d=", d_val)
  )
  data.frame(d = d_val, NeN = round(r$NeN, 6), V = round(r$V, 8))
})

# Now vary Vc_over_c under d=1 (clonal) to show when it matters
vc_vals <- c(1, 2, 5, 10)

sensitivity_Vc <- purrr::map_dfr(vc_vals, function(vc) {
  r <- Ne_mixed_Y2000(
    T_mat     = T_Miz,
    F_vec     = F_vec_Miz,
    D         = D_obs_Miz,
    d         = rep(1, 3),
    Vc_over_c = rep(vc, 3),
    Vk_over_k = rep(1, 3),
    a         = 0,
    L         = L_Miz,
    Ne_target = 5000,
    population = paste0("Miz Vc/c=", vc)
  )
  data.frame(Vc_over_c = vc, NeN = round(r$NeN, 4), V = round(r$V, 6))
})

# Display both sensitivity tables
gt(sensitivity_d) |>
  tab_header(
    title    = md("**Effect of clonal fraction *d* on *N*_e/*N* (Miz, Poisson defaults)**"),
    subtitle = md("Under Poisson variance defaults, *d* has no effect on *N*_e/*N*")
  ) |>
  cols_label(d = md("*d* (clonal fraction)"),
             NeN = md("*N*_e/*N*"),
             V = "V (total variance)") |>
  tab_footnote("V and Ne/N are identical across all values of d under Poisson defaults.") |>
  tab_options(table.width = pct(55))

gt(sensitivity_Vc) |>
  tab_header(
    title    = md("**Effect of super-Poisson clonal variance on *N*_e/*N* (Miz, d=1)**"),
    subtitle = md("When clonal output is more variable than Poisson, Ne/N decreases substantially")
  ) |>
  cols_label(Vc_over_c = md("*V*_c/*c̄* (clonal variance/mean)"),
             NeN = md("*N*_e/*N*"),
             V = "V (total variance)") |>
  tab_footnote(
    md("*V*_c/*c̄* = 1 is Poisson (baseline); > 1 means some genets dominate clonal output,
       increasing genetic drift and reducing *N*_e/*N*.")
  ) |>
  tab_options(table.width = pct(55))

## ----final-summary, message=FALSE, warning=FALSE------------------------------
data.frame(
  Population  = c("Miz", "Nan", "Miz", "Nan"),
  D_type      = c("Observed", "Observed", "Expected", "Expected"),
  L           = c(L_Miz, L_Nan, L_Miz, L_Nan),
  NyN         = round(c(Ne_Miz_obs$NyN, Ne_Nan_obs$NyN,
                         Ne_Miz_exp$NyN, Ne_Nan_exp$NyN), 3),
  NeN         = round(c(Ne_Miz_obs$NeN, Ne_Nan_obs$NeN,
                         Ne_Miz_exp$NeN, Ne_Nan_exp$NeN), 3),
  Ne_at_N200  = round(c(Ne_Miz_obs$NeN, Ne_Nan_obs$NeN,
                         Ne_Miz_exp$NeN, Ne_Nan_exp$NeN) * 200, 1),
  MinN        = c(Ne_Miz_obs$Min_N, Ne_Nan_obs$Min_N,
                  Ne_Miz_exp$Min_N, Ne_Nan_exp$Min_N)
) |>
  gt(groupname_col = "D_type") |>
  tab_header(
    title    = md("**Full replication of Table 4 — Yonezawa et al. (2000)**"),
    subtitle = md("All values computed by `Ne_clonal_Y2000()` with paper *L* values")
  ) |>
  cols_label(
    Population  = "Population",
    L           = md("*L* (years)"),
    NyN         = md("*N*_y/*N*"),
    NeN         = md("*N*_e/*N*"),
    Ne_at_N200  = md("*N*_e at *N* = 200"),
    MinN        = md("Min *N* for *N*_e ≥ 5000")
  ) |>
  fmt_number(columns = c(L, NyN, NeN), decimals = 3) |>
  fmt_number(columns = Ne_at_N200, decimals = 1)     |>
  fmt_integer(columns = MinN)                        |>
  tab_style(
    style     = list(cell_fill(color = "#f5f5f5"), cell_text(weight = "bold")),
    locations = cells_row_groups()
  ) |>
  tab_footnote(
    footnote  = "Ne at N=200 = NeN * 200; illustrates the actual Ne achieved by a
                 census population of 200 individuals.",
    locations = cells_column_labels(columns = Ne_at_N200)
  ) |>
  tab_footnote(
    footnote  = "Min N = minimum census size for Ne >= 5000 (Lande 1995).
                 For Fritillaria, a widespread alpine species, the long-term
                 evolutionary threshold is the appropriate benchmark.",
    locations = cells_column_labels(columns = MinN)
  ) |>
  tab_options(table.width = pct(85))

## ----session-info-------------------------------------------------------------
sessionInfo()

