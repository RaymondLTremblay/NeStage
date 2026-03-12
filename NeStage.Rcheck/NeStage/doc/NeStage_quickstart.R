## ----setup, include=FALSE-----------------------------------------------------

knitr::opts_chunk$set(
  echo    = TRUE,
  eval    = TRUE,
  message = FALSE,
  warning = FALSE
)
library(NeStage)

## ----tmat-example, eval=FALSE-------------------------------------------------
# # Example: 3-stage plant (seedling, juvenile, adult)
# T_mat <- matrix(c(
#   0.30, 0.05, 0.00,   # row 1: transitions INTO stage 1
#   0.40, 0.65, 0.10,   # row 2: transitions INTO stage 2
#   0.00, 0.20, 0.80    # row 3: transitions INTO stage 3
# ), nrow = 3, byrow = TRUE)
# 
# # Quick check: column sums must all be <= 1
# colSums(T_mat)   # should all be <= 1.0

## ----fvec-example, eval=FALSE-------------------------------------------------
# # Example: only adults (stage 3) reproduce
# F_vec <- c(0.0, 0.5, 3.0)

## ----dvec-example, eval=FALSE-------------------------------------------------
# # Example: from field counts
# counts <- c(120, 50, 30)
# D <- counts / sum(counts)
# D        # c(0.60, 0.25, 0.15)
# sum(D)   # must equal 1

## ----clonal-example-----------------------------------------------------------
library(NeStage)

# --- Your data goes here ---
T_mat <- matrix(c(
  0.789, 0.121, 0.054,
  0.007, 0.621, 0.335,
  0.001, 0.258, 0.611
), nrow = 3, byrow = TRUE)

F_vec <- c(0.055, 1.328, 2.398)   # clonal propagules per stage per year
D     <- c(0.935, 0.038, 0.027)   # observed stage fractions (must sum to 1)
L     <- 13.399                    # generation time (years); omit to compute

# --- Run the function ---
result <- Ne_clonal_Y2000(
  T_mat      = T_mat,
  F_vec      = F_vec,
  D          = D,
  L          = L,
  Ne_target  = 5000,        # minimum Ne for long-term viability (Lande 1995)
  population = "My population"
)

print(result)

## ----sexual-example-----------------------------------------------------------
# --- Your data goes here ---
T_mat <- matrix(c(
  0.30, 0.05, 0.00,
  0.40, 0.65, 0.10,
  0.00, 0.20, 0.80
), nrow = 3, byrow = TRUE)

F_vec <- c(0.0, 0.5, 3.0)      # seeds/offspring per individual per stage
D     <- c(0.60, 0.25, 0.15)   # observed stage fractions (must sum to 1)

# --- Run the function ---
result_sex <- Ne_sexual_Y2000(
  T_mat      = T_mat,
  F_vec      = F_vec,
  D          = D,
  # L omitted: computed internally from T_mat and F_vec
  Ne_target  = 5000,
  population = "My sexual population"
)

print(result_sex)

## ----sexual-vk, eval=FALSE----------------------------------------------------
# # Stage 3 adults have high reproductive variance (pollinator-limited)
# # Stages 1 and 2 kept at Poisson default (= 1)
# result_highvar <- Ne_sexual_Y2000(
#   T_mat     = T_mat,
#   F_vec     = F_vec,
#   D         = D,
#   Vk_over_k = c(1, 1, 3),   # stage 3: Vk/k_bar = 3
#   Ne_target = 5000
# )
# print(result_highvar)
# # Ne/N will be lower than the Poisson default above

## ----mixed-example------------------------------------------------------------
# --- Your data goes here ---
T_mat <- matrix(c(
  0.30, 0.05, 0.00,
  0.40, 0.65, 0.10,
  0.00, 0.20, 0.80
), nrow = 3, byrow = TRUE)

F_vec <- c(0.0, 0.5, 3.0)      # total fecundity per stage (sexual + clonal)
D     <- c(0.60, 0.25, 0.15)   # observed stage fractions (must sum to 1)

# d: clonal fraction per stage
#   d_i = 0   -> stage i reproduces entirely sexually
#   d_i = 1   -> stage i reproduces entirely clonally
#   d_i = 0.7 -> 70% of stage i reproduction is clonal
d <- c(0.0, 0.0, 0.7)   # only adults (stage 3) reproduce clonally

# --- Run the function ---
result_mix <- Ne_mixed_Y2000(
  T_mat      = T_mat,
  F_vec      = F_vec,
  D          = D,
  d          = d,
  # Vk_over_k and Vc_over_c default to rep(1, s) -- Poisson variance
  Ne_target  = 5000,
  population = "My mixed population"
)

print(result_mix)

## ----d-examples, eval=FALSE---------------------------------------------------
# # All stages reproduce sexually only (equivalent to Ne_sexual_Y2000)
# d <- c(0, 0, 0)
# 
# # All stages reproduce clonally only
# d <- c(1, 1, 1)
# 
# # Only the largest stage (stage 3) has clonal reproduction
# d <- c(0, 0, 1)
# 
# # Gradual increase in clonal fraction across stages
# d <- c(0.1, 0.4, 0.8)

## ----session-info-------------------------------------------------------------
sessionInfo()

