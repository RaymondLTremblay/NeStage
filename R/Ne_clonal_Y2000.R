# =============================================================================
# Ne_clonal_Y2000.R
# -----------------------------------------------------------------------------
# Variance effective population size for 100% CLONAL stage-structured
# populations, following Yonezawa et al. (2000).
#
# Reference:
#   Yonezawa K., Kinoshita E., Watano Y., and Zentools::showNonASCIIfile("R/Ne_clonal_Y2000.R")toh H. (2000).
#   Formulation and estimation of the effective size of stage-structured
#   populations in Fritillaria camtschatcensis, a perennial herb with a
#   complex life history. Evolution 54(6): 2007-2013.
#   https://doi.org/10.1111/j.0014-3820.2000.tb01243.x
#
# When to use this function:
#   Use Ne_clonal_Y2000() when your study population reproduces EXCLUSIVELY
#   through clonal means (e.g., bulblets, stolons, rhizomes, vegetative
#   fragmentation) and no sexual reproduction contributes meaningfully to
#   recruitment. This is the case for Fritillaria camtschatcensis, the
#   validation species in the paper.
#
#   If your species has ANY sexual reproduction contributing to recruitment,
#   use Ne_mixed_Y2000() instead.
#
# Model assumptions:
#   1. The population is at demographic equilibrium: stage fractions D_i and
#      transition rates u_ji are constant over time.
#   2. All reproduction is clonal (d_i = 1 for all stages i).
#   3. Clonal reproductive output follows a Poisson distribution:
#      (V_c / c_bar)_i = 1 for all stages i.
#      This is the assumption that simplifies the general Eq. 6 to Eq. 10.
#   4. The population mates randomly (or does not mate sexually), so the
#      Hardy-Weinberg deviation a = 0.
#   5. Census size N is approximately constant over the generation time L.
#
# Under these assumptions, Eq. 6 of Yonezawa (2000) simplifies to:
#
#   Ny/N = 1 / (1 - u2_bar)          ... Eq. 11  (annual effective size)
#   Ne/N = 1 / ((1 - u2_bar) * L)    ... Eq. 10  (generation-time effective size)
#
# where:
#   u2_bar = sum_i [ D_i * (u_{.i})^2 ]
#   u_{.i} = sum_j u_{ji}  (total annual survival rate of stage i)
#   L      = generation time (mean age of reproduction, see details below)
#
# IMPORTANT: u2_bar vs u_bar^2
#   These are NOT the same quantity.
#   u_bar   = sum_i D_i * u_{.i}    (weighted mean of survivals)
#   u2_bar  = sum_i D_i * u_{.i}^2  (weighted mean of SQUARED survivals)
#   u2_bar >= u_bar^2 always (Jensen's inequality).
#   The difference (u2_bar - u_bar^2) captures between-stage variance in
#   survival and is the primary driver of genetic drift in this model.
#   Using u_bar^2 instead of u2_bar would be a mathematical error.
#
# Generation time L:
#   L is the population-average mean age of reproduction. It is computed by
#   iterating the survival matrix T to age x_max = 500, then averaging
#   cohort statistics across all starting stages weighted by the stable
#   stage distribution w (dominant right eigenvector of A = T + F_mat).
#   See .compute_L_clonal() below for the implementation.
#   Alternatively, the user may supply L directly from the published source
#   (recommended for replication of specific papers).
#
# Author:  Raymond L. Tremblay
# Version: 1.0.0  (2026-03-05)
# =============================================================================

# -----------------------------------------------------------------------------
# SECTION 1: Input validation
# -----------------------------------------------------------------------------
# These helpers check that the inputs are internally consistent before any
# computation begins. Clear error messages help users diagnose data problems.

.validate_T_mat <- function(T_mat) {
  # T_mat must be a square numeric matrix with no negative or non-finite values.
  # Each column represents one stage; each entry u_{ji} is a transition rate.
  # Column sums must be <= 1 (total annual survival cannot exceed 1).
  if (!is.matrix(T_mat) || !is.numeric(T_mat)) {
    stop("T_mat must be a numeric matrix.")
  }
  if (nrow(T_mat) != ncol(T_mat)) {
    stop("T_mat must be square (nrow == ncol).")
  }
  if (any(!is.finite(T_mat))) {
    stop("T_mat contains non-finite values (NA, NaN, Inf).")
  }
  if (any(T_mat < 0)) {
    stop("T_mat contains negative values. All transition rates must be >= 0.")
  }
  col_sums <- colSums(T_mat)
  if (any(col_sums > 1 + 1e-8)) {
    stop(paste0(
      "Column sums of T_mat exceed 1.0. Total annual survival cannot exceed 1.\n",
      "  Column sums: ",
      paste(round(col_sums, 4), collapse = ", ")
    ))
  }
  invisible(TRUE)
}

.validate_D <- function(D, s) {
  # D must be a non-negative numeric vector of length s that sums to 1.
  # It represents the fraction of the population in each stage.
  if (!is.numeric(D) || length(D) != s) {
    stop(paste0(
      "D must be a numeric vector of length ",
      s,
      " (one entry per stage)."
    ))
  }
  if (any(!is.finite(D)) || any(D < 0)) {
    stop("D must contain finite non-negative values.")
  }
  if (abs(sum(D) - 1) > 1e-6) {
    stop(paste0(
      "D must sum to 1. Current sum = ",
      round(sum(D), 6),
      ".\n  Hint: D <- D / sum(D) to normalise."
    ))
  }
  invisible(TRUE)
}

.validate_F_vec <- function(F_vec, s) {
  # F_vec is the fecundity vector: newborns produced per individual per stage
  # per year. For a purely clonal population these are clonal propagules.
  # All entries must be >= 0; at least one must be > 0 for L computation.
  if (!is.numeric(F_vec) || length(F_vec) != s) {
    stop(paste0(
      "F_vec must be a numeric vector of length ",
      s,
      " (one fecundity per stage)."
    ))
  }
  if (any(!is.finite(F_vec)) || any(F_vec < 0)) {
    stop("F_vec must contain finite non-negative values.")
  }
  if (all(F_vec == 0)) {
    stop(
      "F_vec is all zeros: no reproductive output in any stage. Cannot compute L."
    )
  }
  invisible(TRUE)
}

.validate_L <- function(L) {
  # If L is supplied by the user it must be a single positive finite number.
  if (!is.numeric(L) || length(L) != 1) {
    stop("L must be a single numeric value (generation time in years).")
  }
  if (!is.finite(L) || L <= 0) {
    stop(paste0("L must be a positive finite number. Received: ", L))
  }
  invisible(TRUE)
}


# -----------------------------------------------------------------------------
# SECTION 2: Generation time L
# -----------------------------------------------------------------------------
# L is the population-average mean age of reproduction. This section
# implements the Yonezawa (2000) definition (paper p. 2008):
#
#   L = sum_{x=1}^{x_max} x * m_bar_x * l_bar_x
#       / sum_{x=1}^{x_max} m_bar_x * l_bar_x
#
# where the population-average survivorship and reproduction at age x are:
#
#   l_bar_x = sum_i w_i * sum_j T[j, i]^x   (probability of surviving to age x)
#   m_bar_x = sum_i w_i * sum_j F_j * T[j, i]^x / l_bar_x
#             (mean reproductive output at age x, averaged over stages)
#
# T[j,i]^x  = element (j,i) of T^x = the survival matrix raised to power x,
#              i.e. the probability that an individual starting in stage i
#              at age 0 is found in stage j at age x.
#
# The weighting by the stable stage distribution w (dominant right eigenvector
# of A = T + F_mat) is essential: it gives each starting stage its
# demographically appropriate weight in the population average.

.compute_L_clonal <- function(T_mat, F_vec, x_max = 500L) {
  # Compute generation time L using the Yonezawa (2000) definition.
  #
  # Args:
  #   T_mat : s x s survival/transition matrix
  #   F_vec : length-s fecundity vector (newborns per stage per year)
  #   x_max : maximum age to iterate to (default 500; paper uses 500)
  #
  # Returns: L as a single numeric value (years)

  s <- nrow(T_mat)
  F_mat <- matrix(0, s, s)
  F_mat[1, ] <- F_vec # fecundity goes into row 1 of the full matrix

  # Stable stage distribution: dominant right eigenvector of A = T + F_mat
  A <- T_mat + F_mat
  ev <- eigen(A)
  w <- Re(ev$vectors[, which.max(Re(ev$values))])
  w <- abs(w) / sum(abs(w)) # normalise to sum = 1

  # Iterate T^x from x = 1 to x_max, accumulating numerator and denominator
  Tx <- diag(s) # T^0 = identity matrix
  num <- 0 # accumulates sum(x * m_bar_x * l_bar_x)
  den <- 0 # accumulates sum(m_bar_x * l_bar_x)

  for (x in seq_len(x_max)) {
    Tx <- T_mat %*% Tx # left-multiply to advance one year: T^x = T * T^(x-1)

    # Population-average survivorship and fecundity at age x,
    # weighted by stable stage distribution w
    l_bar_x <- 0
    m_bar_x <- 0
    for (i in seq_len(s)) {
      u_jxi <- Tx[, i] # stage distribution at age x, given start in stage i
      l_bar_x <- l_bar_x + w[i] * sum(u_jxi)
      m_bar_x <- m_bar_x + w[i] * sum(F_vec * u_jxi)
    }

    num <- num + x * m_bar_x * l_bar_x
    den <- den + m_bar_x * l_bar_x
  }

  if (den <= 0) {
    stop(paste0(
      "Generation time L is undefined: no reproductive output detected over ",
      x_max,
      " years.\n",
      "  Check that F_vec has positive entries and T_mat allows survival."
    ))
  }

  as.numeric(num / den)
}


# -----------------------------------------------------------------------------
# SECTION 3: Core computation
# -----------------------------------------------------------------------------
# This section implements the key intermediate quantities and the final
# effective size ratios from Equations 10 and 11 of Yonezawa (2000).

.compute_u2_bar <- function(T_mat, D) {
  # Compute u2_bar = sum_i D_i * (u_{.i})^2
  # where u_{.i} = colSums(T_mat)[i] is the total annual survival of stage i.
  #
  # This is the stage-weighted mean of SQUARED survivals -- the central
  # intermediate quantity in Eqs. 10 and 11. A larger u2_bar means more
  # individuals survive and persist in stages from year to year, which
  # reduces the variance in allele-frequency change and increases Ne.
  u_dot <- colSums(T_mat) # u_{.i} for each stage i
  sum(D * u_dot^2) # weighted mean of squared survivals
}

.compute_Ne_clonal_core <- function(u2_bar, L) {
  # Compute the two effective size ratios from u2_bar and L.
  #
  # Eq. 11:  Ny/N = 1 / (1 - u2_bar)
  #   The ANNUAL effective size ratio: how large is Ne relative to census
  #   size N when measured per year rather than per generation.
  #   A larger u2_bar (more survival, less turnover) gives a larger Ny/N.
  #
  # Eq. 10:  Ne/N = 1 / ((1 - u2_bar) * L)
  #   The GENERATION-TIME effective size ratio: the standard conservation
  #   metric. Dividing by L accounts for the fact that in long-lived species
  #   a "generation" spans many years, diluting the per-year drift.
  #
  # Note: (1 - u2_bar) is the "effective annual turnover" -- the fraction
  # of the genetic pool that is replaced each year. It must be > 0 for
  # Ne to be finite (i.e., some mortality must occur somewhere).

  denom <- 1 - u2_bar

  if (!is.finite(denom) || denom <= 0) {
    stop(paste0(
      "1 - u2_bar = ",
      round(denom, 6),
      " is not positive.\n",
      "  This means u2_bar >= 1, implying no mortality occurs.\n",
      "  Check that at least one column sum of T_mat is < 1."
    ))
  }

  NyN <- 1 / denom # Eq. 11
  NeN <- 1 / (denom * L) # Eq. 10

  list(NyN = NyN, NeN = NeN)
}


# -----------------------------------------------------------------------------
# SECTION 4: Public function
# -----------------------------------------------------------------------------

#' Variance effective population size for 100% clonal stage-structured populations
#'
#' Implements Equations 10 and 11 of Yonezawa et al. (2000) for populations
#' that reproduce exclusively through clonal means (d_i = 1 for all stages)
#' with Poisson-distributed clonal output (V_c / c_bar = 1).
#'
#' @param T_mat  Numeric matrix (s x s). The survival/transition matrix.
#'               Row j, column i gives the probability that an individual
#'               in stage i at year t is in stage j at year t+1.
#'               Row/column order: stages 1, 2, ..., s (stage 1 = juvenile/smallest).
#'               Column sums must be <= 1 (total annual survival per stage).
#'
#' @param F_vec  Numeric vector of length s. Fecundity of each stage: the
#'               mean number of clonal propagules produced per individual per
#'               year. Used only to compute L (generation time) via the
#'               stable stage distribution of A = T_mat + F_mat.
#'               Set F_vec = NULL if supplying L directly.
#'
#' @param D      Numeric vector of length s. The stage frequency distribution:
#'               the proportion of the population in each stage. Must sum to 1.
#'               Use observed field counts or the stable-stage distribution
#'               from eigen analysis (they should be close at equilibrium).
#'
#' @param L      Numeric scalar (optional). Generation time in years.
#'               If NULL (default), L is computed internally from T_mat and
#'               F_vec using the Yonezawa (2000) definition (x_max = 500 yr).
#'               Supply L directly to replicate a specific published value
#'               (e.g., L = 13.399 for Fritillaria Miz population).
#'
#' @param x_max  Integer. Maximum age (years) for the L iteration (default 500).
#'               Ignored if L is supplied directly. Increase for very
#'               long-lived species (> 100 yr lifespan).
#'
#' @param Ne_target  Numeric. Ne conservation threshold. Common choices:
#'               50 (Franklin 1980, avoid inbreeding depression),
#'               500 (long-term quantitative variation), or
#'               5000 (Lande 1995, long-term evolutionary potential).
#'               Default 50.
#' @param census_N  Numeric or NULL. Actual or expected census population
#'               size. If supplied, Ne_at_census = NeN * census_N is
#'               added to the output. Default NULL.
#' @param population  Character string. Optional label for the population,
#'               used in printed output and the returned list.
#'
#' @return A named list with the following elements:
#'   \describe{
#'     \item{population}{Character label (from \code{population} argument)}
#'     \item{model}{"clonal_Y2000" -- identifies the model used}
#'     \item{NyN}{Ny/N -- annual effective size ratio (Eq. 11)}
#'     \item{NeN}{Ne/N -- generation-time effective size ratio (Eq. 10)}
#'     \item{L}{Generation time used in the calculation (years)}
#'     \item{L_source}{"user" if L was supplied; "computed" if derived internally}
#'     \item{u_dot}{Named vector of total annual survival per stage (colSums of T_mat)}
#'     \item{u2_bar}{Stage-weighted mean of squared survivals (key intermediate)}
#'     \item{Min_N}{Minimum census size N for Ne >= Ne_target (Lande 1995)}
#'     \item{Ne_target}{The Ne threshold used for Min_N (default 5000)}
#'   }
#'
#' @examples
#' # ---------------------------------------------------------
#' # Example 1: Replicate Yonezawa Table 4 for population Miz
#' # ---------------------------------------------------------
#' # Transition matrix from Table 2 of Yonezawa et al. (2000)
#' T_Miz <- matrix(c(
#'   0.789, 0.121, 0.054,
#'   0.007, 0.621, 0.335,
#'   0.001, 0.258, 0.611
#' ), nrow = 3, byrow = TRUE)
#'
#' # Observed stage fractions (Table 2)
#' D_Miz <- c(0.935, 0.038, 0.027)
#'
#' # Fecundity vector (Table 2) -- clonal propagules per plant per year
#' F_Miz <- c(0.055, 1.328, 2.398)
#'
#' # Using L from Table 4 directly (recommended for exact replication)
#' result <- Ne_clonal_Y2000(
#'   T_mat      = T_Miz,
#'   F_vec      = F_Miz,
#'   D          = D_Miz,
#'   L          = 13.399,
#'   population = "Miz"
#' )
#' print(result)
#' # Expected: Ny/N = 2.932, Ne/N = 0.219, Min_N = 22,831
#'
#' # ---------------------------------------------------------
#' # Example 2: Let the function compute L internally
#' # ---------------------------------------------------------
#' result2 <- Ne_clonal_Y2000(
#'   T_mat      = T_Miz,
#'   F_vec      = F_Miz,
#'   D          = D_Miz,
#'   population = "Miz (computed L)"
#' )
#' print(result2)
#'
#' # ---------------------------------------------------------
#' # Example 3: Replicate Table 4 for population Nan
#' # ---------------------------------------------------------
#' T_Nan <- matrix(c(
#'   0.748, 0.137, 0.138,
#'   0.006, 0.669, 0.374,
#'   0.001, 0.194, 0.488
#' ), nrow = 3, byrow = TRUE)
#'
#' D_Nan <- c(0.958, 0.027, 0.015)
#' F_Nan <- c(0.138, 2.773, 5.016)
#'
#' result3 <- Ne_clonal_Y2000(
#'   T_mat      = T_Nan,
#'   F_vec      = F_Nan,
#'   D          = D_Nan,
#'   L          = 8.353,
#'   population = "Nan"
#' )
#' print(result3)
#' # Expected: Ny/N = 2.428, Ne/N = 0.291, Min_N = 17,183
#'
#' @references
#' Yonezawa K., Kinoshita E., Watano Y., and Zentoh H. (2000).
#' Formulation and estimation of the effective size of stage-structured
#' populations in Fritillaria camtschatcensis, a perennial herb with a
#' complex life history. \emph{Evolution} \strong{54}(6): 2007-2013.
#' \doi{10.1111/j.0014-3820.2000.tb01244.x}
#'
#' Lande R. (1995). Mutation and conservation.
#' \emph{Conservation Biology} \strong{9}: 728-791.
#'
#' @export
Ne_clonal_Y2000 <- function(
  T_mat,
  F_vec = NULL,
  D,
  L = NULL,
  x_max = 500L,
  Ne_target = 50,
  census_N = NULL,
  population = NULL
) {
  # ------------------------------------------------------------------
  # Step 1: Validate all inputs before any computation
  # ------------------------------------------------------------------
  s <- nrow(T_mat)
  .validate_T_mat(T_mat)
  .validate_D(D, s)

  if (!is.null(F_vec)) {
    .validate_F_vec(F_vec, s)
  }
  if (!is.null(L)) {
    .validate_L(L)
  }

  if (is.null(L) && is.null(F_vec)) {
    stop(
      "Either L or F_vec must be supplied.\n",
      "  Provide L directly (e.g., from a published source), or\n",
      "  provide F_vec so L can be computed internally."
    )
  }

  if (!is.numeric(Ne_target) || Ne_target <= 0) {
    stop("Ne_target must be a positive number (e.g. 50, 500, or 5000).")
  }
  if (!is.null(census_N) && (!is.numeric(census_N) || census_N <= 0)) {
    stop("census_N must be a positive number if supplied.")
  }

  # ------------------------------------------------------------------
  # Step 2: Compute or retrieve generation time L
  # ------------------------------------------------------------------
  if (!is.null(L)) {
    L_use <- as.numeric(L)
    L_source <- "user" # user supplied L directly
  } else {
    L_use <- .compute_L_clonal(T_mat, F_vec, x_max = as.integer(x_max))
    L_source <- "computed" # computed internally via T^x iteration
  }

  # ------------------------------------------------------------------
  # Step 3: Compute u2_bar (stage-weighted mean of squared survivals)
  #
  # u_{.i} = colSums(T_mat)[i] = total annual survival of stage i
  # u2_bar = sum_i D_i * u_{.i}^2
  #
  # This single number drives both Ny/N and Ne/N. It captures how much
  # of the population survives from year to year in a stage-weighted way.
  # ------------------------------------------------------------------
  u_dot <- colSums(T_mat)
  u2_bar <- .compute_u2_bar(T_mat, D)

  # ------------------------------------------------------------------
  # Step 4: Compute Ny/N (Eq. 11) and Ne/N (Eq. 10)
  # ------------------------------------------------------------------
  rates <- .compute_Ne_clonal_core(u2_bar, L_use)

  # ------------------------------------------------------------------
  # Step 5: Compute minimum viable census size
  #
  # The Ne >= 5000 threshold (Lande 1995) is the standard conservation
  # benchmark for maintaining long-term gene diversity. Since Ne/N gives
  # us the ratio, the minimum N is simply: Ne_target / (Ne/N).
  # ------------------------------------------------------------------
  Min_N <- ceiling(Ne_target / rates$NeN)
  Ne_at_census <- if (!is.null(census_N)) {
    round(rates$NeN * census_N, 1)
  } else {
    NULL
  }

  # ------------------------------------------------------------------
  # Step 6: Attach stage names if T_mat has them
  # ------------------------------------------------------------------
  if (!is.null(colnames(T_mat))) {
    names(u_dot) <- colnames(T_mat)
  } else {
    names(u_dot) <- paste0("stage_", seq_len(s))
  }

  # ------------------------------------------------------------------
  # Step 7: Assemble and return results
  # ------------------------------------------------------------------
  result <- list(
    population = if (!is.null(population)) {
      as.character(population)
    } else {
      "unnamed"
    },
    model = "clonal_Y2000",
    NyN = rates$NyN,
    NeN = rates$NeN,
    L = L_use,
    L_source = L_source,
    u_dot = u_dot,
    u2_bar = u2_bar,
    Min_N = Min_N,
    Ne_target = Ne_target,
    census_N = census_N
  )

  class(result) <- c("Ne_clonal_Y2000", "NeStage")
  result
}


# -----------------------------------------------------------------------------
# SECTION 5: Print method
# -----------------------------------------------------------------------------
# A clean print method so that calling print(result) or just typing result
# at the console gives a readable summary rather than a raw list dump.

#' @export
print.Ne_clonal_Y2000 <- function(x, digits = 3, ...) {
  cat("\n")
  cat("=== NeStage: Clonal population Ne (Yonezawa 2000, Eq. 10-11) ===\n")
  cat(sprintf("  Population  : %s\n", x$population))
  cat(sprintf("  Model       : %s\n", x$model))
  cat("\n")
  cat("  --- Stage survival ---\n")
  for (nm in names(x$u_dot)) {
    cat(sprintf("    u_{.%s} = %.4f\n", nm, x$u_dot[nm]))
  }
  cat(sprintf(
    "    u2_bar (stage-weighted mean of squared survivals) = %.6f\n",
    x$u2_bar
  ))
  cat("\n")
  cat("  --- Generation time ---\n")
  cat(sprintf("    L = %.3f yr  [source: %s]\n", x$L, x$L_source))
  cat("\n")
  cat("  --- Effective size ratios ---\n")
  cat(sprintf("    Ny/N (annual, Eq. 11)          = %.3f\n", x$NyN))
  cat(sprintf("    Ne/N (generation-time, Eq. 10) = %.3f\n", x$NeN))
  cat("\n")
  cat("  --- Conservation threshold ---\n")
  cat(sprintf("    Ne target                      = %d\n", x$Ne_target))
  cat(sprintf("    Minimum census size N          = %d\n", x$Min_N))
  cat(sprintf(
    "    (Ne/N = %.3f => need N >= %d for Ne >= %d)\n",
    x$NeN,
    x$Min_N,
    x$Ne_target
  ))
  cat("\n")
  invisible(x)
}


# -----------------------------------------------------------------------------
# SECTION 6: Convenience wrapper -- run both observed and expected D
# -----------------------------------------------------------------------------
# In practice, users often want to compare results using the observed stage
# fractions from the field (D_obs) and the expected (equilibrium) fractions
# from the transition matrix eigenvector (D_exp). The paper reports both
# (Table 4, main values vs. parenthetical values). This wrapper runs both.

#' Run Ne_clonal_Y2000 with both observed and expected stage fractions
#'
#' @param T_mat    Survival/transition matrix (s x s).
#' @param F_vec    Fecundity vector (length s). Used for L and D_exp.
#' @param D_obs    Observed stage fractions from field counts.
#' @param D_exp    Expected (equilibrium) stage fractions. If NULL, derived
#'                 from the dominant eigenvector of A = T_mat + F_mat.
#' @param L        Generation time (years). If NULL, computed internally.
#' @param Ne_target  Numeric. Ne conservation threshold. Common choices:
#'               50 (Franklin 1980, avoid inbreeding depression),
#'               500 (long-term quantitative variation), or
#'               5000 (Lande 1995, long-term evolutionary potential).
#'               Default 50.
#' @param census_N  Numeric or NULL. Actual or expected census population
#'               size. If supplied, Ne_at_census = NeN * census_N is
#'               added to the output. Default NULL.
#' @param population  Character label for the population.
#'
#' @return A list with two Ne_clonal_Y2000 result objects:
#'   \describe{
#'     \item{observed}{Results using D_obs}
#'     \item{expected}{Results using D_exp}
#'   }
#'
#' @examples
#' T_Miz <- matrix(c(
#'   0.789, 0.121, 0.054,
#'   0.007, 0.621, 0.335,
#'   0.001, 0.258, 0.611
#' ), nrow = 3, byrow = TRUE)
#'
#' both <- Ne_clonal_Y2000_both(
#'   T_mat      = T_Miz,
#'   F_vec      = c(0.055, 1.328, 2.398),
#'   D_obs      = c(0.935, 0.038, 0.027),
#'   D_exp      = c(0.921, 0.046, 0.033),
#'   L          = 13.399,
#'   population = "Miz"
#' )
#' print(both$observed)   # main values in Table 4
#' print(both$expected)   # parenthetical values in Table 4
#'
#' @export
Ne_clonal_Y2000_both <- function(
  T_mat,
  F_vec,
  D_obs,
  D_exp = NULL,
  L = NULL,
  Ne_target = 50,
  census_N = NULL,
  population = NULL
) {
  # Derive D_exp from the stable stage distribution if not supplied
  if (is.null(D_exp)) {
    s <- nrow(T_mat)
    F_mat <- matrix(0, s, s)
    F_mat[1, ] <- F_vec
    A <- T_mat + F_mat
    ev <- eigen(A)
    w <- Re(ev$vectors[, which.max(Re(ev$values))])
    D_exp <- abs(w) / sum(abs(w))
    message(
      "D_exp not supplied -- derived from stable stage distribution of A."
    )
  }

  pop_label <- if (!is.null(population)) population else "unnamed"

  list(
    observed = Ne_clonal_Y2000(
      T_mat = T_mat,
      F_vec = F_vec,
      D = D_obs,
      L = L,
      Ne_target = Ne_target,
      census_N = census_N,
      population = paste0(pop_label, " (observed D)")
    ),
    expected = Ne_clonal_Y2000(
      T_mat = T_mat,
      F_vec = F_vec,
      D = D_exp,
      L = L,
      Ne_target = Ne_target,
      census_N = census_N,
      population = paste0(pop_label, " (expected D)")
    )
  )
}
