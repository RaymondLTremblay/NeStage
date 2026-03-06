# =============================================================================
# Ne_sexual_Y2000.R
# -----------------------------------------------------------------------------
# Variance effective population size for 100% SEXUAL stage-structured
# populations, following Yonezawa et al. (2000), Equation 6 with d_i = 0.
#
# Reference:
#   Yonezawa K., Kinoshita E., Watano Y., and Zentoh H. (2000).
#   Formulation and estimation of the effective size of stage-structured
#   populations in Fritillaria camtschatcensis, a perennial herb with a
#   complex life history. Evolution 54(6): 2007-2013.
#   https://doi.org/10.1111/j.0014-3820.2000.tb01244.x
#
# When to use this function:
#   Use Ne_sexual_Y2000() when your study population reproduces EXCLUSIVELY
#   through sexual means (seeds, spores, eggs) and clonal reproduction
#   contributes nothing meaningful to recruitment. Examples include most
#   vertebrates, many annual plants, and outcrossing trees.
#
#   If your species has ANY clonal reproduction contributing to recruitment,
#   use Ne_mixed_Y2000() instead.
#
# How this differs from Ne_clonal_Y2000:
#   In the clonal function, all reproduction is vegetative (d_i = 1), and
#   the formula simplifies to Equations 10-11. Here, all reproduction is
#   sexual (d_i = 0 for all stages), which means:
#
#     1. The Avr(Ad) term in Eq. 6 vanishes entirely because d_i = 0.
#     2. The formula is driven by Avr(S), which depends on the variance of
#        SEXUAL reproductive output (Vk/k_bar)_i per stage.
#     3. The annual effective size Ny is NOT computed — Eq. 11 was derived
#        specifically for the clonal case and does not apply here.
#
# Model assumptions:
#   1. The population is at demographic equilibrium: stage fractions D_i and
#      transition rates u_ji are constant over time.
#   2. All reproduction is sexual (d_i = 0 for all stages i).
#   3. The population mates randomly (no selfing, no inbreeding): a = 0.
#      Supply a non-zero value via the 'a' argument if inbreeding is known.
#   4. Sexual reproductive output follows a Poisson distribution by default:
#      (Vk / k_bar)_i = 1 for all stages i.
#      This means reproductive success is equally distributed among
#      individuals within a stage — no individual is a consistently
#      better or worse reproducer. See details on Vk_over_k below.
#   5. Census size N is approximately constant over the generation time L.
#
# The general Eq. 6 of Yonezawa (2000) with d_i = 0:
#
#   V  = 2(1 + a)(u2_bar - u_bar^2) + (1 - u_bar) * Avr(S)
#
#   where:
#     u_bar   = sum_i D_i * u_{.i}          (population mean annual survival)
#     u2_bar  = sum_i D_i * u_{.i}^2        (stage-weighted mean of squared survivals)
#     S_i     = (1 - a) + (1 + a) * (Vk/k_bar)_i
#     Avr(S)  = sum_i(r_i * S_i) / sum_i(r_i)   (recruitment-weighted mean of S_i)
#     r_i     = F_i * D_i                   (newborn contribution of stage i)
#
#   Ne/N = 2 / (V * L)                      ... Eq. 6  (generation-time effective size)
#
# Under Poisson defaults (a = 0, Vk/k_bar = 1):
#   S_i = 2 for all stages
#   Avr(S) = 2
#   V = 2(u2_bar - u_bar^2) + 2(1 - u_bar)
#
# IMPORTANT — Two sources of genetic drift in this model:
#   The V term has two additive components, each capturing a different
#   biological source of allele-frequency change:
#
#   Component 1: 2(1 + a)(u2_bar - u_bar^2)
#     This captures BETWEEN-STAGE variance in survival. Stages with
#     very different survival rates (e.g., seedlings dying at high rates
#     while adults persist) generate more drift than stages with uniform
#     survival. If all stages had identical survival, this term = 0.
#
#   Component 2: (1 - u_bar) * Avr(S)
#     This captures REPRODUCTIVE VARIANCE among individuals. The factor
#     (1 - u_bar) is the annual recruitment rate — the fraction of the
#     population replaced by new individuals each year. Avr(S) scales
#     this by how unequal reproductive success is. High variance in who
#     reproduces (large Vk/k_bar) means fewer individuals effectively
#     contribute genes, reducing Ne.
#
# Author:  Raymond L. Tremblay
# Version: 1.0.0  (2026-03-05)
# =============================================================================

# -----------------------------------------------------------------------------
# SECTION 1: Input validation
# -----------------------------------------------------------------------------
# These validators are parallel to those in Ne_clonal_Y2000.R.
# They are redefined here so this script is fully self-contained and does
# not depend on Ne_clonal_Y2000.R being sourced first.

.validate_T_mat_sx <- function(T_mat) {
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

.validate_D_sx <- function(D, s) {
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

.validate_F_vec_sx <- function(F_vec, s) {
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
    stop(paste0(
      "F_vec is all zeros: no reproductive output in any stage.\n",
      "  For a sexual population F_vec represents seeds/offspring per ",
      "individual per stage per year.\n",
      "  At least one stage must have F_vec > 0."
    ))
  }
  invisible(TRUE)
}

.validate_L_sx <- function(L) {
  if (!is.numeric(L) || length(L) != 1) {
    stop("L must be a single numeric value (generation time in years).")
  }
  if (!is.finite(L) || L <= 0) {
    stop(paste0("L must be a positive finite number. Received: ", L))
  }
  invisible(TRUE)
}

.validate_Vk_over_k <- function(Vk_over_k, s) {
  # Vk_over_k is the variance-to-mean ratio of sexual reproductive output
  # per stage. It controls how unequal reproductive success is among
  # individuals within each stage.
  #
  # Poisson default: Vk_over_k = 1 for all stages.
  #   Means reproductive success is randomly and equally distributed —
  #   no individual is a consistently better or worse reproducer.
  #
  # Vk_over_k > 1: some individuals reproduce much more than others
  #   (e.g., dominant plants, strong pollinator preference, resource
  #   limitation). This REDUCES Ne because fewer individuals effectively
  #   contribute to the gene pool.
  #
  # Vk_over_k < 1: reproductive success is more equal than random
  #   (rare in nature; possible in managed populations). This INCREASES Ne.
  if (!is.numeric(Vk_over_k) || length(Vk_over_k) != s) {
    stop(paste0(
      "Vk_over_k must be a numeric vector of length ",
      s,
      " (one value per stage).\n",
      "  Default is rep(1, s) for Poisson reproductive variance.\n",
      "  Values > 1 mean unequal reproductive success (reduces Ne).\n",
      "  Values < 1 mean more equal than random (increases Ne)."
    ))
  }
  if (any(!is.finite(Vk_over_k)) || any(Vk_over_k < 0)) {
    stop("Vk_over_k must contain finite non-negative values.")
  }
  invisible(TRUE)
}

.validate_a_sx <- function(a) {
  # a is the deviation from Hardy-Weinberg proportions.
  # a = 0  : random mating (default, recommended for most outcrossing species)
  # a > 0  : inbreeding (consanguineous mating), reduces Ne
  # a = -1 : maximum heterozygosity excess (theoretical lower bound)
  # a < -1 : not biologically meaningful
  if (!is.numeric(a) || length(a) != 1) {
    stop("a must be a single numeric value (Hardy-Weinberg deviation).")
  }
  if (!is.finite(a)) {
    stop("a must be a finite number.")
  }
  if (a < -1 || a > 1) {
    stop(paste0(
      "a must be between -1 and 1.\n",
      "  a = 0  : random mating (default)\n",
      "  a > 0  : inbreeding present\n",
      "  a < 0  : excess heterozygosity\n",
      "  Received: ",
      a
    ))
  }
  invisible(TRUE)
}


# -----------------------------------------------------------------------------
# SECTION 2: Generation time L
# -----------------------------------------------------------------------------
# Same definition as in Ne_clonal_Y2000.R: the population-average mean age
# of reproduction, computed by iterating the survival matrix T to x_max = 500,
# weighted by the stable stage distribution w of A = T + F_mat.
#
# For sexual populations, F_vec represents seeds or offspring produced per
# individual per stage per year. The computation is identical to the clonal
# case because L depends only on the demographic structure (T, F, w),
# not on whether reproduction is sexual or clonal.

.compute_L_sexual <- function(T_mat, F_vec, x_max = 500L) {
  # Compute generation time L using the Yonezawa (2000) definition.
  #
  # Args:
  #   T_mat : s x s survival/transition matrix
  #   F_vec : length-s fecundity vector (offspring per individual per stage per year)
  #   x_max : maximum age to iterate (default 500)
  #
  # Returns: L as a single numeric value (years)

  s <- nrow(T_mat)
  F_mat <- matrix(0, s, s)
  F_mat[1, ] <- F_vec # offspring enter stage 1 (juvenile stage)

  # Stable stage distribution: dominant right eigenvector of A = T + F_mat
  A <- T_mat + F_mat
  ev <- eigen(A)
  w <- Re(ev$vectors[, which.max(Re(ev$values))])
  w <- abs(w) / sum(abs(w)) # normalise to sum = 1

  # Iterate T^x, accumulating the generation time numerator and denominator
  Tx <- diag(s) # T^0 = identity
  num <- 0 # sum(x * m_bar_x * l_bar_x)
  den <- 0 # sum(m_bar_x * l_bar_x)

  for (x in seq_len(x_max)) {
    Tx <- T_mat %*% Tx # T^x = T * T^(x-1)

    l_bar_x <- 0
    m_bar_x <- 0
    for (i in seq_len(s)) {
      u_jxi <- Tx[, i]
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
# SECTION 3: Core computation — Eq. 6 with d_i = 0
# -----------------------------------------------------------------------------

.compute_V_sexual <- function(T_mat, F_vec, D, Vk_over_k, a) {
  # Compute the variance term V in Eq. 6 for a purely sexual population.
  #
  # With d_i = 0, Eq. 6 reduces to:
  #   V = 2(1 + a)(u2_bar - u_bar^2) + (1 - u_bar) * Avr(S)
  #
  # Step-by-step:
  #   1. u_{.i}  = colSums(T_mat)           total annual survival per stage
  #   2. u_bar   = sum_i D_i * u_{.i}       population mean annual survival
  #   3. u2_bar  = sum_i D_i * u_{.i}^2    stage-weighted mean squared survival
  #   4. r_i     = F_i * D_i               newborn contribution per stage
  #   5. S_i     = (1-a) + (1+a)*(Vk/k_bar)_i
  #   6. Avr(S)  = sum(r_i * S_i) / sum(r_i)
  #   7. V       = 2(1+a)(u2_bar - u_bar^2) + (1 - u_bar) * Avr(S)

  s <- length(D)
  u_dot <- colSums(T_mat) # step 1: u_{.i} for each stage
  u_bar <- sum(D * u_dot) # step 2: population mean survival
  u2_bar <- sum(D * u_dot^2) # step 3: mean of squared survivals

  r_i <- F_vec * D # step 4: newborn contributions

  if (sum(r_i) <= 0) {
    stop(paste0(
      "sum(r_i) = sum(F_vec * D) = 0.\n",
      "  The sexual model requires at least one stage with both F_vec > 0\n",
      "  and D > 0 so that Avr(S) can be computed."
    ))
  }

  # Step 5: S_i per stage
  # S_i = (1 - a) + (1 + a) * (Vk/k_bar)_i
  # Under defaults (a=0, Vk/k_bar=1): S_i = 1 + 1 = 2 for all stages.
  # Under inbreeding (a > 0): S_i shifts upward, increasing V and reducing Ne.
  # Under high reproductive variance (Vk/k_bar > 1): S_i increases,
  #   meaning fewer individuals effectively contribute genes -> lower Ne.
  S_i <- (1 - a) + (1 + a) * Vk_over_k

  # Step 6: Avr(S) — recruitment-weighted average of S_i
  # Weighted by r_i because stages that contribute more offspring to the
  # population have more influence on the genetic composition of the next
  # generation.
  Avr_S <- sum(r_i * S_i) / sum(r_i)

  # Step 7: V — the full variance term
  # Component 1: between-stage survival variance
  # Component 2: reproductive variance weighted by annual turnover
  V <- 2 * (1 + a) * (u2_bar - u_bar^2) + (1 - u_bar) * Avr_S

  if (!is.finite(V) || V <= 0) {
    stop(paste0(
      "V = ",
      round(V, 8),
      " is not positive.\n",
      "  V must be > 0 for Ne to be defined.\n",
      "  Check that T_mat has realistic survival rates and F_vec > 0\n",
      "  in at least one stage with D > 0."
    ))
  }

  # Return all intermediates so the user can inspect the computation
  list(
    u_dot = u_dot,
    u_bar = u_bar,
    u2_bar = u2_bar,
    r_i = r_i,
    S_i = S_i,
    Avr_S = Avr_S,
    V = V,
    # Decomposed components for transparency
    V_component1 = 2 * (1 + a) * (u2_bar - u_bar^2),
    V_component2 = (1 - u_bar) * Avr_S
  )
}

.compute_Ne_sexual_core <- function(V_list, L) {
  # Compute Ne/N from the variance term V and generation time L.
  #
  # Eq. 6:  Ne/N = 2 / (V * L)
  #
  # Note: unlike the clonal case, there is no Ny/N here.
  # The annual effective size Ny (Eq. 11) was derived specifically for
  # the clonal-dominant Poisson case and is not defined for sexual populations.

  NeN <- 2 / (V_list$V * L)

  if (!is.finite(NeN) || NeN <= 0) {
    stop(paste0(
      "Ne/N = ",
      NeN,
      " is not a positive finite number.\n",
      "  V = ",
      round(V_list$V, 8),
      ", L = ",
      round(L, 4),
      "\n",
      "  Check inputs: V and L must both be positive."
    ))
  }

  NeN
}


# -----------------------------------------------------------------------------
# SECTION 4: Public function
# -----------------------------------------------------------------------------

#' Variance effective population size for 100% sexual stage-structured populations
#'
#' Implements Equation 6 of Yonezawa et al. (2000) with d_i = 0 (no clonal
#' reproduction) for populations that reproduce exclusively through sexual
#' means. The annual effective size Ny is not computed as it applies only
#' to clonal populations (Eq. 11).
#'
#' @param T_mat  Numeric matrix (s x s). The survival/transition matrix.
#'               Entry [j, i] is the probability that an individual in stage i
#'               at year t is found in stage j at year t+1.
#'               Column sums must be <= 1 (total annual survival per stage).
#'
#' @param F_vec  Numeric vector of length s. Fecundity of each stage: the
#'               mean number of offspring (seeds, juveniles) produced per
#'               individual per year through sexual reproduction.
#'               Used to compute L and r_i = F_i * D_i.
#'               Set F_vec = NULL only if supplying L directly AND the
#'               function does not need r_i (not recommended).
#'
#' @param D      Numeric vector of length s. Stage frequency distribution:
#'               the proportion of the population in each stage. Must sum to 1.
#'
#' @param Vk_over_k  Numeric vector of length s. The variance-to-mean ratio
#'               of sexual reproductive output per stage: (Vk / k_bar)_i.
#'               Default is rep(1, s) — Poisson reproductive variance, meaning
#'               reproductive success is equally distributed among individuals
#'               within each stage. Values > 1 indicate that some individuals
#'               reproduce much more than others (common when pollinators show
#'               strong preferences, or when resources are limiting), which
#'               reduces Ne. Values < 1 indicate more equal than random
#'               reproductive success, which increases Ne.
#'
#' @param a      Numeric scalar. Deviation from Hardy-Weinberg proportions.
#'               Default 0 = random mating (appropriate for most outcrossing
#'               species with no selfing). Supply a positive value if the
#'               population has known inbreeding (e.g., from FIS estimates).
#'               Must be between -1 and 1.
#'
#' @param L      Numeric scalar (optional). Generation time in years.
#'               If NULL (default), L is computed internally from T_mat and
#'               F_vec using the Yonezawa (2000) definition (x_max = 500 yr).
#'               Supply L directly to match a specific published value.
#'
#' @param x_max  Integer. Maximum age for the L iteration (default 500).
#'               Ignored if L is supplied directly.
#'
#' @param Ne_target  Numeric. Ne conservation threshold for computing minimum
#'               viable census size (default 5000, following Lande 1995).
#'
#' @param population  Character string. Optional label for the population.
#'
#' @return A named list with the following elements:
#'   \describe{
#'     \item{population}{Character label}
#'     \item{model}{"sexual_Y2000"}
#'     \item{NeN}{Ne/N — generation-time effective size ratio (Eq. 6)}
#'     \item{L}{Generation time used (years)}
#'     \item{L_source}{"user" or "computed"}
#'     \item{u_dot}{Total annual survival per stage (colSums of T_mat)}
#'     \item{u_bar}{Population mean annual survival}
#'     \item{u2_bar}{Stage-weighted mean of squared survivals}
#'     \item{r_i}{Newborn contributions per stage (F_vec * D)}
#'     \item{S_i}{Per-stage S values used in Avr(S)}
#'     \item{Avr_S}{Recruitment-weighted mean of S_i}
#'     \item{V}{Full variance term from Eq. 6}
#'     \item{V_component1}{Between-stage survival variance component}
#'     \item{V_component2}{Reproductive variance component}
#'     \item{Vk_over_k}{Variance-to-mean ratio of sexual output used}
#'     \item{a}{Hardy-Weinberg deviation used}
#'     \item{Min_N}{Minimum census N for Ne >= Ne_target}
#'     \item{Ne_target}{Ne threshold used}
#'   }
#'
#' @examples
#' # ---------------------------------------------------------
#' # Example 1: Simple 3-stage plant, Poisson defaults
#' # A hypothetical outcrossing perennial herb with three stages:
#' #   Stage 1: seedling/juvenile (high mortality, no reproduction)
#' #   Stage 2: vegetative adult  (moderate survival, low reproduction)
#' #   Stage 3: reproductive adult (high survival, high reproduction)
#' # ---------------------------------------------------------
#' T_plant <- matrix(c(
#'   0.30, 0.05, 0.00,
#'   0.40, 0.65, 0.10,
#'   0.00, 0.20, 0.80
#' ), nrow = 3, byrow = TRUE)
#'
#' F_plant <- c(0.0, 0.5, 3.0)   # seeds per individual per stage per year
#' D_plant <- c(0.60, 0.25, 0.15) # observed stage fractions
#'
#' result <- Ne_sexual_Y2000(
#'   T_mat      = T_plant,
#'   F_vec      = F_plant,
#'   D          = D_plant,
#'   population = "hypothetical herb"
#' )
#' print(result)
#'
#' # ---------------------------------------------------------
#' # Example 2: Same population, high reproductive variance in stage 3
#' # Vk/k_bar = 3 for reproductive adults — pollinator-limited, so only
#' # a few adults contribute most of the seeds in any given year.
#' # This should reduce Ne relative to Example 1.
#' # ---------------------------------------------------------
#' result_highvar <- Ne_sexual_Y2000(
#'   T_mat       = T_plant,
#'   F_vec       = F_plant,
#'   D           = D_plant,
#'   Vk_over_k   = c(1, 1, 3),
#'   population  = "hypothetical herb (high repro variance stage 3)"
#' )
#' print(result_highvar)
#' # Ne/N should be lower than Example 1.
#'
#' # ---------------------------------------------------------
#' # Example 3: Supply L directly from a published source
#' # ---------------------------------------------------------
#' result_Luser <- Ne_sexual_Y2000(
#'   T_mat      = T_plant,
#'   F_vec      = F_plant,
#'   D          = D_plant,
#'   L          = 8.5,
#'   population = "hypothetical herb (published L)"
#' )
#' print(result_Luser)
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
Ne_sexual_Y2000 <- function(
  T_mat,
  F_vec,
  D,
  Vk_over_k = NULL,
  a = 0,
  L = NULL,
  x_max = 500L,
  Ne_target = 5000,
  population = NULL
) {
  # ------------------------------------------------------------------
  # Step 1: Validate all inputs
  # ------------------------------------------------------------------
  s <- nrow(T_mat)
  .validate_T_mat_sx(T_mat)
  .validate_D_sx(D, s)
  .validate_F_vec_sx(F_vec, s)
  .validate_a_sx(a)
  if (!is.null(L)) {
    .validate_L_sx(L)
  }

  # Default Vk_over_k = 1 for all stages (Poisson reproductive variance)
  if (is.null(Vk_over_k)) {
    Vk_over_k <- rep(1, s)
  } else {
    .validate_Vk_over_k(Vk_over_k, s)
  }

  if (!is.numeric(Ne_target) || Ne_target <= 0) {
    stop("Ne_target must be a positive number (default 5000).")
  }

  # ------------------------------------------------------------------
  # Step 2: Compute or retrieve generation time L
  # ------------------------------------------------------------------
  if (!is.null(L)) {
    L_use <- as.numeric(L)
    L_source <- "user"
  } else {
    L_use <- .compute_L_sexual(T_mat, F_vec, x_max = as.integer(x_max))
    L_source <- "computed"
  }

  # ------------------------------------------------------------------
  # Step 3: Compute V and all intermediate quantities
  # ------------------------------------------------------------------
  V_list <- .compute_V_sexual(T_mat, F_vec, D, Vk_over_k, a)

  # ------------------------------------------------------------------
  # Step 4: Compute Ne/N (Eq. 6 with d_i = 0)
  # ------------------------------------------------------------------
  NeN <- .compute_Ne_sexual_core(V_list, L_use)

  # ------------------------------------------------------------------
  # Step 5: Minimum viable census size (Lande 1995 criterion)
  # ------------------------------------------------------------------
  Min_N <- ceiling(Ne_target / NeN)

  # ------------------------------------------------------------------
  # Step 6: Attach stage names if available
  # ------------------------------------------------------------------
  stage_names <- if (!is.null(colnames(T_mat))) {
    colnames(T_mat)
  } else {
    paste0("stage_", seq_len(s))
  }

  names(V_list$u_dot) <- stage_names
  names(V_list$r_i) <- stage_names
  names(V_list$S_i) <- stage_names
  names(Vk_over_k) <- stage_names

  # ------------------------------------------------------------------
  # Step 7: Assemble and return results
  # ------------------------------------------------------------------
  result <- list(
    population = if (!is.null(population)) {
      as.character(population)
    } else {
      "unnamed"
    },
    model = "sexual_Y2000",
    NeN = NeN,
    L = L_use,
    L_source = L_source,
    u_dot = V_list$u_dot,
    u_bar = V_list$u_bar,
    u2_bar = V_list$u2_bar,
    r_i = V_list$r_i,
    S_i = V_list$S_i,
    Avr_S = V_list$Avr_S,
    V = V_list$V,
    V_component1 = V_list$V_component1,
    V_component2 = V_list$V_component2,
    Vk_over_k = Vk_over_k,
    a = a,
    Min_N = Min_N,
    Ne_target = Ne_target
  )

  class(result) <- c("Ne_sexual_Y2000", "NeStage")
  result
}


# -----------------------------------------------------------------------------
# SECTION 5: Print method
# -----------------------------------------------------------------------------

#' @export
print.Ne_sexual_Y2000 <- function(x, digits = 3, ...) {
  cat("\n")
  cat("=== NeStage: Sexual population Ne (Yonezawa 2000, Eq. 6, d=0) ===\n")
  cat(sprintf("  Population  : %s\n", x$population))
  cat(sprintf("  Model       : %s\n", x$model))
  cat("\n")
  cat("  --- Stage survival ---\n")
  for (nm in names(x$u_dot)) {
    cat(sprintf("    u_{.%s} = %.4f\n", nm, x$u_dot[nm]))
  }
  cat(sprintf(
    "    u_bar  (population mean annual survival)          = %.6f\n",
    x$u_bar
  ))
  cat(sprintf(
    "    u2_bar (stage-weighted mean of squared survivals) = %.6f\n",
    x$u2_bar
  ))
  cat("\n")
  cat("  --- Reproductive parameters ---\n")
  cat(sprintf(
    "    a (Hardy-Weinberg deviation)                      = %.4f\n",
    x$a
  ))
  for (nm in names(x$Vk_over_k)) {
    cat(sprintf("    Vk/k_bar (%s) = %.4f\n", nm, x$Vk_over_k[nm]))
  }
  cat(sprintf(
    "    Avr(S) (recruitment-weighted mean S_i)            = %.6f\n",
    x$Avr_S
  ))
  cat("\n")
  cat("  --- Variance decomposition ---\n")
  cat(sprintf(
    "    V component 1 (between-stage survival variance)   = %.6f\n",
    x$V_component1
  ))
  cat(sprintf(
    "    V component 2 (reproductive variance)             = %.6f\n",
    x$V_component2
  ))
  cat(sprintf(
    "    V total                                           = %.6f\n",
    x$V
  ))
  cat("\n")
  cat("  --- Generation time ---\n")
  cat(sprintf("    L = %.3f yr  [source: %s]\n", x$L, x$L_source))
  cat("\n")
  cat("  --- Effective size ratio ---\n")
  cat(sprintf("    Ne/N (generation-time, Eq. 6) = %.3f\n", x$NeN))
  cat("    Note: Ny/N is not defined for purely sexual populations.\n")
  cat("\n")
  cat("  --- Conservation threshold ---\n")
  cat(sprintf("    Ne target              = %d\n", x$Ne_target))
  cat(sprintf("    Minimum census size N  = %d\n", x$Min_N))
  cat(sprintf(
    "    (Ne/N = %.3f => need N >= %d for Ne >= %d)\n",
    x$NeN,
    x$Min_N,
    x$Ne_target
  ))
  cat("\n")
  invisible(x)
}

# NOTE: Sensitivity analysis functions (e.g., varying Vk_over_k, D, L, u_ji)
# have been moved to Ne_sensitivity.R, which will contain sensitivity tools
# for all three model functions (clonal, sexual, mixed) in one place.
