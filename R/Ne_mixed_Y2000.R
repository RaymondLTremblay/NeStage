# =============================================================================
# Ne_mixed_Y2000.R
# -----------------------------------------------------------------------------
# Variance effective population size for stage-structured populations with
# MIXED sexual and clonal reproduction, following Yonezawa et al. (2000),
# full general Equation 6.
#
# Reference:
#   Yonezawa K., Kinoshita E., Watano Y., and Zentoh H. (2000).
#   Formulation and estimation of the effective size of stage-structured
#   populations in Fritillaria camtschatcensis, a perennial herb with a
#   complex life history. Evolution 54(6): 2007-2013.
#   https://doi.org/10.1111/j.0014-3820.2000.tb01244.x
#
# When to use this function:
#   Use Ne_mixed_Y2000() when your study population reproduces through BOTH
#   sexual and clonal means, with each stage potentially having a different
#   mix of the two. Examples include:
#     - Many perennial herbs with both seed and vegetative reproduction
#     - Grasses that produce both seeds and tillers
#     - Corals that spawn sexually and fragment clonally
#     - Some ferns with both spore and rhizome reproduction
#
#   If your species reproduces ONLY clonally, use Ne_clonal_Y2000().
#   If your species reproduces ONLY sexually, use Ne_sexual_Y2000().
#
# Relationship to the other two functions:
#   Ne_mixed_Y2000() implements the full general Eq. 6.
#
#   Ne_sexual_Y2000() is a true special case of Ne_mixed_Y2000():
#     Setting d = rep(0, s) in Ne_mixed_Y2000() zeroes out Avr(Ad)
#     and reproduces Ne_sexual_Y2000() exactly.
#
#   Ne_clonal_Y2000() is NOT algebraically equivalent to Ne_mixed_Y2000()
#   with d = rep(1, s). This is because Ne_clonal_Y2000() implements
#   Equations 10-11, which are derived via a DIFFERENT simplification
#   path from Eq. 6 — one that uses the survival variance directly as
#   V = 2(1 - u2_bar), not the general V formula. Under Poisson defaults
#   (a=0, Vc/c_bar=1), A_i = 2*(1+0)*1 - S_i = 2 - 2 = 0, so V_term3
#   is always zero in Ne_mixed_Y2000() regardless of d. Use
#   Ne_clonal_Y2000() directly for purely clonal populations.
#
# The full general Eq. 6 of Yonezawa (2000):
#
#   V = 2(1+a)(u2_bar - u_bar^2) + (1 - u_bar) * [Avr(S) + Avr(Ad)]
#
#   where:
#     u_bar    = sum_i D_i * u_{.i}         (population mean annual survival)
#     u2_bar   = sum_i D_i * u_{.i}^2       (stage-weighted mean squared survival)
#     r_i      = F_i * D_i                  (newborn contribution of stage i)
#     d_i      = clonal reproduction rate of stage i; 1-d_i = sexual rate
#     S_i      = (1-a) + (1+a) * (Vk/k_bar)_i
#     A_i      = 2*(1+a) * (Vc/c_bar)_i - S_i
#     Avr(S)   = sum_i(r_i * S_i)  / sum_i(r_i)   (recruitment-weighted mean S_i)
#     Avr(Ad)  = sum_i(r_i * A_i * d_i) / sum_i(r_i)
#
#   Ne/N = 2 / (V * L)                      ... Eq. 6
#   Ny/N = 1 / (1 - u2_bar)                 ... Eq. 11 (optional, clonal component)
#
# The three terms in V and what they capture:
#
#   Term 1: 2(1+a)(u2_bar - u_bar^2)
#     Between-stage variance in SURVIVAL. Stages with very different
#     survival rates generate more genetic drift. Zero if all stages
#     have identical survival.
#
#   Term 2: (1 - u_bar) * Avr(S)
#     Variance in SEXUAL reproductive output. Driven by (Vk/k_bar)_i —
#     how unequal sexual reproductive success is within each stage.
#     Larger when few individuals dominate sexual reproduction.
#
#   Term 3: (1 - u_bar) * Avr(Ad)
#     Variance in CLONAL reproductive output. Driven by d_i (how much
#     of reproduction is clonal) and (Vc/c_bar)_i (how unequal clonal
#     output is). Zero when d_i = 0 for all stages (purely sexual).
#
# Author:  Raymond L. Tremblay
# Version: 1.0.0  (2026-03-05)
# =============================================================================


# -----------------------------------------------------------------------------
# SECTION 1: Input validation
# -----------------------------------------------------------------------------

.validate_T_mat_mx <- function(T_mat) {
  if (!is.matrix(T_mat) || !is.numeric(T_mat))
    stop("T_mat must be a numeric matrix.")
  if (nrow(T_mat) != ncol(T_mat))
    stop("T_mat must be square (nrow == ncol).")
  if (any(!is.finite(T_mat)))
    stop("T_mat contains non-finite values (NA, NaN, Inf).")
  if (any(T_mat < 0))
    stop("T_mat contains negative values. All transition rates must be >= 0.")
  col_sums <- colSums(T_mat)
  if (any(col_sums > 1 + 1e-8))
    stop(paste0(
      "Column sums of T_mat exceed 1.0. Total annual survival cannot exceed 1.\n",
      "  Column sums: ", paste(round(col_sums, 4), collapse = ", ")
    ))
  invisible(TRUE)
}

.validate_D_mx <- function(D, s) {
  if (!is.numeric(D) || length(D) != s)
    stop(paste0("D must be a numeric vector of length ", s,
                " (one entry per stage)."))
  if (any(!is.finite(D)) || any(D < 0))
    stop("D must contain finite non-negative values.")
  if (abs(sum(D) - 1) > 1e-6)
    stop(paste0("D must sum to 1. Current sum = ", round(sum(D), 6),
                ".\n  Hint: D <- D / sum(D) to normalise."))
  invisible(TRUE)
}

.validate_F_vec_mx <- function(F_vec, s) {
  if (!is.numeric(F_vec) || length(F_vec) != s)
    stop(paste0("F_vec must be a numeric vector of length ", s,
                " (one fecundity per stage)."))
  if (any(!is.finite(F_vec)) || any(F_vec < 0))
    stop("F_vec must contain finite non-negative values.")
  if (all(F_vec == 0))
    stop("F_vec is all zeros: no reproductive output in any stage.")
  invisible(TRUE)
}

.validate_L_mx <- function(L) {
  if (!is.numeric(L) || length(L) != 1)
    stop("L must be a single numeric value (generation time in years).")
  if (!is.finite(L) || L <= 0)
    stop(paste0("L must be a positive finite number. Received: ", L))
  invisible(TRUE)
}

.validate_d <- function(d, s) {
  # d is the per-stage clonal reproduction rate.
  # d_i = 1 means stage i reproduces entirely clonally.
  # d_i = 0 means stage i reproduces entirely sexually.
  # d_i = 0.7 means 70% of reproductive output from stage i is clonal.
  #
  # This must be supplied by the user — there is no sensible default
  # because the mix of sexual and clonal reproduction is species- and
  # stage-specific and cannot be assumed without field data.
  #
  # Example for a species where only adults reproduce clonally:
  #   d = c(0, 0, 0.8)  — stages 1 and 2 sexual only, stage 3 mostly clonal
  if (!is.numeric(d) || length(d) != s)
    stop(paste0(
      "d must be a numeric vector of length ", s,
      " (one clonal fraction per stage).\n",
      "  d_i = 0 : stage i reproduces entirely sexually\n",
      "  d_i = 1 : stage i reproduces entirely clonally\n",
      "  d_i = 0.7 : 70% of stage i reproduction is clonal\n",
      "  Example: d = c(0, 0, 0.8) for a 3-stage population where\n",
      "           only stage 3 has substantial clonal reproduction."
    ))
  if (any(!is.finite(d)) || any(d < 0) || any(d > 1))
    stop(paste0(
      "All values of d must be between 0 and 1.\n",
      "  d_i = 0 : purely sexual\n",
      "  d_i = 1 : purely clonal\n",
      "  Received: ", paste(round(d, 4), collapse = ", ")
    ))
  invisible(TRUE)
}

.validate_Vc_over_c <- function(Vc_over_c, s) {
  # Vc_over_c is the variance-to-mean ratio of CLONAL reproductive output
  # per stage. Parallel to Vk_over_k for sexual reproduction.
  #
  # Poisson default: Vc_over_c = 1 for all stages.
  #   Clonal output (e.g., bulblets, stolons, tillers) is equally
  #   distributed among individuals within each stage.
  #
  # Vc_over_c > 1: some individuals produce many more clonal propagules
  #   than others (e.g., large dominant plants producing many more
  #   rhizomes). This REDUCES Ne.
  #
  # Vc_over_c < 1: clonal output is more equal than random. Rare in
  #   nature but possible in some highly regulated clonal systems.
  if (!is.numeric(Vc_over_c) || length(Vc_over_c) != s)
    stop(paste0(
      "Vc_over_c must be a numeric vector of length ", s,
      " (one value per stage).\n",
      "  Default is rep(1, s) for Poisson clonal variance.\n",
      "  Values > 1: unequal clonal output (reduces Ne).\n",
      "  Values < 1: more equal than random clonal output (increases Ne)."
    ))
  if (any(!is.finite(Vc_over_c)) || any(Vc_over_c < 0))
    stop("Vc_over_c must contain finite non-negative values.")
  invisible(TRUE)
}

.validate_Vk_over_k_mx <- function(Vk_over_k, s) {
  if (!is.numeric(Vk_over_k) || length(Vk_over_k) != s)
    stop(paste0(
      "Vk_over_k must be a numeric vector of length ", s,
      " (one value per stage).\n",
      "  Default is rep(1, s) for Poisson sexual reproductive variance."
    ))
  if (any(!is.finite(Vk_over_k)) || any(Vk_over_k < 0))
    stop("Vk_over_k must contain finite non-negative values.")
  invisible(TRUE)
}

.validate_a_mx <- function(a) {
  if (!is.numeric(a) || length(a) != 1)
    stop("a must be a single numeric value (Hardy-Weinberg deviation).")
  if (!is.finite(a))
    stop("a must be a finite number.")
  if (a < -1 || a > 1)
    stop(paste0(
      "a must be between -1 and 1.\n",
      "  a = 0  : random mating (default)\n",
      "  a > 0  : inbreeding present\n",
      "  a < 0  : excess heterozygosity\n",
      "  Received: ", a
    ))
  invisible(TRUE)
}


# -----------------------------------------------------------------------------
# SECTION 2: Generation time L
# -----------------------------------------------------------------------------
# Same definition as Ne_clonal_Y2000 and Ne_sexual_Y2000:
# population-average mean age of reproduction, computed by iterating
# the survival matrix T to x_max = 500, weighted by the stable stage
# distribution w of A = T + F_mat.

.compute_L_mixed <- function(T_mat, F_vec, x_max = 500L) {
  s     <- nrow(T_mat)
  F_mat <- matrix(0, s, s)
  F_mat[1, ] <- F_vec

  A  <- T_mat + F_mat
  ev <- eigen(A)
  w  <- Re(ev$vectors[, which.max(Re(ev$values))])
  w  <- abs(w) / sum(abs(w))

  Tx  <- diag(s)
  num <- 0
  den <- 0

  for (x in seq_len(x_max)) {
    Tx <- T_mat %*% Tx

    l_bar_x <- 0
    m_bar_x <- 0
    for (i in seq_len(s)) {
      u_jxi   <- Tx[, i]
      l_bar_x <- l_bar_x + w[i] * sum(u_jxi)
      m_bar_x <- m_bar_x + w[i] * sum(F_vec * u_jxi)
    }

    num <- num + x * m_bar_x * l_bar_x
    den <- den +     m_bar_x * l_bar_x
  }

  if (den <= 0)
    stop(paste0(
      "Generation time L is undefined: no reproductive output detected over ",
      x_max, " years.\n",
      "  Check that F_vec has positive entries and T_mat allows survival."
    ))

  as.numeric(num / den)
}


# -----------------------------------------------------------------------------
# SECTION 3: Core computation — full Eq. 6
# -----------------------------------------------------------------------------

.compute_V_mixed <- function(T_mat, F_vec, D, d, Vk_over_k, Vc_over_c, a) {
  # Compute the full variance term V in Eq. 6.
  #
  # With both sexual and clonal reproduction:
  #   V = 2(1+a)(u2_bar - u_bar^2) + (1 - u_bar) * [Avr(S) + Avr(Ad)]
  #
  # Step by step:
  #   1. u_{.i}  = colSums(T_mat)
  #   2. u_bar   = sum_i D_i * u_{.i}
  #   3. u2_bar  = sum_i D_i * u_{.i}^2
  #   4. r_i     = F_i * D_i
  #   5. S_i     = (1-a) + (1+a) * (Vk/k_bar)_i
  #   6. A_i     = 2*(1+a) * (Vc/c_bar)_i - S_i
  #      NOTE: A_i subtracts S_i ONCE — this is exact per paper p. 2008
  #   7. Avr(S)  = sum(r_i * S_i)       / sum(r_i)
  #   8. Avr(Ad) = sum(r_i * A_i * d_i) / sum(r_i)
  #   9. V       = 2(1+a)(u2_bar - u_bar^2) + (1 - u_bar)(Avr_S + Avr_Ad)

  u_dot  <- colSums(T_mat)           # step 1
  u_bar  <- sum(D * u_dot)           # step 2
  u2_bar <- sum(D * u_dot^2)         # step 3
  r_i    <- F_vec * D                # step 4

  if (sum(r_i) <= 0)
    stop(paste0(
      "sum(r_i) = sum(F_vec * D) = 0.\n",
      "  At least one stage must have both F_vec > 0 and D > 0."
    ))

  # Step 5: S_i — captures sexual reproductive variance
  # Under defaults (a=0, Vk/k_bar=1): S_i = 2 for all stages
  S_i <- (1 - a) + (1 + a) * Vk_over_k

  # Step 6: A_i — captures EXCESS clonal over sexual reproductive variance
  # A_i = 2*(1+a)*(Vc/c_bar)_i - S_i   (definition per paper p.2008)
  #
  # IMPORTANT: Under Poisson defaults (a=0, Vc/c_bar=1, Vk/k_bar=1):
  #   S_i = (1-0) + (1+0)*1 = 2
  #   A_i = 2*(1+0)*1 - 2   = 0
  # A_i = 0 is CORRECT and INTENTIONAL. It means when clonal and sexual
  # output both follow Poisson distributions, changing d_i does NOT affect
  # Ne. V_term3 is non-zero only when Vc_over_c != Vk_over_k or a != 0.
  # Supply e.g. Vc_over_c = c(1,1,2) to see a non-zero V_term3.
  A_i <- 2 * (1 + a) * Vc_over_c - S_i

  # Step 7: Avr(S) — recruitment-weighted mean of S_i
  # Stages contributing more offspring have more influence on gene pool
  Avr_S  <- sum(r_i * S_i)       / sum(r_i)

  # Step 8: Avr(Ad) — recruitment-weighted mean of A_i * d_i
  # This term is zero when d_i = 0 for all stages (purely sexual case)
  # and grows as clonal reproduction increases
  Avr_Ad <- sum(r_i * A_i * d) / sum(r_i)

  # Step 9: Full V
  V_term1 <- 2 * (1 + a) * (u2_bar - u_bar^2)   # between-stage survival variance
  V_term2 <- (1 - u_bar) * Avr_S                 # sexual reproductive variance
  V_term3 <- (1 - u_bar) * Avr_Ad                # clonal reproductive variance
  V       <- V_term1 + V_term2 + V_term3

  if (!is.finite(V) || V <= 0)
    stop(paste0(
      "V = ", round(V, 8), " is not positive.\n",
      "  V must be > 0 for Ne to be defined.\n",
      "  Check T_mat survival rates, F_vec, d, and reproductive variance inputs."
    ))

  list(
    u_dot   = u_dot,
    u_bar   = u_bar,
    u2_bar  = u2_bar,
    r_i     = r_i,
    S_i     = S_i,
    A_i     = A_i,
    Avr_S   = Avr_S,
    Avr_Ad  = Avr_Ad,
    V       = V,
    V_term1 = V_term1,   # between-stage survival variance
    V_term2 = V_term2,   # sexual reproductive variance
    V_term3 = V_term3    # clonal reproductive variance
  )
}


# -----------------------------------------------------------------------------
# SECTION 4: Public function
# -----------------------------------------------------------------------------

#' Variance effective population size for mixed sexual/clonal stage-structured
#' populations
#'
#' Implements the full general Equation 6 of Yonezawa et al. (2000) for
#' populations with both sexual and clonal reproduction. Each stage can
#' have a different mix of the two reproductive modes via the per-stage
#' clonal fraction vector d.
#'
#' Note on relationship to the other two functions:
#' When d = rep(0, s), Ne_mixed_Y2000() reproduces Ne_sexual_Y2000() exactly,
#' because Avr(Ad) vanishes. However, Ne_mixed_Y2000() with d = rep(1, s) does
#' NOT reproduce Ne_clonal_Y2000(), because Ne_clonal_Y2000() implements the
#' simplified Equations 10-11 via a different algebraic path (V = 2(1-u2_bar))
#' that is not the same as the general Eq. 6 with d=1. For purely clonal
#' populations always use Ne_clonal_Y2000() directly.
#'
#' @param T_mat  Numeric matrix (s x s). Survival/transition matrix.
#'               Entry [j, i] = probability that a stage-i individual at
#'               year t is in stage j at year t+1. Column sums <= 1.
#'
#' @param F_vec  Numeric vector (length s). Total fecundity per stage per
#'               year: mean number of propagules (sexual + clonal combined)
#'               produced per individual. Used for L and r_i = F_i * D_i.
#'
#' @param D      Numeric vector (length s). Stage frequency distribution.
#'               Must sum to 1.
#'
#' @param d      Numeric vector (length s). Per-stage clonal reproduction
#'               fraction. MUST be supplied — there is no default.
#'               d_i = 0 : stage i reproduces entirely sexually.
#'               d_i = 1 : stage i reproduces entirely clonally.
#'               d_i = 0.7 : 70% of stage i reproduction is clonal.
#'               Example: d = c(0, 0, 0.8) for a 3-stage population where
#'               stages 1 and 2 are sexual only and stage 3 is mostly clonal.
#'
#' @param Vk_over_k  Numeric vector (length s). Variance-to-mean ratio of
#'               SEXUAL reproductive output per stage (Vk/k_bar)_i.
#'               Default rep(1, s) = Poisson: reproductive success is equally
#'               distributed among individuals within each stage.
#'               Values > 1 reduce Ne (unequal sexual reproductive success).
#'
#' @param Vc_over_c  Numeric vector (length s). Variance-to-mean ratio of
#'               CLONAL reproductive output per stage (Vc/c_bar)_i.
#'               Default rep(1, s) = Poisson: clonal output is equally
#'               distributed among individuals within each stage.
#'               Values > 1 reduce Ne (some individuals produce many more
#'               clonal propagules than others).
#'
#' @param a      Numeric scalar. Hardy-Weinberg deviation. Default 0 =
#'               random mating. Supply a positive value for inbreeding.
#'               Must be between -1 and 1.
#'
#' @param L      Numeric scalar (optional). Generation time in years.
#'               If NULL, computed internally via T^x iteration (x_max = 500).
#'
#' @param x_max  Integer. Maximum age for L iteration (default 500).
#'               Ignored if L is supplied directly.
#'
#' @param Ne_target  Numeric. Ne conservation threshold (default 5000,
#'               Lande 1995). Used to compute minimum viable census size.
#'
#' @param show_Ny  Logical. Whether to compute and return Ny/N (annual
#'               effective size, Eq. 11). Default FALSE. Note that Ny/N
#'               is defined here using the clonal component of the model
#'               and is most meaningful when clonal reproduction dominates.
#'
#' @param population  Character string. Optional population label.
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{population}{Character label}
#'     \item{model}{"mixed_Y2000"}
#'     \item{NeN}{Ne/N — generation-time effective size ratio (Eq. 6)}
#'     \item{NyN}{Ny/N — annual effective size ratio (Eq. 11), or NA if show_Ny = FALSE}
#'     \item{L}{Generation time used (years)}
#'     \item{L_source}{"user" or "computed"}
#'     \item{d}{Per-stage clonal fractions supplied}
#'     \item{u_dot}{Total annual survival per stage}
#'     \item{u_bar}{Population mean annual survival}
#'     \item{u2_bar}{Stage-weighted mean of squared survivals}
#'     \item{r_i}{Newborn contributions per stage}
#'     \item{S_i}{Per-stage S values}
#'     \item{A_i}{Per-stage A values}
#'     \item{Avr_S}{Recruitment-weighted mean of S_i}
#'     \item{Avr_Ad}{Recruitment-weighted mean of A_i * d_i}
#'     \item{V}{Full variance term}
#'     \item{V_term1}{Between-stage survival variance component}
#'     \item{V_term2}{Sexual reproductive variance component}
#'     \item{V_term3}{Clonal reproductive variance component}
#'     \item{Vk_over_k}{Sexual variance ratio used}
#'     \item{Vc_over_c}{Clonal variance ratio used}
#'     \item{a}{Hardy-Weinberg deviation used}
#'     \item{Min_N}{Minimum census N for Ne >= Ne_target}
#'     \item{Ne_target}{Ne threshold used}
#'   }
#'
#' @examples
#' # ---------------------------------------------------------
#' # Example 1: 3-stage perennial herb, stage 3 reproduces both
#' # sexually (seeds) and clonally (rhizomes). Stages 1 and 2
#' # reproduce sexually only.
#' # ---------------------------------------------------------
#' T_herb <- matrix(c(
#'   0.30, 0.05, 0.00,
#'   0.40, 0.65, 0.10,
#'   0.00, 0.20, 0.80
#' ), nrow = 3, byrow = TRUE)
#'
#' result <- Ne_mixed_Y2000(
#'   T_mat      = T_herb,
#'   F_vec      = c(0.0, 0.5, 3.0),
#'   D          = c(0.60, 0.25, 0.15),
#'   d          = c(0.0, 0.0, 0.7),   # stage 3 is 70% clonal
#'   population = "mixed herb"
#' )
#' print(result)
#'
#' # ---------------------------------------------------------
#' # Example 2: Show Ny/N for a predominantly clonal population
#' # ---------------------------------------------------------
#' result_Ny <- Ne_mixed_Y2000(
#'   T_mat      = T_herb,
#'   F_vec      = c(0.0, 0.5, 3.0),
#'   D          = c(0.60, 0.25, 0.15),
#'   d          = c(0.5, 0.8, 0.9),
#'   show_Ny    = TRUE,
#'   population = "predominantly clonal herb"
#' )
#' print(result_Ny)
#'

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
Ne_mixed_Y2000 <- function(T_mat,
                            F_vec,
                            D,
                            d,
                            Vk_over_k  = NULL,
                            Vc_over_c  = NULL,
                            a          = 0,
                            L          = NULL,
                            x_max      = 500L,
                            Ne_target  = 5000,
                            show_Ny    = FALSE,
                            population = NULL) {

  # ------------------------------------------------------------------
  # Step 1: Validate all inputs
  # ------------------------------------------------------------------
  s <- nrow(T_mat)
  .validate_T_mat_mx(T_mat)
  .validate_D_mx(D, s)
  .validate_F_vec_mx(F_vec, s)
  .validate_d(d, s)
  .validate_a_mx(a)
  if (!is.null(L)) .validate_L_mx(L)

  # Defaults: Poisson variance for both sexual and clonal output
  if (is.null(Vk_over_k)) {
    Vk_over_k <- rep(1, s)
  } else {
    .validate_Vk_over_k_mx(Vk_over_k, s)
  }

  if (is.null(Vc_over_c)) {
    Vc_over_c <- rep(1, s)
  } else {
    .validate_Vc_over_c(Vc_over_c, s)
  }

  if (!is.numeric(Ne_target) || Ne_target <= 0)
    stop("Ne_target must be a positive number (default 5000).")

  # ------------------------------------------------------------------
  # Step 2: Compute or retrieve generation time L
  # ------------------------------------------------------------------
  if (!is.null(L)) {
    L_use    <- as.numeric(L)
    L_source <- "user"
  } else {
    L_use    <- .compute_L_mixed(T_mat, F_vec, x_max = as.integer(x_max))
    L_source <- "computed"
  }

  # ------------------------------------------------------------------
  # Step 3: Compute full V and all intermediates
  # ------------------------------------------------------------------
  V_list <- .compute_V_mixed(T_mat, F_vec, D, d, Vk_over_k, Vc_over_c, a)

  # ------------------------------------------------------------------
  # Step 4: Compute Ne/N (Eq. 6)
  # ------------------------------------------------------------------
  NeN <- 2 / (V_list$V * L_use)

  if (!is.finite(NeN) || NeN <= 0)
    stop(paste0("Ne/N = ", NeN, " is not positive. Check V and L."))

  # ------------------------------------------------------------------
  # Step 5: Optionally compute Ny/N (Eq. 11)
  # Ny/N uses u2_bar, the same quantity as in Ne_clonal_Y2000.
  # It is most meaningful when clonal reproduction dominates (d_i near 1).
  # ------------------------------------------------------------------
  NyN <- if (isTRUE(show_Ny)) {
    denom <- 1 - V_list$u2_bar
    if (!is.finite(denom) || denom <= 0) {
      warning("Ny/N is undefined (1 - u2_bar <= 0). Returning NA.")
      NA_real_
    } else {
      1 / denom
    }
  } else {
    NA_real_
  }

  # ------------------------------------------------------------------
  # Step 6: Minimum viable census size
  # ------------------------------------------------------------------
  Min_N <- ceiling(Ne_target / NeN)

  # ------------------------------------------------------------------
  # Step 7: Attach stage names
  # ------------------------------------------------------------------
  stage_names <- if (!is.null(colnames(T_mat))) colnames(T_mat) else
    paste0("stage_", seq_len(s))

  names(V_list$u_dot) <- stage_names
  names(V_list$r_i)   <- stage_names
  names(V_list$S_i)   <- stage_names
  names(V_list$A_i)   <- stage_names
  names(Vk_over_k)    <- stage_names
  names(Vc_over_c)    <- stage_names
  names(d)            <- stage_names

  # ------------------------------------------------------------------
  # Step 8: Assemble and return
  # ------------------------------------------------------------------
  result <- list(
    population   = if (!is.null(population)) as.character(population) else "unnamed",
    model        = "mixed_Y2000",
    NeN          = NeN,
    NyN          = NyN,
    L            = L_use,
    L_source     = L_source,
    d            = d,
    u_dot        = V_list$u_dot,
    u_bar        = V_list$u_bar,
    u2_bar       = V_list$u2_bar,
    r_i          = V_list$r_i,
    S_i          = V_list$S_i,
    A_i          = V_list$A_i,
    Avr_S        = V_list$Avr_S,
    Avr_Ad       = V_list$Avr_Ad,
    V            = V_list$V,
    V_term1      = V_list$V_term1,
    V_term2      = V_list$V_term2,
    V_term3      = V_list$V_term3,
    Vk_over_k    = Vk_over_k,
    Vc_over_c    = Vc_over_c,
    a            = a,
    Min_N        = Min_N,
    Ne_target    = Ne_target,
    show_Ny      = show_Ny
  )

  class(result) <- c("Ne_mixed_Y2000", "NeStage")
  result
}


# -----------------------------------------------------------------------------
# SECTION 5: Print method
# -----------------------------------------------------------------------------

#' @export
print.Ne_mixed_Y2000 <- function(x, digits = 3, ...) {
  cat("\n")
  cat("=== NeStage: Mixed sexual/clonal Ne (Yonezawa 2000, Eq. 6) ===\n")
  cat(sprintf("  Population  : %s\n", x$population))
  cat(sprintf("  Model       : %s\n", x$model))
  cat("\n")
  cat("  --- Clonal fractions per stage (d_i) ---\n")
  for (nm in names(x$d))
    cat(sprintf("    d_%s = %.4f  (%.0f%% clonal, %.0f%% sexual)\n",
                nm, x$d[nm], 100 * x$d[nm], 100 * (1 - x$d[nm])))
  cat("\n")
  cat("  --- Stage survival ---\n")
  for (nm in names(x$u_dot))
    cat(sprintf("    u_{.%s} = %.4f\n", nm, x$u_dot[nm]))
  cat(sprintf("    u_bar  (population mean annual survival)          = %.6f\n",
              x$u_bar))
  cat(sprintf("    u2_bar (stage-weighted mean of squared survivals) = %.6f\n",
              x$u2_bar))
  cat("\n")
  cat("  --- Reproductive parameters ---\n")
  cat(sprintf("    a (Hardy-Weinberg deviation)                      = %.4f\n",
              x$a))
  for (nm in names(x$Vk_over_k))
    cat(sprintf("    Vk/k_bar (%s) = %.4f\n", nm, x$Vk_over_k[nm]))
  for (nm in names(x$Vc_over_c))
    cat(sprintf("    Vc/c_bar (%s) = %.4f\n", nm, x$Vc_over_c[nm]))
  cat(sprintf("    Avr(S)  (sexual reproductive variance term)       = %.6f\n",
              x$Avr_S))
  cat(sprintf("    Avr(Ad) (clonal reproductive variance term)       = %.6f\n",
              x$Avr_Ad))
  cat("\n")
  cat("  --- Variance decomposition ---\n")
  cat(sprintf("    V term 1 (between-stage survival variance)        = %.6f\n",
              x$V_term1))
  cat(sprintf("    V term 2 (sexual reproductive variance)           = %.6f\n",
              x$V_term2))
  cat(sprintf("    V term 3 (clonal reproductive variance)           = %.6f\n",
              x$V_term3))
  if (x$V_term3 == 0 && any(x$d > 0))
    cat("    Note: V term 3 = 0 because Vc/c_bar == Vk/k_bar and a == 0.\n",
        "          Supply Vc_over_c != 1 to model unequal clonal variance.\n")
  cat(sprintf("    V total                                           = %.6f\n",
              x$V))
  cat("\n")
  cat("  --- Generation time ---\n")
  cat(sprintf("    L = %.3f yr  [source: %s]\n", x$L, x$L_source))
  cat("\n")
  cat("  --- Effective size ratios ---\n")
  cat(sprintf("    Ne/N (generation-time, Eq. 6)  = %.3f\n", x$NeN))
  if (isTRUE(x$show_Ny) && is.finite(x$NyN)) {
    cat(sprintf("    Ny/N (annual, Eq. 11)          = %.3f\n", x$NyN))
  } else if (isTRUE(x$show_Ny) && !is.finite(x$NyN)) {
    cat("    Ny/N (annual, Eq. 11)          = undefined\n")
  } else {
    cat("    Ny/N                           = not requested (set show_Ny = TRUE)\n")
  }
  cat("\n")
  cat("  --- Conservation threshold ---\n")
  cat(sprintf("    Ne target              = %d\n",   x$Ne_target))
  cat(sprintf("    Minimum census size N  = %d\n",   x$Min_N))
  cat(sprintf("    (Ne/N = %.3f => need N >= %d for Ne >= %d)\n",
              x$NeN, x$Min_N, x$Ne_target))
  cat("\n")
  invisible(x)
}


# -----------------------------------------------------------------------------
# SECTION 6: Convenience wrapper — observed and expected D in one call
# -----------------------------------------------------------------------------

#' Run Ne_mixed_Y2000 with both observed and expected stage fractions
#'
#' @param T_mat      Survival/transition matrix (s x s).
#' @param F_vec      Fecundity vector (length s).
#' @param D_obs      Observed stage fractions from field counts.
#' @param d          Per-stage clonal fraction vector (length s).
#' @param D_exp      Expected (equilibrium) stage fractions. If NULL,
#'                   derived from stable stage distribution of A.
#' @param Vk_over_k  Sexual variance ratio (default Poisson).
#' @param Vc_over_c  Clonal variance ratio (default Poisson).
#' @param a          Hardy-Weinberg deviation (default 0).
#' @param L          Generation time. If NULL, computed internally.
#' @param Ne_target  Conservation Ne threshold (default 5000).
#' @param show_Ny    Logical. Compute Ny/N? (default FALSE).
#' @param population Character label for the population.
#'
#' @return A list with two Ne_mixed_Y2000 result objects:
#'   \describe{
#'     \item{observed}{Results using D_obs}
#'     \item{expected}{Results using D_exp}
#'   }
#'
#' @export
Ne_mixed_Y2000_both <- function(T_mat,
                                 F_vec,
                                 D_obs,
                                 d,
                                 D_exp      = NULL,
                                 Vk_over_k  = NULL,
                                 Vc_over_c  = NULL,
                                 a          = 0,
                                 L          = NULL,
                                 Ne_target  = 5000,
                                 show_Ny    = FALSE,
                                 population = NULL) {

  if (is.null(D_exp)) {
    s     <- nrow(T_mat)
    F_mat <- matrix(0, s, s); F_mat[1, ] <- F_vec
    A     <- T_mat + F_mat
    ev    <- eigen(A)
    w     <- Re(ev$vectors[, which.max(Re(ev$values))])
    D_exp <- abs(w) / sum(abs(w))
    message("D_exp not supplied — derived from stable stage distribution of A.")
  }

  pop_label <- if (!is.null(population)) population else "unnamed"

  list(
    observed = Ne_mixed_Y2000(
      T_mat      = T_mat, F_vec = F_vec, D = D_obs, d = d,
      Vk_over_k  = Vk_over_k, Vc_over_c = Vc_over_c,
      a = a, L = L, Ne_target = Ne_target, show_Ny = show_Ny,
      population = paste0(pop_label, " (observed D)")
    ),
    expected = Ne_mixed_Y2000(
      T_mat      = T_mat, F_vec = F_vec, D = D_exp, d = d,
      Vk_over_k  = Vk_over_k, Vc_over_c = Vc_over_c,
      a = a, L = L, Ne_target = Ne_target, show_Ny = show_Ny,
      population = paste0(pop_label, " (expected D)")
    )
  )
}
