# =============================================================================
# Ne_sensitivity.R
# -----------------------------------------------------------------------------
# Sensitivity analysis functions for the three NeStage effective population
# size models (Ne_clonal_Y2000, Ne_sexual_Y2000, Ne_mixed_Y2000).
#
# Each function sweeps one parameter across a user-defined range, holding
# all other inputs constant, and returns:
#   (1) a data frame of results ready for custom analysis or plotting
#   (2) a ggplot2 figure showing Ne/N vs. the swept parameter
#   (3) the matched call for reproducibility
#
# Functions in this script:
#   Ne_sensitivity_Vk()  -- sweep (Vk/k_bar)_i for one stage (sexual, mixed)
#   Ne_sensitivity_Vc()  -- sweep (Vc/c_bar)_i for one stage (mixed only)
#   Ne_sensitivity_d()   -- sweep d_i for one stage from 0 to 1 (mixed only)
#   Ne_sensitivity_L()   -- sweep generation time L (all three models)
#
# Dependencies:
#   ggplot2 (>= 3.4.0)
#   Ne_clonal_Y2000.R, Ne_sexual_Y2000.R, Ne_mixed_Y2000.R must be sourced
#
# Author:  Raymond L. Tremblay
# Version: 1.0.0  (2026-03-05)
# =============================================================================

# Check that ggplot2 is available
if (!requireNamespace("ggplot2", quietly = TRUE))
  stop("Package 'ggplot2' is required for Ne_sensitivity functions.\n",
       "  Install it with: install.packages('ggplot2')")


# -----------------------------------------------------------------------------
# Internal helper: resolve stage index from index or name
# -----------------------------------------------------------------------------
.resolve_stage <- function(stage_index, stage_name, s, T_mat) {
  # Accepts either stage_index (integer) or stage_name (character).
  # stage_index takes priority if both are supplied.
  # Returns a validated integer stage index.

  col_names <- if (!is.null(colnames(T_mat))) colnames(T_mat) else
    paste0("stage_", seq_len(s))

  if (!is.null(stage_index)) {
    if (!is.numeric(stage_index) || length(stage_index) != 1 ||
        stage_index < 1 || stage_index > s)
      stop(paste0("stage_index must be a single integer between 1 and ", s, "."))
    return(as.integer(stage_index))
  }

  if (!is.null(stage_name)) {
    idx <- match(stage_name, col_names)
    if (is.na(idx))
      stop(paste0(
        "stage_name '", stage_name, "' not found in column names of T_mat.\n",
        "  Available names: ", paste(col_names, collapse = ", ")
      ))
    return(idx)
  }

  stop("Supply either stage_index (integer) or stage_name (character).")
}


# -----------------------------------------------------------------------------
# Internal helper: shared ggplot2 theme for all sensitivity figures
# -----------------------------------------------------------------------------
.ne_sensitivity_theme <- function() {
  ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", size = 13),
      plot.subtitle = ggplot2::element_text(size = 10, color = "grey40"),
      plot.caption  = ggplot2::element_text(size = 8,  color = "grey50",
                                            hjust = 0),
      axis.title    = ggplot2::element_text(size = 11),
      axis.text     = ggplot2::element_text(size = 10),
      panel.grid.major.y = ggplot2::element_line(color = "grey90",
                                                  linewidth = 0.4),
      legend.position = "none"
    )
}


# =============================================================================
# FUNCTION 1: Ne_sensitivity_Vk
# =============================================================================
# Sweeps (Vk/k_bar)_i — the variance-to-mean ratio of SEXUAL reproductive
# output — for a single focal stage across a user-defined range.
#
# Applies to: Ne_sexual_Y2000 and Ne_mixed_Y2000.
#
# Biological question:
#   How sensitive is Ne to the assumed level of sexual reproductive inequality
#   in the focal stage? Under Poisson (Vk/k_bar = 1), all individuals
#   contribute equally. Larger values mean some individuals dominate sexual
#   reproduction, reducing Ne. This is relevant when:
#     - Pollinator preference strongly favors certain individuals
#     - Resources limit sexual reproduction to a few dominant plants
#     - The user does not have empirical data and wants to bracket uncertainty

#' Sensitivity of Ne/N to sexual reproductive variance in one stage
#'
#' Sweeps (Vk/k_bar)_i across a range for a single focal stage, holding all
#' other inputs constant. Works with Ne_sexual_Y2000 and Ne_mixed_Y2000.
#'
#' @param model_fn   The model function to use: Ne_sexual_Y2000 or
#'                   Ne_mixed_Y2000. Do NOT supply Ne_clonal_Y2000 — sexual
#'                   reproductive variance does not enter that model.
#' @param T_mat      Survival/transition matrix (s x s).
#' @param F_vec      Fecundity vector (length s).
#' @param D          Stage frequency distribution (length s, sums to 1).
#' @param d          Per-stage clonal fraction vector. Required if
#'                   model_fn = Ne_mixed_Y2000; ignored otherwise.
#' @param stage_index  Integer. Which stage to vary. Takes priority over
#'                   stage_name if both are supplied.
#' @param stage_name   Character. Column name of the focal stage in T_mat.
#'                   Used only if stage_index is NULL.
#' @param Vk_range   Numeric vector of (Vk/k_bar) values to test.
#'                   Default seq(0.5, 5, by = 0.5).
#' @param Vk_fixed   Numeric. (Vk/k_bar) for all OTHER stages. Default 1.
#' @param Vc_over_c  Clonal variance vector (for Ne_mixed_Y2000 only).
#'                   Default rep(1, s).
#' @param a          Hardy-Weinberg deviation. Default 0.
#' @param L          Generation time. If NULL, computed internally.
#' @param Ne_target  Ne conservation threshold. Default 5000.
#' @param population Character label.
#'
#' @return A list with:
#'   \describe{
#'     \item{data}{Data frame: Vk_over_k, NeN, Min_N, L, V, V_term1, V_term2,
#'                 and V_term3 (for mixed model)}
#'     \item{plot}{ggplot2 object showing Ne/N vs Vk/k_bar}
#'     \item{call}{Matched call}
#'   }
#'
#' @examples
#' T_herb <- matrix(c(
#'   0.30, 0.05, 0.00,
#'   0.40, 0.65, 0.10,
#'   0.00, 0.20, 0.80
#' ), nrow = 3, byrow = TRUE)
#'
#' sens <- Ne_sensitivity_Vk(
#'   model_fn    = Ne_sexual_Y2000,
#'   T_mat       = T_herb,
#'   F_vec       = c(0.0, 0.5, 3.0),
#'   D           = c(0.60, 0.25, 0.15),
#'   stage_index = 3,
#'   Vk_range    = seq(0.5, 6, by = 0.5),
#'   population  = "hypothetical herb"
#' )
#' print(sens$data)
#' sens$plot
#'
#' @export
Ne_sensitivity_Vk <- function(model_fn,
                               T_mat,
                               F_vec,
                               D,
                               d           = NULL,
                               stage_index = NULL,
                               stage_name  = NULL,
                               Vk_range    = seq(0.5, 5, by = 0.5),
                               Vk_fixed    = 1,
                               Vc_over_c   = NULL,
                               a           = 0,
                               L           = NULL,
                               Ne_target   = 5000,
                               population  = NULL) {

  s      <- nrow(T_mat)
  focal  <- .resolve_stage(stage_index, stage_name, s, T_mat)
  s_name <- if (!is.null(colnames(T_mat))) colnames(T_mat)[focal] else
    paste0("stage_", focal)

  # Check model_fn is not the clonal model
  fn_name <- deparse(substitute(model_fn))
  if (grepl("clonal", fn_name))
    stop("Ne_sensitivity_Vk does not apply to Ne_clonal_Y2000.\n",
         "  Sexual reproductive variance (Vk/k_bar) does not enter the\n",
         "  clonal model. Use Ne_sensitivity_L() for clonal populations.")

  results <- lapply(Vk_range, function(vk) {
    Vk_vec        <- rep(Vk_fixed, s)
    Vk_vec[focal] <- vk

    args <- list(T_mat     = T_mat,
                 F_vec     = F_vec,
                 D         = D,
                 Vk_over_k = Vk_vec,
                 a         = a,
                 L         = L,
                 Ne_target = Ne_target,
                 population = population)

    if (!is.null(d))         args$d         <- d
    if (!is.null(Vc_over_c)) args$Vc_over_c <- Vc_over_c

    r <- do.call(model_fn, args)

    # Build row defensively — field names differ across models:
    # clonal: no V stored directly (derive from u2_bar)
    # sexual: V, V_component1, V_component2
    # mixed:  V, V_term1, V_term2, V_term3
    V_val <- if (!is.null(r$V)) r$V else 2 * (1 - r$u2_bar)
    row <- data.frame(
      Vk_over_k = vk,
      NeN       = r$NeN,
      Min_N     = r$Min_N,
      L         = r$L,
      V         = V_val
    )
    if (!is.null(r$V_term1))      row$V_term1      <- r$V_term1
    if (!is.null(r$V_term2))      row$V_term2      <- r$V_term2
    if (!is.null(r$V_term3))      row$V_term3      <- r$V_term3
    if (!is.null(r$V_component1)) row$V_component1 <- r$V_component1
    if (!is.null(r$V_component2)) row$V_component2 <- r$V_component2
    row
  })

  df <- do.call(rbind, results)

  # Reference value: result at Vk = 1 (Poisson)
  ref_NeN <- df$NeN[which.min(abs(df$Vk_over_k - 1))]

  p <- ggplot2::ggplot(df, ggplot2::aes(x = Vk_over_k, y = NeN)) +
    ggplot2::geom_line(color = "#2166ac", linewidth = 1) +
    ggplot2::geom_point(color = "#2166ac", size = 2.5) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed",
                        color = "grey50", linewidth = 0.7) +
    ggplot2::annotate("text", x = 1, y = max(df$NeN),
                      label = "Poisson\ndefault",
                      hjust = -0.1, vjust = 1, size = 3.2,
                      color = "grey40") +
    ggplot2::geom_hline(yintercept = ref_NeN, linetype = "dotted",
                        color = "#d6604d", linewidth = 0.6) +
    ggplot2::labs(
      title    = paste0("Sensitivity of Ne/N to sexual reproductive variance"),
      subtitle = paste0("Focal stage: ", s_name,
                        if (!is.null(population)) paste0("  |  Population: ", population) else ""),
      x        = expression(paste("Variance-to-mean ratio of sexual output  ", (V[k]/bar(k))[i])),
      y        = expression(N[e]/N),
      caption  = paste0("Other stages fixed at Vk/k_bar = ", Vk_fixed,
                        ".  Red dotted line = Ne/N at Poisson default.")
    ) +
    .ne_sensitivity_theme()

  out <- list(
    data        = df,
    plot        = p,
    call        = match.call(),
    param       = "Vk_over_k",
    param_label = "Vk/k_bar (sexual reproductive variance)",
    ref_value   = 1,
    ref_label   = "Poisson default",
    stage_name  = s_name,
    fn_name     = fn_name,
    population  = if (!is.null(population)) population else "unnamed"
  )
  class(out) <- "Ne_sensitivity"
  out
}


# =============================================================================
# FUNCTION 2: Ne_sensitivity_Vc
# =============================================================================
# Sweeps (Vc/c_bar)_i — the variance-to-mean ratio of CLONAL reproductive
# output — for a single focal stage.
#
# Applies to: Ne_mixed_Y2000 ONLY.
#
# Biological question:
#   How sensitive is Ne to unequal clonal output in the focal stage?
#   Recall that under Poisson defaults, (Vc/c_bar = 1 = Vk/k_bar), the
#   clonal fraction d_i has no effect on Ne. As soon as Vc/c_bar > 1
#   (some individuals produce many more clonal propagules than others),
#   higher d_i begins to reduce Ne. This function quantifies that effect.

#' Sensitivity of Ne/N to clonal reproductive variance in one stage
#'
#' Sweeps (Vc/c_bar)_i across a range for a single focal stage. Mixed model
#' only. Under Poisson defaults (Vc/c_bar = 1 = Vk/k_bar), clonal fraction
#' d_i has no effect on Ne. Non-Poisson clonal variance (Vc/c_bar > 1) is
#' required for d_i to influence Ne.
#'
#' @param T_mat      Survival/transition matrix (s x s).
#' @param F_vec      Fecundity vector (length s).
#' @param D          Stage frequency distribution (length s, sums to 1).
#' @param d          Per-stage clonal fraction vector (length s). Required.
#' @param stage_index  Integer. Which stage to vary. Takes priority over
#'                   stage_name if both are supplied.
#' @param stage_name   Character. Column name of the focal stage in T_mat.
#' @param Vc_range   Numeric vector of (Vc/c_bar) values to test.
#'                   Default seq(0.5, 5, by = 0.5).
#' @param Vc_fixed   Numeric. (Vc/c_bar) for all OTHER stages. Default 1.
#' @param Vk_over_k  Sexual variance vector. Default rep(1, s).
#' @param a          Hardy-Weinberg deviation. Default 0.
#' @param L          Generation time. If NULL, computed internally.
#' @param Ne_target  Ne conservation threshold. Default 5000.
#' @param population Character label.
#'
#' @return A list with:
#'   \describe{
#'     \item{data}{Data frame: Vc_over_c, NeN, Min_N, L, V, V_term1,
#'                 V_term2, V_term3}
#'     \item{plot}{ggplot2 object showing Ne/N vs Vc/c_bar}
#'     \item{call}{Matched call}
#'   }
#'
#' @examples
#' T_herb <- matrix(c(
#'   0.30, 0.05, 0.00,
#'   0.40, 0.65, 0.10,
#'   0.00, 0.20, 0.80
#' ), nrow = 3, byrow = TRUE)
#'
#' sens <- Ne_sensitivity_Vc(
#'   T_mat       = T_herb,
#'   F_vec       = c(0.0, 0.5, 3.0),
#'   D           = c(0.60, 0.25, 0.15),
#'   d           = c(0.0, 0.0, 0.7),
#'   stage_index = 3,
#'   Vc_range    = seq(0.5, 6, by = 0.5),
#'   population  = "mixed herb"
#' )
#' print(sens$data)
#' sens$plot
#'
#' @export
Ne_sensitivity_Vc <- function(T_mat,
                               F_vec,
                               D,
                               d,
                               stage_index = NULL,
                               stage_name  = NULL,
                               Vc_range    = seq(0.5, 5, by = 0.5),
                               Vc_fixed    = 1,
                               Vk_over_k   = NULL,
                               a           = 0,
                               L           = NULL,
                               Ne_target   = 5000,
                               population  = NULL) {

  s      <- nrow(T_mat)
  focal  <- .resolve_stage(stage_index, stage_name, s, T_mat)
  s_name <- if (!is.null(colnames(T_mat))) colnames(T_mat)[focal] else
    paste0("stage_", focal)

  if (is.null(Vk_over_k)) Vk_over_k <- rep(1, s)

  results <- lapply(Vc_range, function(vc) {
    Vc_vec        <- rep(Vc_fixed, s)
    Vc_vec[focal] <- vc

    r <- Ne_mixed_Y2000(
      T_mat      = T_mat,
      F_vec      = F_vec,
      D          = D,
      d          = d,
      Vk_over_k  = Vk_over_k,
      Vc_over_c  = Vc_vec,
      a          = a,
      L          = L,
      Ne_target  = Ne_target,
      population = population
    )

    data.frame(
      Vc_over_c = vc,
      NeN       = r$NeN,
      Min_N     = r$Min_N,
      L         = r$L,
      V         = r$V,
      V_term1   = r$V_term1,
      V_term2   = r$V_term2,
      V_term3   = r$V_term3
    )
  })

  df <- do.call(rbind, results)

  ref_NeN <- df$NeN[which.min(abs(df$Vc_over_c - 1))]

  p <- ggplot2::ggplot(df, ggplot2::aes(x = Vc_over_c, y = NeN)) +
    ggplot2::geom_line(color = "#1b7837", linewidth = 1) +
    ggplot2::geom_point(color = "#1b7837", size = 2.5) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed",
                        color = "grey50", linewidth = 0.7) +
    ggplot2::annotate("text", x = 1, y = max(df$NeN),
                      label = "Poisson\ndefault",
                      hjust = -0.1, vjust = 1, size = 3.2,
                      color = "grey40") +
    ggplot2::geom_hline(yintercept = ref_NeN, linetype = "dotted",
                        color = "#d6604d", linewidth = 0.6) +
    ggplot2::labs(
      title    = "Sensitivity of Ne/N to clonal reproductive variance",
      subtitle = paste0("Focal stage: ", s_name,
                        if (!is.null(population)) paste0("  |  Population: ", population) else ""),
      x        = expression(paste("Variance-to-mean ratio of clonal output  ", (V[c]/bar(c))[i])),
      y        = expression(N[e]/N),
      caption  = paste0("Other stages fixed at Vc/c_bar = ", Vc_fixed,
                        ".  Note: V_term3 = 0 when Vc/c_bar = Vk/k_bar = 1 (Poisson).\n",
                        "  Red dotted line = Ne/N at Poisson default.")
    ) +
    .ne_sensitivity_theme()

  out <- list(
    data        = df,
    plot        = p,
    call        = match.call(),
    param       = "Vc_over_c",
    param_label = "Vc/c_bar (clonal reproductive variance)",
    ref_value   = 1,
    ref_label   = "Poisson default",
    stage_name  = s_name,
    fn_name     = "Ne_mixed_Y2000",
    population  = if (!is.null(population)) population else "unnamed"
  )
  class(out) <- "Ne_sensitivity"
  out
}


# =============================================================================
# FUNCTION 3: Ne_sensitivity_d
# =============================================================================
# Sweeps d_i — the clonal reproduction fraction — for one focal stage from
# 0 (fully sexual) to 1 (fully clonal).
#
# Applies to: Ne_mixed_Y2000 ONLY.
#
# Biological question:
#   How much does the clonal fraction in the focal stage matter for Ne?
#   Under Poisson defaults the answer is: not at all (V_term3 = 0).
#   Under non-Poisson clonal variance the answer is: it depends on how
#   much Vc/c_bar exceeds Vk/k_bar.
#   This function makes that relationship explicit and quantifies it.
#
# Note: supplying Vc_over_c != 1 is needed to see any effect on Ne.
#   Under Poisson defaults all d values give the same Ne/N.

#' Sensitivity of Ne/N to clonal fraction d_i in one stage
#'
#' Sweeps d_i from 0 to 1 for a single focal stage, holding all other inputs
#' constant. Under Poisson defaults (Vc/c_bar = Vk/k_bar = 1), Ne/N is
#' constant across all d values. Supply Vc_over_c != 1 to see a meaningful
#' sensitivity curve.
#'
#' @param T_mat      Survival/transition matrix (s x s).
#' @param F_vec      Fecundity vector (length s).
#' @param D          Stage frequency distribution (length s, sums to 1).
#' @param d_fixed    Numeric vector (length s). Clonal fractions for ALL
#'                   stages. The focal stage value will be overwritten by
#'                   the sweep. Must be length s.
#' @param stage_index  Integer. Which stage to vary. Takes priority over
#'                   stage_name if both are supplied.
#' @param stage_name   Character. Column name of the focal stage in T_mat.
#' @param d_range    Numeric vector of d values to test for focal stage.
#'                   Default seq(0, 1, by = 0.1).
#' @param Vk_over_k  Sexual variance vector. Default rep(1, s).
#' @param Vc_over_c  Clonal variance vector. Default rep(1, s).
#'                   Supply values > 1 to see a non-flat sensitivity curve.
#' @param a          Hardy-Weinberg deviation. Default 0.
#' @param L          Generation time. If NULL, computed internally.
#' @param Ne_target  Ne conservation threshold. Default 5000.
#' @param population Character label.
#'
#' @return A list with:
#'   \describe{
#'     \item{data}{Data frame: d_focal, NeN, Min_N, L, V, V_term1, V_term2,
#'                 V_term3}
#'     \item{plot}{ggplot2 object showing Ne/N vs d_i}
#'     \item{call}{Matched call}
#'   }
#'
#' @examples
#' T_herb <- matrix(c(
#'   0.30, 0.05, 0.00,
#'   0.40, 0.65, 0.10,
#'   0.00, 0.20, 0.80
#' ), nrow = 3, byrow = TRUE)
#'
#' # Under Poisson defaults: flat curve (d_i has no effect)
#' sens_poisson <- Ne_sensitivity_d(
#'   T_mat       = T_herb,
#'   F_vec       = c(0.0, 0.5, 3.0),
#'   D           = c(0.60, 0.25, 0.15),
#'   d_fixed     = c(0.0, 0.0, 0.5),
#'   stage_index = 3,
#'   population  = "mixed herb (Poisson)"
#' )
#' sens_poisson$plot
#'
#' # With non-Poisson clonal variance: d_i matters
#' sens_nonpoisson <- Ne_sensitivity_d(
#'   T_mat       = T_herb,
#'   F_vec       = c(0.0, 0.5, 3.0),
#'   D           = c(0.60, 0.25, 0.15),
#'   d_fixed     = c(0.0, 0.0, 0.5),
#'   stage_index = 3,
#'   Vc_over_c   = c(1, 1, 3),
#'   population  = "mixed herb (Vc/c_bar = 3 in stage 3)"
#' )
#' sens_nonpoisson$plot
#'
#' @export
Ne_sensitivity_d <- function(T_mat,
                              F_vec,
                              D,
                              d_fixed,
                              stage_index = NULL,
                              stage_name  = NULL,
                              d_range     = seq(0, 1, by = 0.1),
                              Vk_over_k   = NULL,
                              Vc_over_c   = NULL,
                              a           = 0,
                              L           = NULL,
                              Ne_target   = 5000,
                              population  = NULL) {

  s      <- nrow(T_mat)
  focal  <- .resolve_stage(stage_index, stage_name, s, T_mat)
  s_name <- if (!is.null(colnames(T_mat))) colnames(T_mat)[focal] else
    paste0("stage_", focal)

  if (!is.numeric(d_fixed) || length(d_fixed) != s)
    stop(paste0("d_fixed must be a numeric vector of length ", s,
                " (the base clonal fractions; focal stage value will be swept)."))

  if (is.null(Vk_over_k)) Vk_over_k <- rep(1, s)
  if (is.null(Vc_over_c)) Vc_over_c <- rep(1, s)

  # Warn user if Poisson defaults will produce a flat curve
  if (all(Vc_over_c == 1) && all(Vk_over_k == 1) && a == 0)
    message(
      "Note: Under Poisson defaults (Vc/c_bar = Vk/k_bar = 1, a = 0),\n",
      "  Ne/N is constant across all d values (V_term3 = 0).\n",
      "  Supply Vc_over_c != 1 to see a meaningful sensitivity curve."
    )

  results <- lapply(d_range, function(dval) {
    d_vec        <- d_fixed
    d_vec[focal] <- dval

    r <- Ne_mixed_Y2000(
      T_mat      = T_mat,
      F_vec      = F_vec,
      D          = D,
      d          = d_vec,
      Vk_over_k  = Vk_over_k,
      Vc_over_c  = Vc_over_c,
      a          = a,
      L          = L,
      Ne_target  = Ne_target,
      population = population
    )

    data.frame(
      d_focal = dval,
      NeN     = r$NeN,
      Min_N   = r$Min_N,
      L       = r$L,
      V       = r$V,
      V_term1 = r$V_term1,
      V_term2 = r$V_term2,
      V_term3 = r$V_term3
    )
  })

  df <- do.call(rbind, results)

  # Mark the supplied d_fixed value for the focal stage
  d_original <- d_fixed[focal]

  p <- ggplot2::ggplot(df, ggplot2::aes(x = d_focal, y = NeN)) +
    ggplot2::geom_line(color = "#762a83", linewidth = 1) +
    ggplot2::geom_point(color = "#762a83", size = 2.5) +
    ggplot2::geom_vline(xintercept = d_original, linetype = "dashed",
                        color = "grey50", linewidth = 0.7) +
    ggplot2::annotate("text", x = d_original, y = max(df$NeN),
                      label = paste0("d = ", d_original, "\n(base value)"),
                      hjust = -0.08, vjust = 1, size = 3.2,
                      color = "grey40") +
    ggplot2::scale_x_continuous(breaks = seq(0, 1, by = 0.2),
                                 limits = c(0, 1)) +
    ggplot2::labs(
      title    = expression(paste("Sensitivity of ", N[e]/N, " to clonal fraction ", d[i])),
      subtitle = paste0("Focal stage: ", s_name,
                        if (!is.null(population)) paste0("  |  Population: ", population) else ""),
      x        = expression(paste("Clonal fraction  ", d[i],
                                  "  (0 = fully sexual,  1 = fully clonal)")),
      y        = expression(N[e]/N),
      caption  = paste0("Other stage d values fixed. Dashed line = base d value (",
                        d_original, ").\n",
                        "  Flat curve = Poisson defaults active (supply Vc_over_c != 1 to see effect).")
    ) +
    .ne_sensitivity_theme()

  out <- list(
    data        = df,
    plot        = p,
    call        = match.call(),
    param       = "d_focal",
    param_label = "d (clonal fraction)",
    ref_value   = d_original,
    ref_label   = "base d value",
    stage_name  = s_name,
    fn_name     = "Ne_mixed_Y2000",
    population  = if (!is.null(population)) population else "unnamed"
  )
  class(out) <- "Ne_sensitivity"
  out
}


# =============================================================================
# FUNCTION 4: Ne_sensitivity_L
# =============================================================================
# Sweeps generation time L across a range for any of the three model
# functions. This is the most universally applicable sensitivity function
# because L must always be supplied or computed, and its estimation is
# often a source of uncertainty.
#
# Biological question:
#   How sensitive is Ne to the assumed or estimated generation time?
#   If L is uncertain (e.g., computed from a single year of data, or
#   taken from a different population), how much does this affect the
#   minimum census size recommendation?

#' Sensitivity of Ne/N to generation time L
#'
#' Sweeps generation time L across a user-defined range for any of the three
#' NeStage model functions. All other inputs are held constant.
#'
#' @param model_fn   The model function: Ne_clonal_Y2000, Ne_sexual_Y2000,
#'                   or Ne_mixed_Y2000.
#' @param T_mat      Survival/transition matrix (s x s).
#' @param F_vec      Fecundity vector (length s).
#' @param D          Stage frequency distribution (length s, sums to 1).
#' @param d          Per-stage clonal fraction vector. Required for
#'                   Ne_mixed_Y2000; ignored otherwise.
#' @param L_range    Numeric vector of L values (years) to test.
#'                   Default seq(2, 30, by = 2).
#' @param L_ref      Numeric. Reference L value to mark on the plot
#'                   (e.g., the computed or published L). If NULL,
#'                   computed internally from T_mat and F_vec.
#' @param Vk_over_k  Sexual variance vector. Default rep(1, s).
#' @param Vc_over_c  Clonal variance vector (mixed model only).
#' @param a          Hardy-Weinberg deviation. Default 0.
#' @param Ne_target  Ne conservation threshold. Default 5000.
#' @param population Character label.
#'
#' @return A list with:
#'   \describe{
#'     \item{data}{Data frame: L, NeN, Min_N, V, V_term1, V_term2,
#'                 and V_term3 (for mixed model)}
#'     \item{plot}{ggplot2 object showing Ne/N vs L}
#'     \item{call}{Matched call}
#'   }
#'
#' @examples
#' # Clonal model — how sensitive is Ne/N to L for Fritillaria Miz?
#' T_Miz <- matrix(c(
#'   0.789, 0.121, 0.054,
#'   0.007, 0.621, 0.335,
#'   0.001, 0.258, 0.611
#' ), nrow = 3, byrow = TRUE)
#'
#' sens <- Ne_sensitivity_L(
#'   model_fn   = Ne_clonal_Y2000,
#'   T_mat      = T_Miz,
#'   F_vec      = c(0.055, 1.328, 2.398),
#'   D          = c(0.935, 0.038, 0.027),
#'   L_range    = seq(5, 25, by = 1),
#'   L_ref      = 13.399,
#'   population = "Fritillaria Miz"
#' )
#' print(sens$data)
#' sens$plot
#'
#' @export
Ne_sensitivity_L <- function(model_fn,
                              T_mat,
                              F_vec,
                              D,
                              d           = NULL,
                              L_range     = seq(2, 30, by = 2),
                              L_ref       = NULL,
                              Vk_over_k   = NULL,
                              Vc_over_c   = NULL,
                              a           = 0,
                              Ne_target   = 5000,
                              population  = NULL) {

  s       <- nrow(T_mat)
  fn_name <- deparse(substitute(model_fn))

  if (is.null(Vk_over_k)) Vk_over_k <- rep(1, s)

  results <- lapply(L_range, function(lval) {

    args <- list(T_mat      = T_mat,
                 F_vec      = F_vec,
                 D          = D,
                 L          = lval,
                 Ne_target  = Ne_target,
                 population = population)

    # Add model-specific arguments
    if (grepl("sexual|mixed", fn_name)) args$Vk_over_k <- Vk_over_k
    if (grepl("mixed",  fn_name)) {
      if (!is.null(d))         args$d         <- d
      if (!is.null(Vc_over_c)) args$Vc_over_c <- Vc_over_c
    }

    r <- do.call(model_fn, args)

    # Ne_clonal_Y2000 does not store V explicitly — derive it from
    # V = 2(1 - u2_bar), the clonal model formula (Eq. 10).
    # Sexual and mixed models store V directly in r$V.
    V_val <- if (!is.null(r$V)) r$V else 2 * (1 - r$u2_bar)

    row <- data.frame(
      L     = lval,
      NeN   = r$NeN,
      Min_N = r$Min_N,
      V     = V_val
    )
    if (!is.null(r$V_term1)) row$V_term1 <- r$V_term1
    if (!is.null(r$V_term2)) row$V_term2 <- r$V_term2
    if (!is.null(r$V_term3)) row$V_term3 <- r$V_term3
    row
  })

  df <- do.call(rbind, results)

  # Compute reference L if not supplied
  if (is.null(L_ref)) {
    s_     <- nrow(T_mat)
    F_mat  <- matrix(0, s_, s_); F_mat[1, ] <- F_vec
    A_     <- T_mat + F_mat
    ev_    <- eigen(A_)
    w_     <- Re(ev_$vectors[, which.max(Re(ev_$values))])
    w_     <- abs(w_) / sum(abs(w_))
    Tx_    <- diag(s_); num_ <- 0; den_ <- 0
    for (x in seq_len(500L)) {
      Tx_ <- T_mat %*% Tx_
      l_  <- 0; m_ <- 0
      for (i in seq_len(s_)) {
        u_ <- Tx_[, i]
        l_ <- l_ + w_[i] * sum(u_)
        m_ <- m_ + w_[i] * sum(F_vec * u_)
      }
      num_ <- num_ + x * m_ * l_
      den_ <- den_ +     m_ * l_
    }
    L_ref <- if (den_ > 0) as.numeric(num_ / den_) else NA_real_
  }

  # Model label for plot
  model_label <- if (grepl("clonal", fn_name)) "Clonal model (Eq. 10)" else
    if (grepl("sexual", fn_name)) "Sexual model (Eq. 6, d=0)" else
    if (grepl("mixed",  fn_name)) "Mixed model (Eq. 6)" else fn_name

  line_color <- if (grepl("clonal", fn_name)) "#d6604d" else
    if (grepl("sexual", fn_name)) "#2166ac" else "#762a83"

  p <- ggplot2::ggplot(df, ggplot2::aes(x = L, y = NeN)) +
    ggplot2::geom_line(color = line_color, linewidth = 1) +
    ggplot2::geom_point(color = line_color, size = 2.5) +
    ggplot2::labs(
      title    = expression(paste("Sensitivity of ", N[e]/N, " to generation time ", L)),
      subtitle = paste0("Model: ", fn_name,
                        if (!is.null(population)) paste0("  |  Population: ", population) else ""),
      x        = "Generation time  L  (years)",
      y        = expression(N[e]/N),
      caption  = paste0("V is held constant across L values.\n",
                        "  Ne/N = 2 / (V * L), so the relationship is hyperbolic.")
    ) +
    .ne_sensitivity_theme()

  # Add reference L line if available
  if (!is.na(L_ref) && L_ref >= min(L_range) && L_ref <= max(L_range)) {
    ref_NeN <- df$NeN[which.min(abs(df$L - L_ref))]
    p <- p +
      ggplot2::geom_vline(xintercept = L_ref, linetype = "dashed",
                          color = "grey50", linewidth = 0.7) +
      ggplot2::annotate("text", x = L_ref, y = max(df$NeN),
                        label = paste0("L = ", round(L_ref, 2), "\n(reference)"),
                        hjust = -0.08, vjust = 1, size = 3.2,
                        color = "grey40") +
      ggplot2::geom_hline(yintercept = ref_NeN, linetype = "dotted",
                          color = "grey50", linewidth = 0.6)
  }

  out <- list(
    data        = df,
    plot        = p,
    call        = match.call(),
    param       = "L",
    param_label = "L (generation time, years)",
    ref_value   = L_ref,
    ref_label   = "reference L",
    stage_name  = NA,
    fn_name     = fn_name,
    population  = if (!is.null(population)) population else "unnamed"
  )
  class(out) <- "Ne_sensitivity"
  out
}


# =============================================================================
# S3 PRINT METHOD: print.Ne_sensitivity
# =============================================================================
# Prints a full interpretive summary when the user calls print(sens).
# Also displays the plot automatically.
#
# Elasticity methodology:
#   Caswell H. (2001). Matrix Population Models, 2nd ed. Sinauer Associates.
#   Morris W.F. & Doak D.F. (2002). Quantitative Conservation Biology.
#     Sinauer Associates. Chapter 9.
#
# Local elasticity (at the reference point):
#   E_local = (d NeN / d theta) * (theta_ref / NeN_ref)
#   Approximated by finite difference between the two points bracketing
#   the reference value in the swept range.
#
# Global elasticity (mean across full range):
#   E_global = mean( |delta NeN / NeN| / |delta theta / theta| )
#   Computed as the mean of consecutive proportional changes.
#
# Empirical Ne/N benchmarks:
#   Frankham R. (1995). Effective population size/adult population size
#   ratios in wildlife: a review. Genetical Research 66: 95-107.
#   Reports median Ne/N ~ 0.10-0.11 across 102 species (range 0.01-0.99).

#' @export
print.Ne_sensitivity <- function(x, ...) {

  df  <- x$data
  col <- x$param       # name of the swept parameter column in df

  # ------------------------------------------------------------------
  # 1. Reference point values
  # ------------------------------------------------------------------
  ref_val  <- x$ref_value
  ref_idx  <- which.min(abs(df[[col]] - ref_val))
  ref_NeN  <- df$NeN[ref_idx]
  ref_MinN <- df$Min_N[ref_idx]

  # ------------------------------------------------------------------
  # 2. Range summary
  # ------------------------------------------------------------------
  max_NeN  <- max(df$NeN);  max_idx  <- which.max(df$NeN)
  min_NeN  <- min(df$NeN);  min_idx  <- which.min(df$NeN)
  max_par  <- df[[col]][max_idx]
  min_par  <- df[[col]][min_idx]
  max_MinN <- df$Min_N[max_idx]
  min_MinN <- df$Min_N[min_idx]

  # ------------------------------------------------------------------
  # 3. Local elasticity at the reference point
  # Method: Caswell (2001), Morris & Doak (2002) Ch. 9
  # E = (dNeN/dtheta) * (theta_ref / NeN_ref)
  # Approximated by central finite difference where possible,
  # or forward/backward difference at the boundary.
  # ------------------------------------------------------------------
  n <- nrow(df)
  if (ref_idx > 1 && ref_idx < n) {
    # Central difference
    dNeN   <- df$NeN[ref_idx + 1]  - df$NeN[ref_idx - 1]
    dtheta <- df[[col]][ref_idx + 1] - df[[col]][ref_idx - 1]
  } else if (ref_idx == 1) {
    dNeN   <- df$NeN[2]   - df$NeN[1]
    dtheta <- df[[col]][2] - df[[col]][1]
  } else {
    dNeN   <- df$NeN[n]   - df$NeN[n - 1]
    dtheta <- df[[col]][n] - df[[col]][n - 1]
  }
  E_local <- if (abs(dtheta) > 0 && ref_NeN > 0 && abs(ref_val) > 0)
    (dNeN / dtheta) * (ref_val / ref_NeN)
  else NA_real_

  # ------------------------------------------------------------------
  # 4. Global elasticity (mean proportional sensitivity across range)
  # ------------------------------------------------------------------
  prop_dNeN   <- abs(diff(df$NeN)    / df$NeN[-n])
  prop_dtheta <- abs(diff(df[[col]]) / df[[col]][-n])
  valid       <- prop_dtheta > 0
  E_global    <- if (any(valid))
    mean(prop_dNeN[valid] / prop_dtheta[valid])
  else NA_real_

  # ------------------------------------------------------------------
  # 5. Conservation ratio: Min_N at worst vs reference point
  # ------------------------------------------------------------------
  worst_MinN  <- max(df$Min_N)
  worst_par   <- df[[col]][which.max(df$Min_N)]
  cons_ratio  <- worst_MinN / ref_MinN

  # ------------------------------------------------------------------
  # 6. Frankham (1995) benchmark comparison
  # ------------------------------------------------------------------
  frankham_median <- 0.10  # median Ne/N across 102 species
  frankham_range  <- c(0.01, 0.99)

  bench_msg <- if (ref_NeN < frankham_median) {
    paste0("    Ne/N at reference (", round(ref_NeN, 3), ") is BELOW the ",
           "empirical median of ~0.10\n",
           "    across 102 wildlife species (Frankham 1995). This population\n",
           "    may be more genetically vulnerable than average.")
  } else if (ref_NeN < 3 * frankham_median) {
    paste0("    Ne/N at reference (", round(ref_NeN, 3), ") is near the ",
           "empirical median of ~0.10\n",
           "    across 102 wildlife species (Frankham 1995).")
  } else {
    paste0("    Ne/N at reference (", round(ref_NeN, 3), ") is ABOVE the ",
           "empirical median of ~0.10\n",
           "    across 102 wildlife species (Frankham 1995). This population\n",
           "    has relatively high effective size.")
  }

  # ------------------------------------------------------------------
  # 7. Print the summary
  # ------------------------------------------------------------------
  cat("\n")
  cat("=== NeStage Sensitivity Analysis ===\n")
  cat(sprintf("  Parameter swept  : %s\n", x$param_label))
  if (!is.na(x$stage_name))
    cat(sprintf("  Focal stage      : %s\n", x$stage_name))
  cat(sprintf("  Model            : %s\n", x$fn_name))
  cat(sprintf("  Population       : %s\n", x$population))
  cat("\n")

  cat("  --- Range swept ---\n")
  cat(sprintf("    %s from %.3g to %.3g  (%d values)\n",
              col, min(df[[col]]), max(df[[col]]), n))
  cat("\n")

  cat("  --- Ne/N summary ---\n")
  cat(sprintf("    At reference (%s = %s)  : Ne/N = %.3f  |  Min N = %s\n",
              col, round(ref_val, 3), ref_NeN,
              format(ref_MinN, big.mark = ",")))
  cat(sprintf("    Maximum Ne/N  (%s = %s)  : Ne/N = %.3f  |  Min N = %s\n",
              col, round(max_par, 3), max_NeN,
              format(max_MinN, big.mark = ",")))
  cat(sprintf("    Minimum Ne/N  (%s = %s)  : Ne/N = %.3f  |  Min N = %s\n",
              col, round(min_par, 3), min_NeN,
              format(min_MinN, big.mark = ",")))
  cat(sprintf("    Ne/N range across sweep  : %.3f\n", max_NeN - min_NeN))
  cat("\n")

  cat("  --- Elasticity (Caswell 2001; Morris & Doak 2002) ---\n")
  cat(sprintf("    Local elasticity at reference point : E = %s\n",
              if (is.na(E_local)) "undefined" else round(E_local, 3)))
  cat(sprintf("    Global elasticity (mean over range) : E = %s\n",
              if (is.na(E_global)) "undefined" else round(E_global, 3)))
  cat("\n")

  # Interpretation of elasticity
  if (!is.na(E_local)) {
    abs_E <- abs(E_local)
    cat("    Interpretation (local elasticity):\n")
    cat(sprintf("      A 10%% increase in %s at the reference point\n", col))
    cat(sprintf("      produces a %.1f%% change in Ne/N.\n",
                abs(E_local) * 10))
    cat(sprintf("      Direction: %s\n",
                if (E_local < 0) "NEGATIVE — higher variance reduces Ne/N (typical)"
                else "POSITIVE — higher value increases Ne/N (atypical)"))
    cat("\n")
  }

  cat("  --- Conservation implications ---\n")
  cat(sprintf("    Min N at reference point             : %s\n",
              format(ref_MinN, big.mark = ",")))
  cat(sprintf("    Min N at worst case (%s = %s)  : %s\n",
              col, round(worst_par, 3),
              format(worst_MinN, big.mark = ",")))
  cat(sprintf("    Ratio worst/reference                : %.1fx more individuals\n",
              cons_ratio))
  cat(sprintf("    Ne target used                       : %s\n",
              format(round(ref_MinN * ref_NeN), big.mark = ",", scientific = FALSE)))
  cat("\n")

  cat("  --- Comparison to empirical Ne/N (Frankham 1995) ---\n")
  cat(bench_msg, "\n")
  cat(sprintf("    Frankham (1995) empirical range: %.2f -- %.2f\n",
              frankham_range[1], frankham_range[2]))
  cat("\n")

  cat("  --- References ---\n")
  cat("    Caswell H. (2001). Matrix Population Models, 2nd ed.\n")
  cat("      Sinauer Associates. [elasticity methodology]\n")
  cat("    Morris W.F. & Doak D.F. (2002). Quantitative Conservation\n")
  cat("      Biology. Sinauer Associates. Ch. 9. [sensitivity analysis]\n")
  cat("    Frankham R. (1995). Effective population size/adult population\n")
  cat("      size ratios in wildlife: a review.\n")
  cat("      Genetical Research 66: 95-107. [Ne/N empirical benchmarks]\n")
  cat("    Yonezawa et al. (2000). Evolution 54(6): 2007-2013. [model]\n")
  cat("\n")

  # ------------------------------------------------------------------
  # 8. Display the plot automatically
  # ------------------------------------------------------------------
  print(x$plot)

  invisible(x)
}
