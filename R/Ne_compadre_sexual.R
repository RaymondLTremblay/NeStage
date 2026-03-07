# =============================================================================
# Ne_compadre_sexual.R
# -----------------------------------------------------------------------------
# Two functions for computing variance effective population size (Ne) across
# sexually reproducing plant species using COMPADRE matrix population models.
#
#   Ne_compadre_sexual_filter()  -- fetch + filter COMPADRE, one species/row
#   Ne_compadre_sexual()         -- run Ne_sexual_Y2000() across the filtered db
#
# Workflow:
#   db  <- Ne_compadre_sexual_filter()          # requires internet, ~1 min
#   out <- Ne_compadre_sexual(db)               # Ne analysis, ~seconds
#   out$summary                                 # tidy results table
#   out$results[["Cirsium pitcheri"]]           # full result for one species
#
# Dependencies (all on CRAN unless noted):
#   Rcompadre  -- COMPADRE access and matrix accessors
#   NeStage    -- Ne_sexual_Y2000()   [github.com/RaymondLTremblay/NeStage]
#
# Author:  Raymond L. Tremblay  <raymond.tremblay@upr.edu>
# ORCID:   0000-0002-8588-4372
# Version: 1.1.0  (2026-03-06)
# Reference:
#   Yonezawa et al. (2000) Evolution 54(6): 2007-2013.
#   doi:10.1111/j.0014-3820.2000.tb01243.x
# =============================================================================


# Column names used in subset() calls on CompadreDB objects.
# These are database metadata columns, not R variables -- suppress R CMD check notes.
utils::globalVariables(c(
  "MatrixComposite", "MatrixTreatment", "MatrixCaptivity",
  "AnnualPeriodicity", "MatrixSplit", "MatrixFec", "Kingdom",
  "check_NA_U", "check_NA_F", "check_NA_C",
  "check_zero_U", "check_zero_F"
))

# Column names used in subset() calls on CompadreDB objects.
# These are database metadata columns, not R variables -- suppress R CMD check notes.
utils::globalVariables(c(
  "MatrixComposite", "MatrixTreatment", "MatrixCaptivity",
  "AnnualPeriodicity", "MatrixSplit", "MatrixFec", "Kingdom",
  "check_NA_U", "check_NA_F", "check_NA_C",
  "check_zero_U", "check_zero_F"
))

# =============================================================================
# PRIVATE HELPERS
# (prefixed with . so they do not pollute the user namespace)
# =============================================================================

# Stable stage distribution D from A = matU + matF
# Returns a vector summing to 1 (dominant right eigenvector, absolute values).
.stable_stage_D <- function(matU, matF) {
  A  <- matU + matF
  ev <- eigen(A)
  # dominant eigenvalue index
  idx <- which.max(Re(ev$values))
  w   <- Re(ev$vectors[, idx])
  w   <- abs(w)
  if (sum(w) == 0) stop("eigenvector sums to zero")
  w / sum(w)
}

# TRUE if every element of matC is zero (NA-safe: NA treated as non-zero)
.is_pure_sexual <- function(matC) {
  if (any(is.na(matC))) return(FALSE)
  all(matC == 0)
}


# =============================================================================
# FUNCTION 1: Ne_compadre_sexual_filter
# =============================================================================

#' Filter COMPADRE for sexually reproducing plant species
#'
#' Fetches the full COMPADRE Plant Matrix Database and applies a sequential
#' set of biological and quality filters. Returns one mean matrix per species
#' (the first available if multiple populations exist), ready for use with
#' \code{\link{Ne_compadre_sexual}}.
#'
#' @section Filters applied (in order):
#' \enumerate{
#'   \item \code{MatrixComposite == "Mean"} -- mean matrices only
#'   \item \code{MatrixTreatment == "Unmanipulated"} -- no experimental treatment
#'   \item \code{MatrixCaptivity == "W"} -- wild populations only
#'   \item \code{AnnualPeriodicity == 1} -- annual projection interval
#'   \item \code{MatrixSplit == "Divided"} -- U, F, C submatrices available
#'   \item \code{MatrixFec == "Yes"} -- fecundity was measured
#'   \item \code{Kingdom == "Plantae"} -- plants only
#'   \item No NA values in matU, matF, or matC (\code{cdb_flag})
#'   \item No all-zero U column sums (\code{check_zero_U})
#'   \item No fecundity in matF (\code{check_zero_F} removes these)
#'   \item \code{matC} all zeros -- purely sexual (no clonal reproduction)
#'   \item One row per species (first mean matrix retained if duplicates)
#' }
#'
#' @param verbose Logical. Print step-by-step filtering counts and a family
#'   summary? Default \code{TRUE}.
#'
#' @return A \code{CompadreDB} object with one row per species. Attach it to
#'   \code{\link{Ne_compadre_sexual}} to compute Ne. You can also inspect the
#'   object directly -- it behaves like a data frame with an extra \code{mat}
#'   list column.
#'
#' @examples
#' \dontrun{
#' db <- Ne_compadre_sexual_filter()
#' nrow(db)              # species retained
#' table(db$Family)      # species per family
#' }
#'
#' @seealso \code{\link{Ne_compadre_sexual}}
#' @export
Ne_compadre_sexual_filter <- function(verbose = TRUE) {

  # ---- package checks -------------------------------------------------------
  for (pkg in c("Rcompadre")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Package '", pkg, "' is required. Install with:\n",
           "  install.packages('", pkg, "')")
  }

  msg <- function(...) if (verbose) message(...)

  # ---- Step 1: fetch --------------------------------------------------------
  msg("Fetching COMPADRE database (requires internet connection)...")
  db <- Rcompadre::cdb_fetch("compadre")
  msg("  Total matrices downloaded: ", nrow(db))

  # ---- Step 2: metadata filters ---------------------------------------------
  db <- subset(db, MatrixComposite  == "Mean")
  msg("  After MatrixComposite == 'Mean':           ", nrow(db))

  db <- subset(db, MatrixTreatment  == "Unmanipulated")
  msg("  After MatrixTreatment == 'Unmanipulated':  ", nrow(db))

  db <- subset(db, MatrixCaptivity  == "W")
  msg("  After MatrixCaptivity == 'W':              ", nrow(db))

  # AnnualPeriodicity is stored as numeric 1 in modern Rcompadre
  db <- subset(db, AnnualPeriodicity == 1)
  msg("  After AnnualPeriodicity == 1:              ", nrow(db))

  db <- subset(db, MatrixSplit == "Divided")
  msg("  After MatrixSplit == 'Divided':            ", nrow(db))

  db <- subset(db, MatrixFec == "Yes")
  msg("  After MatrixFec == 'Yes':                  ", nrow(db))

  db <- subset(db, Kingdom == "Plantae")
  msg("  After Kingdom == 'Plantae':                ", nrow(db))

  # ---- Step 3: matrix quality flags ----------------------------------------
  # cdb_flag() adds logical check_* columns; TRUE = problem found
  db <- Rcompadre::cdb_flag(db)

  db <- subset(db, check_NA_U    == FALSE)
  db <- subset(db, check_NA_F    == FALSE)
  db <- subset(db, check_NA_C    == FALSE)
  db <- subset(db, check_zero_U  == FALSE)   # no all-zero survival columns
  db <- subset(db, check_zero_F  == FALSE)   # matF must have non-zero entries
  msg("  After matrix quality flags (no NA, no zero-U, no zero-F): ", nrow(db))

  # ---- Step 4: purely sexual (matC all zeros) -------------------------------
  matC_list  <- Rcompadre::matC(db)
  is_sexual  <- vapply(matC_list, .is_pure_sexual, logical(1))
  db         <- db[is_sexual, ]
  msg("  After keeping matC == 0 (purely sexual):   ", nrow(db))

  # ---- Step 5: one row per species ------------------------------------------
  species_vec <- as.character(as.data.frame(db)$SpeciesAccepted)
  db          <- db[!duplicated(species_vec), ]
  msg("  After deduplication (one per species):     ", nrow(db), " species")

  # ---- Step 6: family summary -----------------------------------------------
  if (verbose) {
    fam_tab <- sort(table(as.character(db$Family)), decreasing = TRUE)
    message("\nFamily representation:")
    print(fam_tab)
    message("\n  Total species : ", nrow(db))
    message("  Total families: ", length(fam_tab))
  }

  db
}


# =============================================================================
# FUNCTION 2: Ne_compadre_sexual
# =============================================================================

#' Compute Ne/N for sexually reproducing plant species from COMPADRE
#'
#' Applies \code{Ne_sexual_Y2000()} to every species in a filtered
#' \code{CompadreDB} object. The stage frequency distribution \eqn{D} is
#' derived from the stable stage distribution of \eqn{A = matU + matF}
#' (the theoretical equilibrium distribution). Poisson reproductive variance
#' is assumed by default (\eqn{V_k / \bar{k} = 1}), which is the most
#' conservative estimate when individual-level variance data are unavailable.
#'
#' @param db A \code{CompadreDB} object produced by
#'   \code{\link{Ne_compadre_sexual_filter}}, or any \code{CompadreDB} object
#'   where each row has valid (no NA, non-zero) \code{matU} and \code{matF}.
#' @param Ne_target Numeric. Conservation Ne threshold for computing minimum
#'   census size \eqn{N_{min} = N_{e,target} / (N_e/N)}. Default \code{5000}.
#' @param Vk_over_k Numeric scalar or \code{NULL}. Ratio of variance to mean
#'   offspring number per stage. \code{NULL} (default) assumes Poisson
#'   (\code{Vk/k = 1}) for all stages. A value > 1 (e.g. 2) represents
#'   over-dispersed reproduction (higher variance reduces Ne).
#' @param verbose Logical. Print one line per species showing Ne/N and Min_N?
#'   Default \code{FALSE}.
#'
#' @return A named list with two elements:
#' \describe{
#'   \item{\code{summary}}{A \code{data.frame} with one row per species and
#'     columns: \code{Species}, \code{Family}, \code{Order}, \code{Country},
#'     \code{Lat}, \code{Lon}, \code{MatrixDimension}, \code{lambda},
#'     \code{NeN}, \code{Min_N}, \code{L}, \code{V}, \code{u_bar},
#'     \code{u2_bar}, \code{Ne_target}. Sorted by Family then Species.}
#'   \item{\code{results}}{A named list of full \code{Ne_sexual_Y2000} result
#'     objects, one per species. Names are the accepted species names.}
#' }
#'
#' @examples
#' \dontrun{
#' # Full workflow
#' db  <- Ne_compadre_sexual_filter()
#' out <- Ne_compadre_sexual(db, Ne_target = 5000)
#'
#' # Summary table
#' head(out$summary)
#'
#' # Individual species result
#' out$results[["Cirsium pitcheri"]]
#'
#' # Species requiring the most individuals for Ne >= 5000
#' out$summary[order(-out$summary$Min_N), c("Species", "Family", "NeN", "Min_N")]
#'
#' # Plot Ne/N by family
#' library(ggplot2)
#' ggplot(out$summary,
#'        aes(x = reorder(Family, NeN, median), y = NeN)) +
#'   geom_boxplot() + coord_flip() +
#'   labs(x = NULL, y = expression(N[e]/N))
#' }
#'
#' @seealso \code{\link{Ne_compadre_sexual_filter}}, \code{\link{Ne_sexual_Y2000}}
#' @export
Ne_compadre_sexual <- function(db,
                               Ne_target  = 5000,
                               Vk_over_k  = NULL,
                               verbose    = FALSE) {

  if (!requireNamespace("Rcompadre", quietly = TRUE))
    stop("Package 'Rcompadre' is required.")

  n <- nrow(db)
  if (n == 0L)
    stop("db has 0 rows. Run Ne_compadre_sexual_filter() first.")

  message("Running Ne_sexual_Y2000() on ", n, " species...")

  # Extract matrix lists and metadata once (avoids repeated S4 dispatch)
  matU_list <- Rcompadre::matU(db)
  matF_list <- Rcompadre::matF(db)
  meta      <- as.data.frame(db)

  # Detect Lat/Lon column name variants across COMPADRE releases
  lat_col <- intersect(c("Lat", "LatDeg"), names(meta))[1]
  lon_col <- intersect(c("Lon", "LonDeg"), names(meta))[1]

  results_list <- vector("list", n)
  summary_rows <- vector("list", n)
  sp_names     <- as.character(meta$SpeciesAccepted)

  for (i in seq_len(n)) {

    sp  <- sp_names[i]
    TU  <- matU_list[[i]]
    TF  <- matF_list[[i]]
    s   <- nrow(TU)

    # -- D: stable stage distribution from A = matU + matF ------------------
    D <- tryCatch(
      .stable_stage_D(TU, TF),
      error = function(e) {
        if (verbose) message("  SKIP [D error] ", sp, ": ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(D)) next

    # -- lambda for reporting ------------------------------------------------
    lambda_val <- tryCatch({
      A  <- TU + TF
      ev <- eigen(A)
      round(max(Re(ev$values)), 4)
    }, error = function(e) NA_real_)

    # -- F_vec: fecundity from row 1 of matF ---------------------------------
    F_vec <- TF[1, ]

    # -- Vk_over_k: Poisson default or user scalar --------------------------
    Vk <- if (is.null(Vk_over_k)) rep(1, s) else rep(as.numeric(Vk_over_k)[1], s)

    # -- Run Ne_sexual_Y2000 -------------------------------------------------
    result <- tryCatch(
      NeStage::Ne_sexual_Y2000(
        T_mat      = TU,
        F_vec      = F_vec,
        D          = D,
        Vk_over_k  = Vk,
        Ne_target  = Ne_target,
        population = sp
      ),
      error = function(e) {
        if (verbose) message("  SKIP [Ne error] ", sp, ": ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(result)) next

    if (verbose)
      message(sprintf("  OK  %-45s  Ne/N = %6.4f  Min N = %d",
                      sp, result$NeN, result$Min_N))

    results_list[[i]] <- result

    summary_rows[[i]] <- data.frame(
      Species         = sp,
      Family          = as.character(meta$Family[i]),
      Order           = as.character(meta$Order[i]),
      Country         = as.character(meta$Country[i]),
      Lat             = if (!is.na(lat_col)) meta[[lat_col]][i] else NA_real_,
      Lon             = if (!is.na(lon_col)) meta[[lon_col]][i] else NA_real_,
      MatrixDimension = s,
      lambda          = lambda_val,
      NeN             = round(result$NeN,    4),
      Min_N           = result$Min_N,
      L               = round(result$L,      3),
      V               = round(result$V,      6),
      u_bar           = round(result$u_bar,  4),
      u2_bar          = round(result$u2_bar, 4),
      Ne_target       = Ne_target,
      stringsAsFactors = FALSE
    )
  }

  # -- Assemble output --------------------------------------------------------
  valid        <- !vapply(results_list, is.null, logical(1))
  results_ok   <- results_list[valid]
  summary_df   <- do.call(rbind, summary_rows[valid])
  rownames(summary_df) <- NULL

  # Name the results list by species
  names(results_ok) <- summary_df$Species

  # Sort by Family then Species
  ord        <- order(summary_df$Family, summary_df$Species)
  summary_df <- summary_df[ord, ]
  results_ok <- results_ok[ord]
  rownames(summary_df) <- NULL

  n_ok   <- sum(valid)
  n_skip <- n - n_ok
  message("Done.  ", n_ok, " species computed, ", n_skip, " skipped.")
  if (n_skip > 0L)
    message("  Re-run with verbose = TRUE to see per-species skip reasons.")

  list(summary = summary_df,
       results = results_ok)
}
