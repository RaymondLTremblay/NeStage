# =============================================================================
# NeStage core (scripts mode)
# =============================================================================

# ---- version & numerics ------------------------------------------------------
ns_version <- function() "v0.6.0"
.ns_eps_fecundity <- 1e-10
.ns_eps_stasis    <- 1e-6
.ns_eps_div       <- 1e-12

# ---- validation --------------------------------------------------------------
nestage_check_matrix <- function(A) {
  stopifnot(is.matrix(A), is.numeric(A), nrow(A) == ncol(A))
  if (any(!is.finite(A))) stop("Matrix contains non-finite values.")
  invisible(TRUE)
}

# ---- eigen helpers -----------------------------------------------------------
.ns_dom_right <- function(A) {
  ee <- eigen(A)
  idx <- which.max(Re(ee$values))
  val <- Re(ee$values[idx])
  vec <- Re(ee$vectors[, idx])
  if (!is.finite(sum(vec)) || sum(vec) == 0) vec[] <- 1
  w <- vec / sum(vec)
  if (sum(w) < 0) w <- -w
  list(lambda = val, w = w)
}

.ns_dom_left <- function(A) {
  ee <- eigen(t(A))
  idx <- which.max(Re(ee$values))
  v  <- Re(ee$vectors[, idx])
  if (!is.finite(sum(v))) v[] <- 1
  if (sum(v) < 0) v <- -v
  v
}

.ns_normalize_v <- function(v, w) {
  denom <- sum(v * w)
  if (!is.finite(denom) || abs(denom) < .ns_eps_div) return(v)
  v / denom
}

# ---- pieces for proxy estimators --------------------------------------------
# Whole-population annual survival: u_bar = sum_j w_j * (sum_{i>=2} A[i, j])
.ns_u_bar <- function(A, w) {
  k <- nrow(A)
  if (k < 2) return(0)
  u_j   <- colSums(A[2:k, , drop = FALSE])
  u_bar <- sum(u_j * w)
  as.numeric(u_bar)
}

# Stage-weighted generation time proxy: L = sum_i i * v_i * w_i, with v^T w = 1
.ns_generation_time <- function(A) {
  d <- .ns_dom_right(A); w <- d$w
  v <- .ns_normalize_v(.ns_dom_left(A), w)
  denom <- sum(v * w)
  if (!is.finite(denom) || abs(denom) < .ns_eps_div) return(NA_real_)
  idx <- seq_along(w)
  as.numeric(sum(idx * v * w))
}

# ---------------------- TEMPORARY PROXY ESTIMATORS ---------------------------
# Replace this block with your v0.6.0 math when ready.
.ns_estimate_core <- function(A, compute_L = TRUE) {
  d <- .ns_dom_right(A); w <- d$w
  u_bar <- .ns_u_bar(A, w)

  # Ny/N (proxy)
  den <- 1 - u_bar
  NyN <- if (!is.finite(den) || abs(den) < .ns_eps_div) NA_real_ else 1 / den

  # L (proxy)
  L <- if (isTRUE(compute_L)) .ns_generation_time(A) else NA_real_

  # Ne/N (proxy)
  NeN <- if (!isTRUE(compute_L) || !is.finite(L) || L <= .ns_eps_div ||
             !is.finite(den) || abs(den) < .ns_eps_div) {
    NA_real_
  } else {
    1 / (den * L)
  }

  list(NyN = NyN, NeN = NeN, L = L, lambda = d$lambda)
}
# -------------------- END TEMPORARY PROXY ESTIMATORS -------------------------

# >>>>>>>>>>>>>>>>>>>>>>>>  PASTE v0.6.0 HERE  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# When you are ready, replace `.ns_estimate_core()` above with your exact
# v0.6.0 routines that compute NyN, NeN and (when requested) L from A.
# Keep the function signature and returned list the same.

# ---- sex-only transform (stable) --------------------------------------------
.ns_make_sex_only <- function(A_obs, repro_stage, frac_repro) {
  k <- nrow(A_obs)
  A <- matrix(0, k, k)
  repro_stage <- as.integer(repro_stage)
  repro_stage <- repro_stage[repro_stage >= 1 & repro_stage <= k]
  if (length(repro_stage) == 0) stop("Invalid repro_stage indices.")

  for (r in repro_stage) {
    f_obs <- A_obs[1, r]
    A[1, r] <- max(f_obs * frac_repro, .ns_eps_fecundity)
    A[r, r] <- max(A_obs[r, r], .ns_eps_stasis)   # retain small stasis
  }
  A
}

# ---- public API --------------------------------------------------------------
nestage_exact <- function(A_obs, meta = list()) {
  nestage_check_matrix(A_obs)
  est <- .ns_estimate_core(A_obs, compute_L = TRUE)
  list(model   = "Exact-Observed",
       NyN     = est$NyN,
       NeN     = est$NeN,
       L       = est$L,
       lambda  = est$lambda,
       version = ns_version(),
       meta    = meta)
}

nestage_expected <- function(A_exp, meta = list()) {
  nestage_check_matrix(A_exp)
  est <- .ns_estimate_core(A_exp, compute_L = TRUE)
  list(model   = "Exact-Expected",
       NyN     = est$NyN,
       NeN     = est$NeN,
       L       = est$L,
       lambda  = est$lambda,
       version = ns_version(),
       meta    = meta)
}

nestage_sex_only <- function(A_obs, repro_stage, frac_repro = 0.10, meta = list()) {
  nestage_check_matrix(A_obs)
  stopifnot(is.numeric(repro_stage), all(repro_stage %in% seq_len(nrow(A_obs))))
  stopifnot(is.numeric(frac_repro), frac_repro > 0, frac_repro < 1)

  A_sex <- .ns_make_sex_only(A_obs, repro_stage, frac_repro)
  est   <- .ns_estimate_core(A_sex, compute_L = TRUE)

  list(model   = sprintf("Sex-only (%.0f%% stage)", 100 * frac_repro),
       NyN     = est$NyN,
       NeN     = est$NeN,
       L       = est$L,
       lambda  = est$lambda,
       version = ns_version(),
       meta    = c(meta, list(repro_stage = repro_stage, frac_repro = frac_repro)))
}

nestage_run_all <- function(A_obs, A_exp = NULL, repro_stage, frac_repro = 0.10, meta = list()) {
  out <- list(
    nestage_exact(A_obs, meta),
    if (!is.null(A_exp)) nestage_expected(A_exp, meta) else NULL,
    nestage_sex_only(A_obs, repro_stage, frac_repro, meta)
  )
  out <- out[!vapply(out, is.null, logical(1))]

  td <- do.call(rbind, lapply(out, function(x)
    data.frame(population = if (!is.null(x$meta$population)) x$meta$population else NA_character_,
               model      = x$model,
               NyN        = x$NyN,
               NeN        = x$NeN,
               L          = x$L,
               lambda     = x$lambda,
               version    = x$version,
               stringsAsFactors = FALSE)))
  rownames(td) <- NULL
  td
}

nestage_compare <- function(df) {
  by_pop <- split(df, df$population)
  out <- lapply(by_pop, function(dd) {
    eo <- dd[dd$model == "Exact-Observed", , drop = FALSE]
    so <- dd[grepl("^Sex-only", dd$model), , drop = FALSE]
    ex <- dd[dd$model == "Exact-Expected", , drop = FALSE]
    bind <- function(a, b, tag) {
      if (nrow(a) == 1 && nrow(b) == 1) {
        data.frame(
          population = a$population,
          compare    = tag,
          d_NyN      = b$NyN - a$NyN,
          fc_NyN     = b$NyN / a$NyN,
          d_NeN      = b$NeN - a$NeN,
          fc_NeN     = b$NeN / a$NeN,
          d_L        = b$L - a$L,
          stringsAsFactors = FALSE
        )
      } else NULL
    }
    rbind(bind(eo, so, "Sex-only vs Exact-Observed"),
          if (nrow(ex) == 1) bind(eo, ex, "Exact-Expected vs Exact-Observed"))
  })
  do.call(rbind, out)
}

# ---- diagnostics -------------------------------------------------------------
ns_assert_finite <- function(df) {
  if (!any(is.finite(df$NyN))) stop("All NyN are NA/Non-finite (check estimators).")
  if (!any(is.finite(df$NeN))) stop("All NeN are NA/Non-finite (check estimators).")
  if (any(grepl("^Sex-only", df$model))) {
    dL <- df[grepl("^Sex-only", df$model), , drop = FALSE]
    if (!any(is.finite(dL$L))) stop("All L in Sex-only rows are NA/Non-finite (check generation-time).")
  }
  invisible(TRUE)
}

# ---- tiny demo matrices (Yonezawa-style) ------------------------------------
yonezawa_example_matrices <- function(repro_stage = 3) {
  stopifnot(repro_stage %in% 1:3)
  A_obs <- matrix(0, 3, 3)
  A_obs[1, 3] <- 0.9
  A_obs[1, 1] <- 0.10; A_obs[2, 2] <- 0.20; A_obs[3, 3] <- 0.60
  A_obs[2, 1] <- 0.35; A_obs[3, 2] <- 0.25
  A_obs[1, 2] <- 0.05; A_obs[2, 3] <- 0.10

  A_exp <- matrix(0, 3, 3)
  A_exp[1, 3] <- 1.0
  A_exp[1, 1] <- 0.12; A_exp[2, 2] <- 0.22; A_exp[3, 3] <- 0.58
  A_exp[2, 1] <- 0.30; A_exp[3, 2] <- 0.28
  A_exp[1, 2] <- 0.06; A_exp[2, 3] <- 0.09

  list(A_obs = A_obs, A_exp = A_exp, repro_stage = repro_stage)
}
