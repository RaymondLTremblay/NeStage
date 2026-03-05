# =============================================================================
# NeStage core
# -----------------------------------------------------------------------------
# v0.6.1 (2026-03-03)
# - Adds estimator dispatch with "paper_v1" (exact Yonezawa Table 4) and "proxy"
# - Keeps public API the same; no breaking changes
# - Fixes lambda assembly: A_full <- U; A_full[1,] <- F  (row 1 = fecundity only)
# =============================================================================

# ---- version & numerics ------------------------------------------------------
ns_version <- function() "v0.6.1"
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
  ee  <- eigen(A)
  idx <- which.max(Re(ee$values))
  val <- Re(ee$values[idx])
  vec <- Re(ee$vectors[, idx])
  if (!is.finite(sum(vec)) || sum(vec) == 0) vec[] <- 1
  w <- vec / sum(vec)
  if (sum(w) < 0) w <- -w
  list(lambda = val, w = w)
}

.ns_dom_left <- function(A) {
  ee  <- eigen(t(A))
  idx <- which.max(Re(ee$values))
  v   <- Re(ee$vectors[, idx])
  if (!is.finite(sum(v))) v[] <- 1
  if (sum(v) < 0) v <- -v
  v
}

.ns_normalize_v <- function(v, w) {
  denom <- sum(v * w)
  if (!is.finite(denom) || abs(denom) < .ns_eps_div) return(v)
  v / denom
}

# ---- proxy estimators (kept for backward compat.) ---------------------------
# Whole-population annual survival: u_bar = sum_j w_j * (sum_{i>=2} A[i, j])
.ns_u_bar <- function(A, w) {
  k <- nrow(A)
  if (k < 2) return(0)
  u_j   <- colSums(A[2:k, , drop = FALSE])
  u_bar <- sum(u_j * w)
  as.numeric(u_bar)
}

# Stage-weighted L proxy: sum_i i * v_i * w_i, v^T w = 1
.ns_generation_time <- function(A) {
  d <- .ns_dom_right(A); w <- d$w
  v <- .ns_normalize_v(.ns_dom_left(A), w)
  denom <- sum(v * w)
  if (!is.finite(denom) || abs(denom) < .ns_eps_div) return(NA_real_)
  idx <- seq_along(w)
  as.numeric(sum(idx * v * w))
}

# TEMPORARY proxy Ny/N and Ne/N (placeholder math)
.ns_estimate_core <- function(A, compute_L = TRUE) {
  d <- .ns_dom_right(A); w <- d$w
  u_bar <- .ns_u_bar(A, w)
  den <- 1 - u_bar
  NyN <- if (!is.finite(den) || abs(den) < .ns_eps_div) NA_real_ else 1 / den
  L   <- if (isTRUE(compute_L)) .ns_generation_time(A) else NA_real_
  NeN <- if (!isTRUE(compute_L) || !is.finite(L) || L <= .ns_eps_div ||
             !is.finite(den) || abs(den) < .ns_eps_div) {
    NA_real_
  } else 1 / (den * L)
  list(NyN = NyN, NeN = NeN, L = L, lambda = d$lambda)
}

# ---- exact Yonezawa (paper_v1) backend -------------------------------------
# Utilities and data for exact replication of Table 4.

# Table-2 inputs (U, F, D_obs, D_exp, L) for Miz / Nan (Yonezawa 2000)
get_table2_inputs <- function(pop = c("Miz","Nan")) {
  pop <- match.arg(pop)
  if (pop == "Miz") {
    U <- matrix(c(
      0.789,0.121,0.054,
      0.007,0.621,0.335,
      0.001,0.258,0.611
    ), nrow=3, byrow=TRUE,
    dimnames=list(c("one_leaf","multileaf_nonfl","multileaf_flower"),
                  c("one_leaf","multileaf_nonfl","multileaf_flower")))
    F <- c(0.055,1.328,2.398)
    D_obs <- c(0.935,0.038,0.027)
    D_exp <- c(0.921,0.046,0.033)
    L <- 13.399
  } else {
    U <- matrix(c(
      0.748,0.137,0.138,
      0.006,0.669,0.374,
      0.001,0.194,0.488
    ), nrow=3, byrow=TRUE,
    dimnames=list(c("one_leaf","multileaf_nonfl","multileaf_flower"),
                  c("one_leaf","multileaf_nonfl","multileaf_flower")))
    F <- c(0.138,2.773,5.016)
    D_obs <- c(0.958,0.027,0.015)
    D_exp <- c(0.951,0.034,0.015)
    L <- 8.353
  }
  list(U=U, F=F, D_obs=D_obs, D_exp=D_exp, L=L,
       stages=c("one_leaf","multileaf_nonfl","multileaf_flower"))
}

# Generation time L by the paper’s age-by-stage definition
.ns_compute_L_paper <- function(U, F, x_max = 500L) {
  stopifnot(nrow(U) == ncol(U), length(F) == nrow(U))
  k <- nrow(U)
  Upow <- diag(k)
  num <- 0; den <- 0
  for (x in 1:x_max) {
    Upow <- Upow %*% U
    u_jx <- Upow[, 1]      # newborn -> stage j at age x
    l_x  <- sum(u_jx)
    if (l_x <= 0) next
    m_x <- sum(F * u_jx) / l_x
    num <- num + x * m_x * l_x
    den <- den +       m_x * l_x
  }
  if (den == 0) return(NA_real_)
  as.numeric(num / den)
}

# Exact Ny/N and Ne/N per Yonezawa Table 4:
# over_u2 = sum_i D_i * u_i^2, with u_i = colSums(U)
# Ny/N    = 1 / (1 - over_u2)
# Ne/N    = 1 / ( L * (1 - over_u2) )
.ns_paper_exact <- function(U, F, D, L) {
  u <- colSums(U)
  over_u2 <- sum(D * (u^2))
  NyN <- 1 / (1 - over_u2)
  NeN <- 1 / (L * (1 - over_u2))

  # Assemble A for lambda reporting (row 1 = fecundity only)
  A_full <- U
  A_full[1, ] <- F
  e <- eigen(A_full)
  lam <- Re(e$values[which.max(Re(e$values))])

  list(NyN = as.numeric(NyN),
       NeN = as.numeric(NeN),
       L   = as.numeric(L),
       lambda = as.numeric(lam))
}

# Dispatcher: estimator = "paper_v1" (exact) or proxy
.ns_estimate_core_paper_v1 <- function(A, meta) {
  # Prefer meta$table2_pop, else use meta$U, meta$F, meta$Ds
  if (!is.null(meta$table2_pop)) {
    t2 <- get_table2_inputs(meta$table2_pop)
    U  <- t2$U; F <- t2$F
    if (is.null(meta$Ds)) meta$Ds <- t2$D_obs
    if (is.null(meta$L))  meta$L  <- t2$L
  } else {
    if (is.null(meta$U) || is.null(meta$F)) {
      if (is.matrix(A)) {
        # Heuristic fallback if only A provided:
        F <- as.numeric(A[1, ])
        U <- A; U[1, ] <- 0
        warning("paper_v1 inferred U,F from A. Provide meta$U and meta$F for exact replication.")
      } else {
        stop("paper_v1 requires meta$table2_pop OR (meta$U, meta$F).")
      }
    } else {
      U <- meta$U; F <- meta$F
    }
  }

  if (is.null(meta$Ds)) stop("paper_v1 requires meta$Ds (stage fractions).")
  D <- as.numeric(meta$Ds); D <- D / sum(D)

  # Optional sexual-only mode (if requested via meta)
  F_use <- F
  if (!is.null(meta$mode) && identical(meta$mode, "sexual_only")) {
    if (!is.null(meta$F_sexual)) {
      if (length(meta$F_sexual) != length(F)) stop("meta$F_sexual length must match F.")
      F_use <- as.numeric(meta$F_sexual)
    } else if (!is.null(meta$sexual_fraction)) {
      sf <- as.numeric(meta$sexual_fraction)
      if (length(sf) == 1L) sf <- rep(sf, length(F))
      if (length(sf) != length(F)) stop("meta$sexual_fraction must be length 1 or length(F).")
      F_use <- sf * F
    } else {
      stop("sexual_only: provide meta$F_sexual or meta$sexual_fraction.")
    }
  }

  L_use <- if (!is.null(meta$L)) as.numeric(meta$L) else .ns_compute_L_paper(U, F_use)
  .ns_paper_exact(U, F_use, D, L_use)
}

# Central dispatcher (proxy by default; exact on request)
.ns_estimate_core_dispatch <- function(A, meta, compute_L = TRUE) {
  est <- if (!is.null(meta$estimator)) meta$estimator else "proxy"
  if (identical(est, "paper_v1")) {
    return(.ns_estimate_core_paper_v1(A, meta))
  } else {
    return(.ns_estimate_core(A, compute_L = TRUE))
  }
}

# ---- sex-only transform (proxy path only) -----------------------------------
.ns_make_sex_only <- function(A_obs, repro_stage, frac_repro) {
  k <- nrow(A_obs)
  A <- matrix(0, k, k)
  repro_stage <- as.integer(repro_stage)
  repro_stage <- repro_stage[repro_stage >= 1 & repro_stage <= k]
  if (length(repro_stage) == 0) stop("Invalid repro_stage indices.")
  for (r in repro_stage) {
    f_obs   <- A_obs[1, r]
    A[1, r] <- max(f_obs * frac_repro, .ns_eps_fecundity)
    A[r, r] <- max(A_obs[r, r], .ns_eps_stasis)   # small stasis for stability
  }
  A
}

# ---- public API --------------------------------------------------------------
nestage_exact <- function(A_obs, meta = list()) {
  if (!is.null(A_obs)) nestage_check_matrix(A_obs)
  est <- .ns_estimate_core_dispatch(A_obs, meta, compute_L = TRUE)
  list(model   = "Exact-Observed",
       NyN     = est$NyN,
       NeN     = est$NeN,
       L       = est$L,
       lambda  = est$lambda,
       version = ns_version(),
       meta    = meta)
}

nestage_expected <- function(A_exp, meta = list()) {
  if (!is.null(A_exp)) nestage_check_matrix(A_exp)
  est <- .ns_estimate_core_dispatch(A_exp, meta, compute_L = TRUE)
  list(model   = "Exact-Expected",
       NyN     = est$NyN,
       NeN     = est$NeN,
       L       = est$L,
       lambda  = est$lambda,
       version = ns_version(),
       meta    = meta)
}

nestage_sex_only <- function(
    A_obs, repro_stage, frac_repro = 0.10, meta = list()
) {
  # Proxy-path transform (works when estimator != "paper_v1")
  nestage_check_matrix(A_obs)
  stopifnot(is.numeric(repro_stage), all(repro_stage %in% seq_len(nrow(A_obs))))
  stopifnot(is.numeric(frac_repro), frac_repro > 0, frac_repro < 1)

  # If user asked for paper_v1 sex-only, pass via meta (ignores A_obs)
  if (!is.null(meta$estimator) && identical(meta$estimator, "paper_v1")) {
    # In paper_v1 mode, use meta$U, meta$F and sexual-only flags
    est <- .ns_estimate_core_dispatch(A_obs, meta, compute_L = TRUE)
    return(list(model   = sprintf("Sex-only (paper_v1)"),
                NyN     = est$NyN,
                NeN     = est$NeN,
                L       = est$L,
                lambda  = est$lambda,
                version = ns_version(),
                meta    = meta))
  }

  # Otherwise proxy-transform the observed A
  A_sex <- .ns_make_sex_only(A_obs, repro_stage, frac_repro)
  est   <- .ns_estimate_core_dispatch(A_sex, meta, compute_L = TRUE)
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
    rbind(
      bind(eo, so, "Sex-only vs Exact-Observed"),
      if (nrow(ex) == 1) bind(eo, ex, "Exact-Expected vs Exact-Observed")
    )
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

# ---- tiny demo matrices (proxy demo) ----------------------------------------
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
