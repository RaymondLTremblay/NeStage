# =============================================================================
# NeStage examples (runnable)
# =============================================================================

# --- Example 1: Reproduce Yonezawa Table 4 (Miz / Nan) ----------------------
nestage_example_table4 <- function(verbose = TRUE) {
  # Miz (main: observed Ds)
  miz <- get_table2_inputs("Miz")
  out_miz_main <- nestage_exact(
    A_obs = NULL,
    meta = list(estimator = "paper_v1",
                table2_pop = "Miz",
                Ds = miz$D_obs,
                population = "Miz")
  )
  # Miz (parenthetical: expected Ds)
  out_miz_exp <- nestage_exact(
    A_obs = NULL,
    meta = list(estimator = "paper_v1",
                table2_pop = "Miz",
                Ds = miz$D_exp,
                population = "Miz (expected)")
  )

  # Nan
  nan <- get_table2_inputs("Nan")
  out_nan_main <- nestage_exact(
    A_obs = NULL,
    meta = list(estimator = "paper_v1",
                table2_pop = "Nan",
                Ds = nan$D_obs,
                population = "Nan")
  )
  out_nan_exp <- nestage_exact(
    A_obs = NULL,
    meta = list(estimator = "paper_v1",
                table2_pop = "Nan",
                Ds = nan$D_exp,
                population = "Nan (expected)")
  )

  if (isTRUE(verbose)) {
    cat("\n=== Table 4 replication ===\n")
    cat("Miz: L=", round(out_miz_main$L,3),
        " Ny/N=", round(out_miz_main$NyN,3), " (", round(out_miz_exp$NyN,3), ")",
        " Ne/N=", round(out_miz_main$NeN,3), " (", round(out_miz_exp$NeN,3), ")",
        "\n", sep = "")
    cat("Nan: L=", round(out_nan_main$L,3),
        " Ny/N=", round(out_nan_main$NyN,3), " (", round(out_nan_exp$NyN,3), ")",
        " Ne/N=", round(out_nan_main$NeN,3), " (", round(out_nan_exp$NeN,3), ")",
        "\n", sep = "")
  }

  invisible(list(miz_main = out_miz_main, miz_exp = out_miz_exp,
                 nan_main = out_nan_main, nan_exp = out_nan_exp))
}

# --- Example 2: Proxy demo with toy matrices ---------------------------------
nestage_example_proxy_demo <- function() {
  ex <- yonezawa_example_matrices(repro_stage = 3)
  td <- nestage_run_all(
    A_obs       = ex$A_obs,
    A_exp       = ex$A_exp,
    repro_stage = ex$repro_stage,
    frac_repro  = 0.10,
    meta        = list(population = "demo")
  )
  print(td)
  ns_assert_finite(td)
  cmp <- nestage_compare(td); print(cmp)
  invisible(list(td = td, cmp = cmp))
}

# --- Example 3: Paper_v1 sexual-only (U,F with sexual_fraction) --------------
nestage_example_sex_only_paper <- function(pop = c("Miz","Nan")) {
  pop <- match.arg(pop)
  t2 <- get_table2_inputs(pop)

  # Example: only 10% of stage-3 reproduction is sexual
  out_sex <- nestage_exact(
    A_obs = NULL,
    meta = list(estimator = "paper_v1",
                table2_pop = pop,
                Ds = t2$D_obs,
                mode = "sexual_only",
                sexual_fraction = c(0, 0, 0.10),
                population = paste0(pop, " (sex-only 10% stage3)"))
  )
  print(out_sex)
  invisible(out_sex)
}
