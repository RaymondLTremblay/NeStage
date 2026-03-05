## scripts/03_yonezawa_demo.R

yonezawa_example_matrices <- function(repro_stage = 3) {
  stopifnot(repro_stage %in% 1:3)

  # Observed
  A_obs <- matrix(0, 3, 3)
  A_obs[1, 3] <- 0.9
  A_obs[1, 1] <- 0.10; A_obs[2, 2] <- 0.20; A_obs[3, 3] <- 0.60
  A_obs[2, 1] <- 0.35; A_obs[3, 2] <- 0.25
  A_obs[1, 2] <- 0.05; A_obs[2, 3] <- 0.10

  # Expected (perturbed)
  A_exp <- matrix(0, 3, 3)
  A_exp[1, 3] <- 1.0
  A_exp[1, 1] <- 0.12; A_exp[2, 2] <- 0.22; A_exp[3, 3] <- 0.58
  A_exp[2, 1] <- 0.30; A_exp[3, 2] <- 0.28
  A_exp[1, 2] <- 0.06; A_exp[2, 3] <- 0.09

  list(A_obs = A_obs, A_exp = A_exp)
}

run_yonezawa_demo <- function(save_plots = TRUE) {
  ex <- yonezawa_example_matrices(3)

  res <- nestage_run_all(
    A_obs = ex$A_obs,
    A_exp = ex$A_exp,
    repro_stage = 3,
    frac_repro = 0.10,
    meta = list(population = "Yonezawa-demo")
  )
  print(res)

  cmp <- nestage_compare(res)
  print(cmp)

  # Map demo label to the paper naming (choose Miz or Nan for the demo)
  pop_map <- c("Yonezawa-demo" = "Miz", "Miz" = "Miz", "Nan" = "Nan")
  res_named <- nestage_rename_populations(res, pop_map, order_levels = TRUE)

  p1 <- nestage_plot_ratios_gg(res_named)
  p2 <- nestage_plot_generation_time_gg(res_named)

  if (isTRUE(save_plots)) {
    ggsave("outputs/ratios_yonezawa_demo.png", p1, width = 9, height = 4.5, dpi = 300)
    ggsave("outputs/L_yonezawa_demo.png", p2, width = 6.5, height = 4.2, dpi = 300)
  }

  invisible(list(res = res, cmp = cmp, p1 = p1, p2 = p2))
}
