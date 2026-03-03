# scripts/01_run_single_population.R

# Load NeStage (already installed from your local source)
library(NeStage)

# TODO: Replace with your real observed/expected matrices
# A_Miz_obs <- ...
# A_Miz_exp <- ...
# A_Nan_obs <- ...
# A_Nan_exp <- ...

# Example: run Miz
miz <- nestage_run_all(
  A_obs = A_Miz_obs,
  A_exp = A_Miz_exp,
  repro_stage = 3,
  frac_repro = 0.10,
  meta = list(population = "Miz")
)

# Example: run Nan
nan <- nestage_run_all(
  A_obs = A_Nan_obs,
  A_exp = A_Nan_exp,
  repro_stage = 3,
  frac_repro = 0.10,
  meta = list(population = "Nan")
)

res_all <- rbind(miz, nan)
print(res_all)

# Compare models
cmp <- nestage_compare(res_all)
print(cmp)

# Plot quick looks
nestage_plot_ratios(res_all)
nestage_plot_generation_time(res_all)

# Log the run
nestage_log_run("Ran Miz and Nan with frac_repro=0.10, repro_stage=3.")
