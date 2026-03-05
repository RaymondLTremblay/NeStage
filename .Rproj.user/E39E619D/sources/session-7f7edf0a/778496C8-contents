# scripts/02_batch_all_populations.R

library(NeStage)

# Helper: load a .RData or .rds file that contains A_obs and optionally A_exp
load_matrix_pair <- function(path) {
  e <- new.env(parent = emptyenv())
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    obj <- readRDS(path)
    if (is.list(obj)) {
      e$A_obs <- obj$A_obs
      e$A_exp <- obj$A_exp
      e$population <- obj$population
    }
  } else if (grepl("\\.RData$|\\.rda$", path, ignore.case = TRUE)) {
    load(path, envir = e)
  } else {
    stop("Unsupported file type: ", path)
  }
  e
}

# Configure
in_dir <- "data/matrices"           # change to your folder with matrices
out_csv <- "outputs/nestage_summary.csv"
dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)

files <- list.files(in_dir, full.names = TRUE)
res_list <- list()

for (fp in files) {
  e <- load_matrix_pair(fp)
  pop <- if (!is.null(e$population)) e$population else tools::file_path_sans_ext(basename(fp))
  df <- nestage_run_all(
    A_obs = e$A_obs,
    A_exp = if (exists("A_exp", envir = e)) e$A_exp else NULL,
    repro_stage = 3,
    frac_repro = 0.10,
    meta = list(population = pop)
  )
  res_list[[pop]] <- df
}

res_all <- do.call(rbind, res_list)
utils::write.csv(res_all, out_csv, row.names = FALSE)
print(res_all)

nestage_log_run(sprintf("Batch run on %d populations from %s", length(res_list), in_dir))
