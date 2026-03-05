# =============================================================================
# NeStage logging (very light)
# =============================================================================

# Toggle with options(nestage.verbose = TRUE/FALSE)
.ns_verbose <- function() {
  opt <- getOption("nestage.verbose")
  if (is.null(opt)) TRUE else isTRUE(opt)
}

ns_log_info <- function(...) {
  if (.ns_verbose()) cat(sprintf("[NeStage][INFO] %s\n", paste0(..., collapse = "")))
}

ns_log_warn <- function(...) {
  cat(sprintf("[NeStage][WARN] %s\n", paste0(..., collapse = "")))
}

ns_log_error <- function(...) {
  stop(sprintf("[NeStage][ERROR] %s", paste0(..., collapse = "")), call. = FALSE)
}
