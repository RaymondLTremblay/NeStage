# Overwrite R/nestage-utils.R with a clean, validated version
utils_lines <- c(
  "# =============================================================================",
  "# NeStage utils",
  "# =============================================================================",
  "",
  "# Simple checkers --------------------------------------------------------------",
  "ns_is_square_matrix <- function(A) {",
  "  is.matrix(A) && nrow(A) == ncol(A)",
  "}",
  "",
  "ns_assert_vector <- function(x, len = NULL) {",
  "  stopifnot(is.numeric(x), is.vector(x))",
  "  if (!is.null(len)) stopifnot(length(x) == len)",
  "  invisible(TRUE)",
  "}",
  "",
  "# Build a simple Yonezawa-style A from survival (u) and fecundity (F)",
  "# Row 1 = fecundity only; sub-diagonal = u[1:(k-1)], terminal stasis = u[k]",
  "build_A_yonezawa <- function(u, F) {",
  "  ns_assert_vector(u); ns_assert_vector(F); stopifnot(length(u) == length(F))",
  "  k <- length(u)",
  "  A <- matrix(0, k, k)",
  "  A[1, ] <- F",
  "  if (k >= 2) {",
  "    for (j in 1:(k - 1)) A[j + 1, j] <- u[j]",
  "  }",
  "  A[k, k] <- u[k]",
  "  A",
  "}",
  "",
  "# Assemble Lefkovitch/Leslie A from U and F (row 1 = F)",
  "ns_assemble_A <- function(U, F) {",
  "  stopifnot(ns_is_square_matrix(U), length(F) == nrow(U))",
  "  A <- U",
  "  A[1, ] <- F",
  "  A",
  "}",
  "",
  "# Rounding helpers -------------------------------------------------------------",
  "ns_round_df <- function(df, digits = 3, cols = NULL) {",
  "  if (is.null(cols)) {",
  "    cols <- vapply(df, is.numeric, logical(1))",
  "  } else {",
  "    stopifnot(all(cols %in% names(df)))",
  "    tmp <- logical(ncol(df)); names(tmp) <- names(df)",
  "    tmp[cols] <- TRUE",
  "    cols <- tmp",
  "  }",
  "  out <- df",
  "  for (nm in names(df)[cols]) out[[nm]] <- round(df[[nm]], digits)",
  "  out",
  "}"
)

dir.create("R", showWarnings = FALSE)
con <- file("R/nestage-utils.R", open = "w+", encoding = "UTF-8")
writeLines(utils_lines, con, useBytes = TRUE)
close(con)

# Verify parsing
tryCatch(
  { parse(file = "R/nestage-utils.R"); message("✅ nestage-utils.R parses cleanly.") },
  error = function(e) message("❌ Parse error: ", e$message)
)

# Source it
source("R/nestage-utils.R", echo = TRUE, keep.source = TRUE)
