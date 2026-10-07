# -----------------------------------------------------------------------------
# Internal helpers: placing the fecundity vector into the projection matrix
# -----------------------------------------------------------------------------
# Every NeStage model builds the full projection matrix as A = T_mat + F_mat,
# where T_mat is survival/transition and F_mat carries reproduction. The
# fecundity vector F_vec gives newborns produced *by* each stage (the columns).
# Those newborns all enter a single recruitment stage -- the row of F_mat they
# land in.
#
# Historically that row was hard-wired to 1 ("F_mat[1, ] <- F_vec"), which
# assumes the first stage is the recruitment stage. That is wrong when the
# first stage is dormancy (e.g. a seed bank): recruits actually enter a later,
# active stage. `recruit_row` makes the destination explicit, while the
# validator guarantees fecundity stays confined to exactly one known row.

.validate_recruit_row <- function(recruit_row, s) {
  # recruit_row: the single stage (matrix row) that newborns enter.
  # Must name exactly one valid stage so fecundity occupies one row only.
  if (!is.numeric(recruit_row) || length(recruit_row) != 1L) {
    stop("recruit_row must be a single integer (the stage that newborns enter).")
  }
  if (!is.finite(recruit_row) || recruit_row != round(recruit_row)) {
    stop("recruit_row must be a whole number.")
  }
  recruit_row <- as.integer(recruit_row)
  if (recruit_row < 1L || recruit_row > s) {
    stop(sprintf(
      "recruit_row = %d is outside the valid range 1:%d (number of stages).",
      recruit_row,
      s
    ))
  }
  invisible(TRUE)
}

.build_F_mat <- function(F_vec, s, recruit_row = 1L) {
  # Build the s x s fecundity matrix from the length-s fecundity vector.
  # Newborns produced by each stage (columns) all enter stage `recruit_row`
  # (the row). By construction fecundity is confined to that one row -- the
  # invariant the rest of the package relies on when splitting A = T + F.
  .validate_recruit_row(recruit_row, s)
  recruit_row <- as.integer(recruit_row)
  F_mat <- matrix(0, s, s)
  F_mat[recruit_row, ] <- F_vec
  F_mat
}
