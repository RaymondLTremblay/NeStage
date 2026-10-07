# Tests for the configurable recruitment row (fecundity placement guard).
# The first stage is not always recruitment -- it can be dormancy (a seed
# bank). recruit_row controls which single row of F_mat receives F_vec, and
# the guard ensures fecundity stays confined to exactly that one row.

# Source the current package source files (no installed build required).
for (f in c(
  "R/utils_fecundity.R",
  "R/Ne_clonal_Y2000.R",
  "R/Ne_sexual_Y2000.R",
  "R/Ne_mixed_Y2000.R"
)) {
  if (file.exists(f)) source(f)
}

test_that(".build_F_mat places fecundity in the requested row only", {
  F_vec <- c(0, 0, 2.5, 4.0)
  s <- 4

  F1 <- .build_F_mat(F_vec, s, recruit_row = 1L)
  expect_equal(F1[1, ], F_vec)
  expect_equal(sum(F1[-1, ]), 0) # nothing outside row 1

  F2 <- .build_F_mat(F_vec, s, recruit_row = 2L)
  expect_equal(F2[2, ], F_vec)
  expect_equal(sum(F2[-2, ]), 0) # fecundity confined to row 2 only
  expect_equal(sum(F1[1, ]), sum(F2[2, ])) # same total, different row
})

test_that(".validate_recruit_row rejects invalid rows (the guard)", {
  expect_error(.validate_recruit_row(0L, 4), "valid range")
  expect_error(.validate_recruit_row(5L, 4), "valid range")
  expect_error(.validate_recruit_row(2.5, 4), "whole number")
  expect_error(.validate_recruit_row(c(1L, 2L), 4), "single integer")
  expect_silent(.validate_recruit_row(2L, 4))
})

test_that("Ne_clonal_Y2000 default is unchanged and equals recruit_row = 1", {
  # 3-stage clonal example: stage 1 = dormant seed bank, 2 = juvenile,
  # 3 = adult. Survival/transition matrix (columns = from, rows = to).
  T_mat <- matrix(
    c(
      0.30, 0.00, 0.00, # to dormant
      0.40, 0.50, 0.00, # to juvenile
      0.00, 0.30, 0.80 # to adult
    ),
    nrow = 3, byrow = TRUE
  )
  F_vec <- c(0, 0, 6) # only adults reproduce
  D <- c(0.5, 0.3, 0.2)

  base <- Ne_clonal_Y2000(T_mat, F_vec, D, population = "test")
  r1 <- Ne_clonal_Y2000(T_mat, F_vec, D, population = "test", recruit_row = 1L)
  expect_equal(base$NeN, r1$NeN)
  expect_equal(base$L, r1$L)

  # Sending recruits to the dormant bank (row 1) vs the juvenile stage (row 2)
  # changes the stable stage distribution and hence generation time L.
  r2 <- Ne_clonal_Y2000(T_mat, F_vec, D, population = "test", recruit_row = 2L)
  expect_false(isTRUE(all.equal(r1$L, r2$L)))

  # Out-of-range recruit_row is caught before any computation.
  expect_error(
    Ne_clonal_Y2000(T_mat, F_vec, D, recruit_row = 4L),
    "valid range"
  )
})
