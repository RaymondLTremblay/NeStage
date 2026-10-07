# Generation time L must reproduce Table 4 of Yonezawa et al. (2000)
# when computed internally (L = NULL), not only when supplied by the user.

T_Miz <- matrix(c(
  0.789, 0.121, 0.054,
  0.007, 0.621, 0.335,
  0.001, 0.258, 0.611
), nrow = 3, byrow = TRUE)
F_Miz <- c(0.055, 1.328, 2.398)

T_Nan <- matrix(c(
  0.748, 0.137, 0.138,
  0.006, 0.669, 0.374,
  0.001, 0.194, 0.488
), nrow = 3, byrow = TRUE)
F_Nan <- c(0.138, 2.773, 5.016)

test_that("internal L matches Yonezawa (2000) Table 4 for all three models", {
  for (fn in list(.compute_L_clonal, .compute_L_sexual, .compute_L_mixed)) {
    expect_equal(round(fn(T_Miz, F_Miz), 3), 13.399)
    expect_equal(round(fn(T_Nan, F_Nan), 3), 8.353)
  }
})

test_that("Ne_clonal_Y2000 computes L = 13.399 for Miz when L is not supplied", {
  out <- Ne_clonal_Y2000(
    T_mat = T_Miz,
    F_vec = F_Miz,
    D     = c(0.935, 0.038, 0.027)
  )
  expect_equal(out$L_source, "computed")
  expect_equal(round(out$L, 3), 13.399)
})

test_that("L counts survival once (single-stage closed form)", {
  # One stage, annual survival p, fecundity f: l_x * m_x = f * p^x, so
  # L = sum(x p^x) / sum(p^x) = 1 / (1 - p) for x = 1, 2, ...
  p <- 0.8
  L <- .compute_L_clonal(matrix(p), 1, x_max = 2000L)
  expect_equal(L, 1 / (1 - p), tolerance = 1e-8)
})
