# Reproduce Yonezawa et al. (2000) Table 4 for the two Fritillaria
# camtschatcensis populations (Miz, Nan) through the current clonal API.
#
# Inputs are Table 2 of the paper (also shown in ?Ne_clonal_Y2000 examples).
# The "observed" column uses the field stage fractions (D_obs) and reproduces
# the paper exactly. The "expected" column uses the stable stage distribution
# derived from the dominant eigenvector of A = T + F_mat (via the _both
# wrapper). Note: the paper's parenthetical column used its own published
# expected fractions, which differ very slightly from the eigen-derived ones
# (e.g. Nan Ny/N: paper 2.446 vs eigen 2.421), so the expected column is
# locked to current behaviour, not the paper's printed parentheticals.

T_Miz <- matrix(c(
  0.789, 0.121, 0.054,
  0.007, 0.621, 0.335,
  0.001, 0.258, 0.611
), nrow = 3, byrow = TRUE)
D_Miz <- c(0.935, 0.038, 0.027)
F_Miz <- c(0.055, 1.328, 2.398)

T_Nan <- matrix(c(
  0.748, 0.137, 0.138,
  0.006, 0.669, 0.374,
  0.001, 0.194, 0.488
), nrow = 3, byrow = TRUE)
D_Nan <- c(0.958, 0.027, 0.015)
F_Nan <- c(0.138, 2.773, 5.016)

test_that("Table 4 observed column is reproduced exactly (Miz, Nan)", {
  miz <- Ne_clonal_Y2000(T_Miz, F_Miz, D_Miz, L = 13.399, population = "Miz")
  expect_equal(round(miz$L, 3), 13.399)
  expect_equal(round(miz$NyN, 3), 2.932)
  expect_equal(round(miz$NeN, 3), 0.219)

  nan <- Ne_clonal_Y2000(T_Nan, F_Nan, D_Nan, L = 8.353, population = "Nan")
  expect_equal(round(nan$L, 3), 8.353)
  expect_equal(round(nan$NyN, 3), 2.428)
  expect_equal(round(nan$NeN, 3), 0.291)
})

test_that("Expected column (eigen-derived D) matches current behaviour", {
  miz <- Ne_clonal_Y2000_both(T_Miz, F_Miz, D_obs = D_Miz, L = 13.399)
  expect_equal(round(miz$observed$NeN, 3), 0.219)
  expect_equal(round(miz$expected$NyN, 3), 2.973)
  expect_equal(round(miz$expected$NeN, 3), 0.222)

  nan <- Ne_clonal_Y2000_both(T_Nan, F_Nan, D_obs = D_Nan, L = 8.353)
  expect_equal(round(nan$observed$NeN, 3), 0.291)
  expect_equal(round(nan$expected$NyN, 3), 2.421)
  expect_equal(round(nan$expected$NeN, 3), 0.290)
})

test_that("Internally computed L is positive and finite", {
  # NOTE: the internal generation-time estimate (eigenvector-weighted T^x
  # method) does NOT equal the paper's reported L (e.g. Miz internal ~4.7 vs
  # paper 13.399). The two use different L definitions, which is why exact
  # Table 4 replication supplies L explicitly. Here we only check the internal
  # estimate is well-formed.
  miz <- Ne_clonal_Y2000(T_Miz, F_Miz, D_Miz, population = "Miz")
  expect_true(is.finite(miz$L) && miz$L > 0)
})
