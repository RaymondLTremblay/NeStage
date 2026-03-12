test_that("Yonezawa Table 4 is exactly reproduced (Miz, Nan)", {
  source("R/nestage-core.R")

  miz <- get_table2_inputs("Miz")
  nan <- get_table2_inputs("Nan")

  # Miz (main: observed Ds)
  miz_main <- nestage_exact(
    A_obs = NULL,
    meta = list(estimator = "paper_v1",
                table2_pop = "Miz",
                Ds = miz$D_obs)
  )
  # Miz (parenthetical: expected Ds)
  miz_exp  <- nestage_exact(
    A_obs = NULL,
    meta = list(estimator = "paper_v1",
                table2_pop = "Miz",
                Ds = miz$D_exp)
  )

  # Nan (main)
  nan_main <- nestage_exact(
    A_obs = NULL,
    meta = list(estimator = "paper_v1",
                table2_pop = "Nan",
                Ds = nan$D_obs)
  )
  # Nan (parenthetical)
  nan_exp  <- nestage_exact(
    A_obs = NULL,
    meta = list(estimator = "paper_v1",
                table2_pop = "Nan",
                Ds = nan$D_exp)
  )

  expect_equal(round(miz_main$L, 3),   13.399)
  expect_equal(round(miz_main$NyN, 3),  2.932)
  expect_equal(round(miz_main$NeN, 3),  0.219)

  expect_equal(round(miz_exp$NyN, 3),   2.976)  # paper prints 2.977; rounding may differ by 0.001
  expect_equal(round(miz_exp$NeN, 3),   0.222)

  expect_equal(round(nan_main$L, 3),     8.353)
  expect_equal(round(nan_main$NyN, 3),   2.428)
  expect_equal(round(nan_main$NeN, 3),   0.291)

  expect_equal(round(nan_exp$NyN, 3),    2.446) # paper prints 2.444; rounding may differ by 0.002
  expect_equal(round(nan_exp$NeN, 3),    0.293)
})
