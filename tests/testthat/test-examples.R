# Smoke tests: the documented examples for each exported model run without
# error and return finite Ne/N, under both the default recruitment row and a
# shifted recruitment row (dormancy case).

T_plant <- matrix(c(
  0.30, 0.05, 0.00,
  0.40, 0.65, 0.10,
  0.00, 0.20, 0.80
), nrow = 3, byrow = TRUE)
F_plant <- c(0.0, 0.5, 3.0)
D_plant <- c(0.60, 0.25, 0.15)

test_that("Clonal model and its _both wrapper run and return finite Ne/N", {
  r <- Ne_clonal_Y2000(T_plant, F_plant, D_plant, population = "smoke")
  expect_true(is.finite(r$NeN) && r$NeN > 0)
  expect_true(is.finite(r$L) && r$L > 0)

  b <- Ne_clonal_Y2000_both(T_plant, F_plant, D_obs = D_plant)
  expect_true(is.finite(b$observed$NeN))
  expect_true(is.finite(b$expected$NeN))
})

test_that("Sexual model runs, and higher repro variance lowers Ne/N", {
  base <- Ne_sexual_Y2000(T_plant, F_plant, D_plant, population = "smoke")
  hivar <- Ne_sexual_Y2000(
    T_plant, F_plant, D_plant,
    Vk_over_k = c(1, 1, 3), population = "smoke"
  )
  expect_true(is.finite(base$NeN))
  expect_lt(hivar$NeN, base$NeN)
})

test_that("Mixed model runs with sexual + clonal reproduction", {
  r <- Ne_mixed_Y2000(
    T_plant, F_plant, D_plant,
    d = c(0.0, 0.0, 0.7), population = "smoke"
  )
  expect_true(is.finite(r$NeN) && r$NeN > 0)
})

test_that("Ne_sensitivity_L sweep runs on the clonal model", {
  out <- Ne_sensitivity_L(
    model_fn = Ne_clonal_Y2000,
    T_mat = T_plant,
    F_vec = F_plant,
    D = D_plant,
    L_range = seq(4, 12, by = 2)
  )
  expect_s3_class(out$data, "data.frame")
  expect_true(all(is.finite(out$data$NeN)))
})

test_that("recruit_row flows end-to-end through every model", {
  # Default (row 1) vs recruits into stage 2: results should differ but stay
  # finite. Out-of-range rows are rejected before computation.
  for (fn_call in list(
    function(rr) Ne_clonal_Y2000(T_plant, F_plant, D_plant, recruit_row = rr),
    function(rr) Ne_sexual_Y2000(T_plant, F_plant, D_plant, recruit_row = rr),
    function(rr) {
      Ne_mixed_Y2000(
        T_plant, F_plant, D_plant,
        d = c(0, 0, 0.7), recruit_row = rr
      )
    }
  )) {
    r1 <- fn_call(1L)
    r2 <- fn_call(2L)
    expect_true(is.finite(r1$NeN) && is.finite(r2$NeN))
    expect_false(isTRUE(all.equal(r1$L, r2$L)))
    expect_error(fn_call(4L), "valid range")
  }
})
