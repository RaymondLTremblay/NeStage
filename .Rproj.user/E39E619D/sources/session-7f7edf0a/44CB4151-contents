test_that("Examples run without error", {
  source("R/nestage-core.R")
  source("R/nestage-examples.R")

  expect_silent({
    out <- nestage_example_table4(FALSE)
    stopifnot(is.list(out), length(out) == 4)
  })

  expect_silent({
    ex <- nestage_example_proxy_demo()
    stopifnot(is.list(ex))
  })
})
