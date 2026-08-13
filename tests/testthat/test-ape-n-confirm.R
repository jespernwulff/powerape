test_that("ape_n confirmation tightens the answer and reports itself", {
  d <- set_ape(ape_dgp("probit", focal = pa_var("t", "binary", p = 0.5),
                       baseline = 0.30), 0.10)
  an <- ape_n(d, power = 0.80, claim = "detect", nsim = 1000, seed = 11)
  expect_true(an$confirmed)
  expect_true(any(an$history$stage == "confirm"))
  expect_identical(an$nsim_confirm, 4000L)
  ## the returned power/MCSE come from the high-precision confirmation run
  expect_lt(an$mcse, 0.009)
  ## exact finite-sample power at the returned n is close to the goal
  ex <- exact_power_sat(an$n, 0.5, 0.30, 0.40, 0.95, "detect")
  expect_gt(ex, 0.80 - 0.015)
  expect_output(print(an), "confirmed")
})

test_that("confirm = FALSE restores the fast single-stage behaviour", {
  d <- set_ape(ape_dgp("probit", focal = pa_var("t", "binary", p = 0.5),
                       baseline = 0.30), 0.10)
  an <- ape_n(d, power = 0.80, claim = "detect", nsim = 800, seed = 12,
              confirm = FALSE)
  expect_false(an$confirmed)
  expect_false(any(an$history$stage == "confirm"))
  expect_true(is.na(an$nsim_confirm))
  expect_output(print(an), "larger nsim")
})
