# Seeded calls run under with_seed() (utils.R), which must leave the
# caller's random-number stream exactly where it was.

test_that("a seeded ape_power() call does not disturb the caller's RNG stream", {
  d <- set_ape(ape_dgp("probit", focal = pa_var("t", "binary", p = 0.5),
                       baseline = 0.30), 0.10)
  set.seed(123)
  a <- runif(3)
  set.seed(123)
  invisible(ape_power(d, 200, claim = "detect", nsim = 20, seed = 9))
  b <- runif(3)
  expect_identical(a, b)
})
