# The equivalence claim needs BOTH interval bounds inside the margin, so the
# searches of ape_n() and ape_mde() must invert the two-bound (TOST) power,
# not the one-sided map z + qnorm(p) the detect and minimum claims use.
# Through 1.11.1 they inverted the one-sided map: at a truth of 0 the first
# proposal sat at about 72% of the requirement, the Newton steps oscillated,
# and the exact power at the answer fell more than the confirmation
# tolerance below the goal in 3-5% of runs (validation-gaps review,
# 2026-10-14, item A1). Exact oracle: exact_power_sat() (helper-exact.R).

test_that("req_delta inverts the one-sided and the two-bound power maps", {
  rd <- powerape:::req_delta
  z <- qnorm(0.95)
  for (p in c(0.5, 0.8, 0.9, 0.95)) {
    expect_identical(rd(p, z), z + qnorm(p))                 # one-sided claims
    expect_equal(rd(p, z, 1), z + qnorm((1 + p) / 2), tolerance = 1e-12)
    for (r in c(1.2, 2, 5, 50)) {
      d <- rd(p, z, r)
      expect_equal(pnorm(d - z) + pnorm(r * d - z) - 1, p, tolerance = 1e-9)
    }
  }
  ## the far bound matters less as it recedes: decreasing in r, one-sided in the limit
  ds <- vapply(c(1, 1.5, 3, 10), function(r) rd(0.8, z, r), numeric(1))
  expect_true(all(diff(ds) < 0))
  expect_equal(rd(0.8, z, 1e6), z + qnorm(0.8), tolerance = 1e-9)
})

test_that("req_margin is the TOST inverse in the margin", {
  rm <- powerape:::req_margin
  z <- qnorm(0.95)
  s <- 0.02
  for (t in c(0, 0.01, 0.03)) {
    m <- rm(0.8, s, t, z)
    expect_equal(pnorm((m - t) / s - z) + pnorm((m + t) / s - z) - 1, 0.8,
                 tolerance = 1e-9)
  }
  expect_equal(rm(0.8, s, 0, z), s * (z + qnorm(0.9)), tolerance = 1e-12)
})

test_that("ape_n's equivalence search starts at the TOST requirement and lands on it", {
  ## saturated design, truth 0, SESOI .10: the TOST normal formula needs
  ## n = .84 ((z + qnorm(.9)) / .10)^2 = 719.4; the one-sided map proposed 519
  d0 <- set_ape(ape_dgp("probit", focal = pa_var("t", "binary", p = 0.5),
                        baseline = 0.30), 0)
  an <- ape_n(d0, power = 0.80, claim = "equivalence", sesoi = 0.10,
              nsim = 500, seed = 21)
  n_tost <- 0.84 * ((qnorm(0.95) + qnorm(0.90)) / 0.10)^2
  first <- an$history$n[an$history$stage == "search"][1L]
  expect_lt(abs(first / n_tost - 1), 0.04)
  expect_lte(sum(an$history$stage == "search"), 2L)
  ex <- exact_power_sat(an$n, 0.5, 0.30, 0.30, 0.90, "equivalence", sesoi = 0.10)
  expect_gt(ex, 0.80 - 0.015)
  expect_lt(ex, 0.83)

  ## a truth of .02 against the same SESOI: the far bound sits 1.5 times as
  ## far as the near one, so both still matter
  d2 <- set_ape(ape_dgp("probit", focal = pa_var("t", "binary", p = 0.5),
                        baseline = 0.30), 0.02)
  a2 <- ape_n(d2, power = 0.80, claim = "equivalence", sesoi = 0.10,
              nsim = 500, seed = 22)
  c2 <- sqrt((0.30 * 0.70 + 0.32 * 0.68) / 0.5)
  n2 <- (c2 * powerape:::req_delta(0.80, qnorm(0.95), 0.12 / 0.08) / 0.08)^2
  first2 <- a2$history$n[a2$history$stage == "search"][1L]
  expect_lt(abs(first2 / n2 - 1), 0.04)
  ex2 <- exact_power_sat(a2$n, 0.5, 0.30, 0.32, 0.90, "equivalence", sesoi = 0.10)
  expect_gt(ex2, 0.80 - 0.015)
  expect_lt(ex2, 0.83)
})

test_that("ape_mde's equivalence margin search lands on the exact margin from both sides", {
  d0 <- set_ape(ape_dgp("probit", focal = pa_var("t", "binary", p = 0.5),
                        baseline = 0.30), 0)
  m <- ape_mde(d0, n = 720, claim = "equivalence", nsim = 500, seed = 23)
  expect_true(m$is_margin)
  first <- m$history$value[m$history$stage == "search"][1L]
  m_tost <- sqrt(0.84 / 720) * (qnorm(0.95) + qnorm(0.90))     # .1000
  expect_lt(abs(first / m_tost - 1), 0.03)
  expect_lte(sum(m$history$stage == "search"), 2L)
  ex <- exact_power_sat(720, 0.5, 0.30, 0.30, 0.90, "equivalence", sesoi = m$mde)
  expect_gt(ex, 0.80 - 0.015)
  expect_lt(ex, 0.83)
})
