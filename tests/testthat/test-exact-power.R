# The strongest test in the suite: in the saturated case the engine's entire
# decision rule is a deterministic function of binomial counts, so its power
# has an EXACT finite-sample value by triple binomial enumeration -- no
# approximation of any kind. Simulated power must match it within Monte Carlo
# error. This validates the DGP draw, the ML fit, the APE estimator, the
# delta-method SE, the CI, and the claim logic against arithmetic; the
# rare-outcome case below also certifies the separation convention, which
# the common designs never exercise (their boundary mass is ~1e-80).
# (exact_power_sat lives in helper-exact.R; mirrored in
# validation/run-validation.R, V1.)

test_that("simulated power equals the exact enumerated power (saturated case)", {
  n <- 400
  nsim <- 4000
  d10 <- set_ape(ape_dgp("probit", focal = pa_var("t", "binary", p = 0.5),
                         baseline = 0.30), 0.10)
  d00 <- set_ape(ape_dgp("probit", focal = pa_var("t", "binary", p = 0.5),
                         baseline = 0.30), 0)

  cases <- list(
    list(dgp = d10, claim = "detect", sesoi = NULL, p1 = 0.40, seed = 1),
    list(dgp = d10, claim = "minimum", sesoi = 0.04, p1 = 0.40, seed = 2),
    list(dgp = d00, claim = "equivalence", sesoi = 0.07, p1 = 0.30, seed = 3)
  )
  for (cs in cases) {
    ## the claim's own interval level under the default alpha = .05:
    ## 95% for detection, 90% for the one-sided SESOI claims (1.8.0)
    cl <- powerape:::claim_levels(cs$claim)$conf_claim
    ex <- exact_power_sat(n, 0.5, 0.30, cs$p1, cl, cs$claim,
                          if (is.null(cs$sesoi)) NA else cs$sesoi)
    pw <- ape_power(cs$dgp, n = n, claim = cs$claim, sesoi = cs$sesoi,
                    nsim = nsim, seed = cs$seed)
    mcse <- sqrt(max(ex * (1 - ex), 1e-6) / nsim)
    expect_lt(abs(pw$power - ex), 4 * mcse + 1e-4,
              label = sprintf("|simulated - exact| for claim %s", cs$claim))
  }
})

test_that("separation convention matches exact enumeration where it matters", {
  ## rare outcome, protective effect: .05 -> .01 at n = 300. Over a fifth of
  ## simulated studies have a treated arm without events; the ML effect does
  ## not exist there. "fail" (default) scores them as failed fits, "keep"
  ## with their degenerate Wald intervals (the arm contributes zero
  ## variance); the two conventions differ by ~21 points of exact power.
  ## Negative direction: enumerate with the arms swapped.
  d <- set_ape(ape_dgp("probit", focal = pa_var("t", "binary", p = 0.5),
                       baseline = 0.05), -0.04)
  nsim <- 4000
  for (conv in c("fail", "keep")) {
    ex <- exact_power_sat(300, 0.5, 0.01, 0.05, 0.95, "detect",
                          separation = conv)
    pw <- suppressWarnings(ape_power(d, n = 300, claim = "detect", nsim = nsim,
                                     seed = 11, separation = conv))
    mcse <- sqrt(ex * (1 - ex) / nsim)
    expect_lt(abs(pw$power - ex), 4 * mcse,
              label = sprintf("|simulated - exact| under separation = '%s'", conv))
    expect_gt(pw$separated, 0.15)
  }
  expect_lt(exact_power_sat(300, 0.5, 0.01, 0.05, 0.95, "detect"), 0.36)
  expect_gt(exact_power_sat(300, 0.5, 0.01, 0.05, 0.95, "detect",
                            separation = "keep"), 0.55)
})

test_that("probit and logit give identical saturated-case power on the same seed", {
  ## With no covariates both DGPs produce identical outcome probabilities
  ## (p0 = .30, p1 = .40), so with the same seed the drawn counts and the
  ## saturated MLE fitted probabilities -- hence every CI and decision --
  ## coincide exactly across links.
  dp <- set_ape(ape_dgp("probit", focal = pa_var("t", "binary", p = 0.5),
                        baseline = 0.30), 0.10)
  dl <- set_ape(ape_dgp("logit", focal = pa_var("t", "binary", p = 0.5),
                        baseline = 0.30), 0.10)
  pp <- ape_power(dp, n = 500, claim = "detect", nsim = 1000, seed = 7)
  pl <- ape_power(dl, n = 500, claim = "detect", nsim = 1000, seed = 7)
  expect_equal(pp$power, pl$power, tolerance = 1e-9)
})
