# The strongest test in the suite: in the saturated case the engine's entire
# decision rule is a deterministic function of binomial counts, so its power
# has an EXACT finite-sample value by triple binomial enumeration -- no
# approximation of any kind. Simulated power must match it within Monte Carlo
# error. This validates the DGP draw, the ML fit, the APE estimator, the
# delta-method SE, the CI, the claim logic, and the failed-fit handling, all
# at once, against arithmetic. (exact_power_sat lives in helper-exact.R;
# mirrored in validation/run-validation.R, V1.)

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
    ex <- exact_power_sat(n, 0.5, 0.30, cs$p1, 0.95, cs$claim,
                          if (is.null(cs$sesoi)) NA else cs$sesoi)
    pw <- ape_power(cs$dgp, n = n, claim = cs$claim, sesoi = cs$sesoi,
                    nsim = nsim, seed = cs$seed)
    mcse <- sqrt(max(ex * (1 - ex), 1e-6) / nsim)
    expect_lt(abs(pw$power - ex), 4 * mcse + 1e-4,
              label = sprintf("|simulated - exact| for claim %s", cs$claim))
  }
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
