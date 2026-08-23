# Claim-semantics concordance with TOSTER -- a DEFAULT-CONVENTION audit.
# powerape's default alpha = .05 must match TOSTER at ITS default alpha with
# no remapping: the minimum-effect claim is a one-sided test at 5% and the
# equivalence claim a TOST at 5%, both read off the 90% interval (Lakens,
# 2017; Lakens et al., 2018; Riesthuis, 2024 for equivalence). Before 1.8.0
# this test matched TOSTER only at alpha = .025 (the 95%-interval mapping)
# and even guarded that the .05 mapping did not match -- i.e. it verified
# the engine's rule rather than auditing the rule's convention. Lesson kept
# here on purpose: a convention-matched check cannot flag a convention
# choice; the audit has to pin the DEFAULT to the field's default.
# TOSTER is pooled/fixed-n, our rule unpooled/random assignment, so exact
# equality is not expected; ~0.01 agreement is.

test_that("minimum-effect claim at the default alpha matches TOSTER at ITS default alpha", {
  skip_if_not_installed("TOSTER")
  cl <- powerape:::claim_levels("minimum")$conf_claim
  expect_equal(cl, 0.90)
  ex <- exact_power_sat(712, 0.5, 0.30, 0.40, cl, "minimum", 0.04)
  tm <- TOSTER::power_twoprop(p1 = 0.40, p2 = 0.30, n = 356, null = 0.04,
                              alpha = 0.05, alternative = "one.sided")$power
  expect_lt(abs(ex - tm), 0.01)
  ## the conservative 95%-interval mapping (one-sided 2.5%) must NOT match
  ## the default -- that was the pre-1.8.0 behavior
  t25 <- TOSTER::power_twoprop(p1 = 0.40, p2 = 0.30, n = 356, null = 0.04,
                               alpha = 0.025, alternative = "one.sided")$power
  expect_gt(abs(ex - t25), 0.05)
})

test_that("equivalence claim at the default alpha matches TOSTER's TOST at ITS default", {
  skip_if_not_installed("TOSTER")
  cl <- powerape:::claim_levels("equivalence")$conf_claim
  expect_equal(cl, 0.90)
  ex <- exact_power_sat(700, 0.5, 0.30, 0.30, cl, "equivalence", 0.10)
  tm <- TOSTER::power_twoprop(p1 = 0.30, p2 = 0.30, n = 350, null = 0.10,
                              alpha = 0.05, alternative = "equivalence")$power
  expect_lt(abs(ex - tm), 0.015)
  t25 <- TOSTER::power_twoprop(p1 = 0.30, p2 = 0.30, n = 350, null = 0.10,
                               alpha = 0.025, alternative = "equivalence")$power
  expect_gt(abs(ex - t25), 0.04)
})

test_that("detection stays the two-sided 5% test (95% interval)", {
  expect_equal(powerape:::claim_levels("detect")$conf_claim, 0.95)
  expect_equal(powerape:::claim_levels("detect")$alpha, 0.05)
})

test_that("conf override reproduces the conservative pre-1.8.0 convention", {
  skip_if_not_installed("TOSTER")
  lv <- powerape:::claim_levels("minimum", conf = 0.95)
  expect_equal(lv$alpha, 0.025)
  ex <- exact_power_sat(712, 0.5, 0.30, 0.40, 0.95, "minimum", 0.04)
  t25 <- TOSTER::power_twoprop(p1 = 0.40, p2 = 0.30, n = 356, null = 0.04,
                               alpha = 0.025, alternative = "one.sided")$power
  expect_lt(abs(ex - t25), 0.01)
  ## and the public functions carry alpha/conf through
  d <- set_ape(ape_dgp("probit", focal = pa_var("t", "binary", p = 0.5),
                       baseline = 0.30), 0.10)
  pw90 <- ape_power(d, 400, claim = "minimum", sesoi = 0.04, nsim = 60, seed = 1)
  pw95 <- ape_power(d, 400, claim = "minimum", sesoi = 0.04, conf = 0.95,
                    nsim = 60, seed = 1)
  expect_equal(pw90$conf, 0.90); expect_equal(pw90$alpha, 0.05)
  expect_equal(pw95$conf, 0.95); expect_equal(pw95$alpha, 0.025)
  expect_gte(pw90$power, pw95$power)
})
