# Claim-semantics concordance with TOSTER: powerape's 95%-CI minimum-effect
# and equivalence claims must match TOSTER's shifted-null / TOST z-tests at
# alpha = 0.025 (the documented CI-level correspondence) -- an external check
# of the claim definitions that the exact enumeration alone cannot provide.
# TOSTER is pooled/fixed-n, our rule unpooled/random assignment, so exact
# equality is not expected; ~0.01 agreement is.

test_that("minimum-effect claim matches TOSTER's shifted-null test at the matched alpha", {
  skip_if_not_installed("TOSTER")
  ex <- exact_power_sat(712, 0.5, 0.30, 0.40, 0.95, "minimum", 0.04)
  tm <- TOSTER::power_twoprop(p1 = 0.40, p2 = 0.30, n = 356, null = 0.04,
                              alpha = 0.025, alternative = "one.sided")$power
  expect_lt(abs(ex - tm), 0.01)
  ## the naive alpha = .05 mapping must NOT match (guards the documentation)
  tn <- TOSTER::power_twoprop(p1 = 0.40, p2 = 0.30, n = 356, null = 0.04,
                              alpha = 0.05, alternative = "one.sided")$power
  expect_gt(abs(ex - tn), 0.05)
})

test_that("equivalence claim matches TOSTER's TOST at the matched alpha", {
  skip_if_not_installed("TOSTER")
  ex <- exact_power_sat(700, 0.5, 0.30, 0.30, 0.95, "equivalence", 0.10)
  tm <- TOSTER::power_twoprop(p1 = 0.30, p2 = 0.30, n = 350, null = 0.10,
                              alpha = 0.025, alternative = "equivalence")$power
  expect_lt(abs(ex - tm), 0.015)
  tn <- TOSTER::power_twoprop(p1 = 0.30, p2 = 0.30, n = 350, null = 0.10,
                              alpha = 0.05, alternative = "equivalence")$power
  expect_gt(abs(ex - tn), 0.04)
})
