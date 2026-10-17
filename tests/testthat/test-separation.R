# Separation in the cells that identify the estimand (1.11.0; DESIGN.md
# section 16): flagged per replication, scored as failed by default
# ("fail", the field's reference analysis refuses such fits) or kept with the
# degenerate Wald interval ("keep", R's glm + marginaleffects), and reported
# either way. The exact-enumeration certificate of both conventions lives in
# test-exact-power.R.

test_that("cell_check flags empty and full cells, and only those", {
  cc <- powerape:::cell_check
  t <- rep(0:1, each = 150)
  y <- c(rep(0, 150), rep(1, 8), rep(0, 142))          # 0/150 vs 8/150
  expect_true(cc(y, t)$sep)
  expect_identical(cc(y, t)$min_cell, 0L)
  y_full <- c(rep(1, 150), rep(1, 8), rep(0, 142))     # control arm all events
  expect_true(cc(y_full, t)$sep)
  set.seed(1)
  y_ok <- c(rbinom(150, 1, 0.1), rbinom(150, 1, 0.2))
  expect_false(cc(y_ok, t)$sep)
  expect_equal(cc(y_ok, t)$min_cell, min(sum(y_ok[1:150]), sum(y_ok[151:300])))
  ## AIE layout: one empty focal-by-moderator cell is enough
  m <- rep(c(0, 1), times = 150)
  y_aie <- y_ok
  y_aie[t == 1 & m == 1] <- 0
  expect_true(cc(y_aie, t, m)$sep)
  expect_false(cc(y_aie, t)$sep)                      # the arms alone are fine
  ## no binary cell identifies a continuous focal's effect
  expect_false(cc(y, NULL, NULL)$sep)
  expect_true(is.na(cc(y, NULL, NULL)$min_cell))
})

test_that("a separated fit converges silently; the Hauck-Donner backup sees it", {
  t <- rep(0:1, each = 150)
  y <- c(rep(0, 150), rep(1, 8), rep(0, 142))
  X <- cbind(1, t)
  for (lk in c("probit", "logit")) {
    ft <- powerape:::fit_index_model(X, y, lk)
    expect_true(ft$ok)                                # glm.fit "converges"
    expect_true(powerape:::hauck_donner(ft$fit, X, 2L))
  }
  set.seed(2)
  y_ok <- c(rbinom(150, 1, 0.1), rbinom(150, 1, 0.2))
  ft <- powerape:::fit_index_model(X, y_ok, "probit")
  expect_false(powerape:::hauck_donner(ft$fit, X, 2L))
})

test_that("both conventions flag the same replications; only scoring differs", {
  d <- set_ape(ape_dgp("probit", focal = pa_var("t", "binary", p = 0.5),
                       baseline = 0.02), 0.03)
  pf <- suppressWarnings(ape_power(d, n = 100, claim = "detect", nsim = 300,
                                   seed = 3))
  pk <- suppressWarnings(ape_power(d, n = 100, claim = "detect", nsim = 300,
                                   seed = 3, separation = "keep"))
  expect_identical(pf$separation, "fail")
  expect_identical(pk$separation, "keep")
  expect_equal(pf$separated, pk$separated)            # same draws, same flags
  expect_gt(pf$separated, 0.2)
  expect_gte(pf$outcomes[["failed"]], pf$separated)   # counted as failed
  expect_lt(pk$outcomes[["failed"]], pk$separated)    # kept
  expect_lt(pf$power, pk$power)
  expect_equal(sum(pf$outcomes), 1, tolerance = 1e-12)
  ## re-thresholding the stored draws respects the convention
  expect_equal(power_at(pf, 0.05)$power, pf$power)
  expect_equal(power_at(pk, 0.05)$power, pk$power)
  expect_output(print(pf), "had separation")
  expect_output(print(pk), "separation = \"keep\"")
  txt <- unclass(power_statement(pf))
  expect_match(txt, "separation")
  expect_match(unclass(power_statement(pk)), "retained with their Wald")
})

test_that("ape_power warns on separation and sparse cells, and only then", {
  d <- set_ape(ape_dgp("probit", focal = pa_var("t", "binary", p = 0.5),
                       baseline = 0.02), 0.03)
  expect_warning(ape_power(d, n = 100, claim = "detect", nsim = 100, seed = 4),
                 "separation")
  d10 <- set_ape(ape_dgp("probit", focal = pa_var("t", "binary", p = 0.5),
                         baseline = 0.30), 0.10)
  expect_no_warning(pw <- ape_power(d10, n = 400, claim = "detect", nsim = 100,
                                    seed = 4))
  expect_identical(pw$separated, 0)
  expect_gt(pw$min_cell, 10)
  ## sparse but (almost) never separated: ~7 events in the small arm
  ds <- set_ape(ape_dgp("probit", focal = pa_var("t", "binary", p = 0.2),
                        baseline = 0.05), 0.04)
  expect_warning(ape_power(ds, n = 400, claim = "detect", nsim = 100, seed = 5),
                 "averages")
})

test_that("AIE, panel, and IV routes flag separation in their identifying cells", {
  da <- set_aie(ape_dgp("probit", focal = pa_var("t", "binary", p = 0.5),
                        moderator = pa_var("m", "binary", p = 0.2),
                        baseline = 0.03),
                0.02, main_focal = 0.01, main_moderator = 0.01)
  pa <- suppressWarnings(ape_power(da, n = 400, claim = "detect", nsim = 100,
                                   seed = 6))
  expect_gt(pa$separated, 0.05)
  expect_warning(ape_power(da, n = 400, claim = "detect", nsim = 60, seed = 6),
                 "focal-by-moderator cell")

  dp <- set_ape(ape_dgp_panel(focal = pa_var("t", "binary", p = 0.5, icc = 1),
                              n_periods = 2, rho = 0.2, baseline = 0.02,
                              n_int = 2e4), 0.03)
  pp <- suppressWarnings(ape_power(dp, n = 60, claim = "detect", nsim = 80,
                                   seed = 7))
  expect_gt(pp$separated, 0.05)
  expect_gte(pp$outcomes[["failed"]], pp$separated)

  di <- set_ape(ape_dgp_iv(focal = pa_var("t", "binary", p = 0.15),
                           instruments = pa_var("s", "normal"),
                           endogeneity = 0.3, iv_strength = 0.3,
                           baseline = 0.02, n_int = 2e4), 0.01)
  pi <- suppressWarnings(ape_power(di, n = 300, claim = "detect", nsim = 60,
                                   seed = 8))
  expect_gt(pi$separated, 0.05)
})

test_that("fail requires the larger n in a rare, unbalanced design", {
  ## .05 -> .01 with a 10% treated group: under "keep" the separated
  ## studies (zero treated events) count as detections, so the requirement
  ## is roughly half of what the documented rule needs
  d <- set_ape(ape_dgp("probit", focal = pa_var("t", "binary", p = 0.10),
                       baseline = 0.05), -0.04)
  nf <- suppressWarnings(ape_n(d, power = 0.8, claim = "detect", nsim = 300,
                               seed = 9, confirm = FALSE))
  nk <- suppressWarnings(ape_n(d, power = 0.8, claim = "detect", nsim = 300,
                               seed = 9, confirm = FALSE, separation = "keep"))
  expect_gt(nf$n, 1.5 * nk$n)
  expect_identical(nf$separation, "fail")
  expect_true(is.finite(nf$separated))
})

test_that("the AIE cell check reads all four focal-by-moderator cells (standard and IV routes)", {
  ## replay the engine's draws by hand (fits consume no random numbers) and
  ## compare the reported min_cell with the four-cell minimum. This guards
  ## the min_cell < 10 sparse-cell warning; a call that dropped the
  ## moderator cells survived the suite before this test (validation-gaps
  ## review, 2026-10-14).
  cell_min <- function(y, g) {
    nn <- tabulate(g + 1L, nbins = 4L)
    ev <- tabulate(g[y == 1] + 1L, nbins = 4L)
    keep <- nn > 0L
    min(pmin(ev, nn - ev)[keep])
  }
  da <- set_aie(ape_dgp("probit", focal = pa_var("t", "binary", p = 0.5),
                        moderator = pa_var("m", "binary", p = 0.2),
                        baseline = 0.03, n_int = 2e4),
                0.02, main_focal = 0.01, main_moderator = 0.01)
  nsim <- 40
  s <- powerape:::sim_ci(da, 400, nsim, 0.95, seed = 6, separation = "keep")
  bt <- powerape:::beta_true(da)
  set.seed(6)
  four <- arms <- numeric(nsim)
  for (r in seq_len(nsim)) {
    xx <- powerape:::draw_x(da, 400)
    y <- rbinom(400, 1L, da$G(drop(xx$X %*% bt)))
    four[r] <- cell_min(y, 2L * as.integer(xx$X[, 2L]) + as.integer(xx$X[, 3L]))
    arms[r] <- cell_min(y, as.integer(xx$X[, 2L]))
  }
  expect_equal(s$min_cell, four)
  ## the moderator cells bind (otherwise the comparison could not tell them apart)
  expect_true(any(four < arms))

  di <- set_aie(ape_dgp_iv(focal = pa_var("t", "binary", p = 0.4),
                           moderator = pa_var("m", "binary", p = 0.2),
                           instruments = pa_var("z", "binary", p = 0.5),
                           endogeneity = 0.4, iv_strength = 0.25,
                           baseline = 0.05, n_int = 2e4, seed_int = 4L),
                0.03, main_focal = 0.03, main_moderator = 0.01)
  si <- powerape:::sim_ci(di, 500, 20, 0.95, seed = 7, separation = "keep")
  set.seed(7)
  four_iv <- arms_iv <- numeric(20)
  for (r in seq_len(20)) {
    xx <- powerape:::draw_x(di, 500)
    four_iv[r] <- cell_min(xx$y, 2L * as.integer(xx$d) + as.integer(xx$m))
    arms_iv[r] <- cell_min(xx$y, as.integer(xx$d))
  }
  expect_equal(si$min_cell, four_iv)
  expect_true(any(four_iv < arms_iv))
})
