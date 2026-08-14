# Control-function estimator: regression tests against Stata 18.5's
# cfprobit (StataNow), verified live on 2026-08-14. Coefficient targets
# agree to ML convergence tolerance (~1e-6 relative); SE targets to ~1e-6
# (the stacked vcov conventions were identified to ~1e-7 -- observed score
# derivatives, generated regressors as fixed instruments in the
# first-stage Jacobian, no small-sample factors).

cf_fixture <- local({
  set.seed(20260814)
  n <- 1500
  x1 <- rnorm(n); z1 <- rnorm(n); z2 <- rnorm(n)
  cid <- rep(1:75, each = 20)
  rho <- 0.5
  v <- rnorm(n); u <- rho * v + sqrt(1 - rho^2) * rnorm(n)
  y2 <- 0.3 + 0.4 * x1 + 0.5 * z1 + 0.3 * z2 + v
  y0 <- as.integer(-0.3 + 0.5 * y2 + 0.4 * x1 + u > 0)
  y2b <- as.integer(0.2 + 0.4 * x1 + 0.6 * z1 + 0.4 * z2 + v > 0)
  y0b <- as.integer(-0.3 + 0.7 * y2b + 0.4 * x1 + u > 0)
  y2f <- pnorm(0.2 + 0.4 * x1 + 0.5 * z1 + 0.5 * v)
  y0f <- as.integer(-0.3 + 0.8 * y2f + 0.4 * x1 + u > 0)
  y2p <- rpois(n, exp(0.2 + 0.3 * x1 + 0.4 * z1 + 0.3 * v))
  y0p <- as.integer(-0.6 + 0.25 * y2p + 0.4 * x1 + u > 0)
  m <- rbinom(n, 1, 0.4)
  y0i <- as.integer(-0.4 + 0.6 * y2b + 0.3 * m + 0.4 * y2b * m + 0.4 * x1 + u > 0)
  list(x1 = x1, z1 = z1, z2 = z2, cid = cid, y2 = y2, y0 = y0,
       y2b = y2b, y0b = y0b, y2f = y2f, y0f = y0f, y2p = y2p, y0p = y0p,
       m = m, y0i = y0i)
})

pa_se_cf <- function(e, V) sqrt(drop(t(e$jac) %*% V %*% e$jac))

test_that("CF linear first stage matches cfprobit exactly (robust + cluster)", {
  fx <- cf_fixture
  X1 <- cbind(1, fx$x1, fx$z1, fx$z2)
  f <- powerape:::cf_fit(fx$y0, fx$y2, Z = cbind(fx$x1), X1 = X1,
                         first = "linear")
  ## Stata: theta (cons, y2, x1, cf) and robust SEs
  expect_equal(unname(f$theta), c(-0.3350216, 0.6105196, 0.4359074, 0.5545095),
               tolerance = 1e-5)
  expect_equal(unname(sqrt(diag(f$V_main))),
               c(0.0499661, 0.0759271, 0.0559708, 0.0861535), tolerance = 1e-5)
  e <- powerape:::cf_ape_est(f$theta, f$W, "normal", f$cf_col)
  expect_equal(e$ape, 0.1412121, tolerance = 1e-6)
  expect_equal(pa_se_cf(e, f$V_main), 0.0163633, tolerance = 1e-5)

  fc <- powerape:::cf_fit(fx$y0, fx$y2, Z = cbind(fx$x1), X1 = X1,
                          first = "linear", id = fx$cid)
  expect_equal(unname(sqrt(diag(fc$V_main))),
               c(0.0500113, 0.0723990, 0.0542260, 0.0885456), tolerance = 1e-5)
  ec <- powerape:::cf_ape_est(fc$theta, fc$W, "normal", fc$cf_col)
  expect_equal(pa_se_cf(ec, fc$V_main), 0.0160701, tolerance = 1e-5)
})

test_that("CF probit and poisson first stages match cfprobit", {
  fx <- cf_fixture
  X1 <- cbind(1, fx$x1, fx$z1, fx$z2)
  f <- powerape:::cf_fit(fx$y0b, fx$y2b, Z = cbind(fx$x1), X1 = X1,
                         first = "probit")
  expect_equal(unname(f$theta), c(-0.3293232, 0.7454643, 0.4767555, 0.5144434),
               tolerance = 1e-4)
  expect_equal(unname(sqrt(diag(f$V_main))),
               c(0.0989886, 0.1628248, 0.0476460, 0.1135337), tolerance = 1e-5)
  e <- powerape:::cf_ape_est(f$theta, f$W, "binary", f$cf_col)
  expect_equal(e$ape, 0.2486321, tolerance = 1e-5)
  expect_equal(pa_se_cf(e, f$V_main), 0.0602349, tolerance = 1e-5)

  fp <- powerape:::cf_fit(fx$y0p, fx$y2p, Z = cbind(fx$x1),
                          X1 = cbind(1, fx$x1, fx$z1), first = "poisson")
  ep <- powerape:::cf_ape_est(fp$theta, fp$W, "normal", fp$cf_col)
  expect_equal(ep$ape, 0.0882267, tolerance = 1e-5)
  expect_equal(pa_se_cf(ep, fp$V_main), 0.0168097, tolerance = 1e-5)
})

test_that("CF fprobit matches cfprobit's fit; margins gap is Stata's remake bug", {
  fx <- cf_fixture
  X1 <- cbind(1, fx$x1, fx$z1)
  f <- powerape:::cf_fit(fx$y0f, fx$y2f, Z = cbind(fx$x1), X1 = X1,
                         first = "fprobit")
  expect_equal(unname(f$theta), c(-0.1306443, 0.5796992, 0.5199608, 2.5120075),
               tolerance = 1e-5)
  expect_equal(unname(sqrt(diag(f$V_main))),
               c(0.1497125, 0.2529033, 0.0579644, 0.2358091), tolerance = 1e-5)
  ## powerape's ASF uses the estimation-time generalized residuals
  e <- powerape:::cf_ape_est(f$theta, f$W, "normal", f$cf_col)
  expect_equal(e$ape, 0.1618274, tolerance = 1e-6)
  ## Stata's margins (0.0438960) instead rebuilds the CF with the BINARY
  ## branch of _remake_cfs -- lambda(q) for every y > 0 -- which is
  ## inconsistent with its own estimator for fractional y. Attribution:
  fs <- powerape:::cf_first(fx$y2f, X1, "fprobit")
  q <- drop(X1 %*% fs$pi)
  c_remake <- ifelse(fx$y2f > 0, dnorm(q) / pnorm(q), -dnorm(q) / pnorm(-q))
  idx_r <- drop(cbind(1, fx$y2f, fx$x1, c_remake) %*% f$theta)
  expect_equal(mean(dnorm(idx_r)) * f$theta[[2L]], 0.0438960, tolerance = 1e-5)
})

test_that("CF AIE with interact(m) matches cfprobit (robust, cluster, z-x-m)", {
  fx <- cf_fixture
  X1m <- cbind(1, fx$x1, fx$m, fx$z1, fx$z2)
  f <- powerape:::cf_fit(fx$y0i, fx$y2b, m = fx$m, Z = cbind(fx$x1),
                         X1 = X1m, first = "probit")
  e <- powerape:::cf_aie_est(f$theta, f$W, "binary", f$cf_col, f$cfm_col)
  expect_equal(e$ape, 0.1336071, tolerance = 1e-4)
  expect_equal(pa_se_cf(e, f$V_main), 0.1030098, tolerance = 1e-5)

  fc <- powerape:::cf_fit(fx$y0i, fx$y2b, m = fx$m, Z = cbind(fx$x1),
                          X1 = X1m, first = "probit", id = fx$cid)
  ec <- powerape:::cf_aie_est(fc$theta, fc$W, "binary", fc$cf_col, fc$cfm_col)
  expect_equal(pa_se_cf(ec, fc$V_main), 0.1155685, tolerance = 1e-5)

  ## instrument-by-moderator interactions in the first stage
  X1i <- cbind(1, fx$x1, fx$m, fx$z1, fx$z2, fx$z1 * fx$m, fx$z2 * fx$m)
  fi <- powerape:::cf_fit(fx$y0i, fx$y2b, m = fx$m, Z = cbind(fx$x1),
                          X1 = X1i, first = "probit")
  ei <- powerape:::cf_aie_est(fi$theta, fi$W, "binary", fi$cf_col, fi$cfm_col)
  expect_equal(ei$ape, 0.1471267, tolerance = 1e-4)
  expect_equal(pa_se_cf(ei, fi$V_main), 0.1031082, tolerance = 1e-5)
})

test_that("generalized residual is the probit score; first-stage moments vanish", {
  fx <- cf_fixture
  X1 <- cbind(1, fx$x1, fx$z1, fx$z2)
  fs <- powerape:::cf_first(fx$y2b, X1, "probit")
  ## moments c x X1 sum to ~0 at the first-stage optimum (glm.fit's
  ## deviance-based stopping leaves gradients of order 1e-3 at n = 1500)
  expect_lt(max(abs(colSums(fs$c * X1))), 0.01)
  ## analytic dc/dq matches numerical differentiation
  q <- drop(X1 %*% fs$pi)
  h <- 1e-6
  num <- (powerape:::probit_score_bits(q + h, fx$y2b)$s -
            powerape:::probit_score_bits(q - h, fx$y2b)$s) / (2 * h)
  expect_equal(fs$dcdq, num, tolerance = 1e-5)
})

test_that("vcov_sandwich matches sandwich::vcovHC(type = 'HC0')", {
  skip_if_not_installed("sandwich")
  fx <- cf_fixture
  df <- data.frame(y = fx$y0, y2 = fx$y2, x1 = fx$x1)
  gm <- glm(y ~ y2 + x1, binomial("probit"), data = df)
  V <- powerape:::vcov_sandwich(gm, model.matrix(gm), "probit")
  ## the bread is identical to machine precision; the meat differs only
  ## because sandwich's estfun uses glm's stored IRLS working weights
  ## (one iteration stale) while ours is the exact analytic score at the
  ## final coefficients -- agreement to ~0.3% on this fit
  expect_equal(unname(powerape:::vcov_from_glmfit(gm)),
               unname(sandwich::bread(gm) / nrow(df)), tolerance = 1e-10)
  expect_equal(unname(V), unname(sandwich::vcovHC(gm, type = "HC0")),
               tolerance = 0.01)
})
