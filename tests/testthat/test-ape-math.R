test_that("saturated-case APE equals the difference in proportions and its SE the unpooled formula", {
  set.seed(1)
  n <- 4000
  d <- rbinom(n, 1, 0.5)
  y <- rbinom(n, 1, 0.30 + 0.10 * d)
  X <- cbind(1, d)

  ft <- powerape:::fit_index_model(X, y, "probit")
  expect_true(ft$ok)
  est <- powerape:::ape_binary_est(ft$fit$coefficients, X, 2L, "probit")

  p1 <- mean(y[d == 1])
  p0 <- mean(y[d == 0])
  expect_equal(est$ape, p1 - p0, tolerance = 1e-6)

  V <- powerape:::vcov_from_glmfit(ft$fit)
  se <- sqrt(drop(t(est$jac) %*% V %*% est$jac))
  n1 <- sum(d == 1)
  n0 <- n - n1
  se_ref <- sqrt(p1 * (1 - p1) / n1 + p0 * (1 - p0) / n0)
  expect_equal(se, se_ref, tolerance = 1e-5)
})

test_that("the same equivalences hold for logit", {
  set.seed(2)
  n <- 4000
  d <- rbinom(n, 1, 0.5)
  y <- rbinom(n, 1, 0.25 + 0.15 * d)
  X <- cbind(1, d)

  ft <- powerape:::fit_index_model(X, y, "logit")
  expect_true(ft$ok)
  est <- powerape:::ape_binary_est(ft$fit$coefficients, X, 2L, "logit")

  p1 <- mean(y[d == 1])
  p0 <- mean(y[d == 0])
  expect_equal(est$ape, p1 - p0, tolerance = 1e-6)

  V <- powerape:::vcov_from_glmfit(ft$fit)
  se <- sqrt(drop(t(est$jac) %*% V %*% est$jac))
  n1 <- sum(d == 1)
  n0 <- n - n1
  se_ref <- sqrt(p1 * (1 - p1) / n1 + p0 * (1 - p0) / n0)
  expect_equal(se, se_ref, tolerance = 1e-5)
})

## The two tests below came out of a mutation pilot (validation-gaps
## review, 2026-10-14): a wrong logit second derivative, a zeroed covariate
## block of the gradient, the probit density in a logit Jacobian, and HC1
## in place of HC0 all survived the suite before them.

test_that("every analytic delta-method gradient equals its finite-difference gradient", {
  fd_jac <- function(f, beta, h = 1e-6) {
    vapply(seq_along(beta), function(j) {
      bp <- beta; bm <- beta; bp[j] <- bp[j] + h; bm[j] <- bm[j] - h
      (f(bp)$ape - f(bm)$ape) / (2 * h)
    }, numeric(1))
  }
  relerr <- function(f, beta) {
    nu <- fd_jac(f, beta)
    max(abs(f(beta)$jac - nu)) / max(abs(nu), 1e-8)
  }
  set.seed(101)
  n <- 600
  d <- rbinom(n, 1, 0.45); w <- rnorm(n); m <- rbinom(n, 1, 0.4)
  v <- rnorm(n, 0.3, 1.2); z1 <- rnorm(n); z2 <- rbinom(n, 1, 0.3)
  err <- c()
  types <- list(c("binary", "binary"), c("normal", "binary"),
                c("binary", "normal"), c("normal", "normal"))
  for (lk in c("probit", "logit")) {
    X <- cbind(1, d, z1, z2)
    err[paste(lk, "APE binary")] <- relerr(
      function(bb) powerape:::ape_binary_est(bb, X, 2L, lk), c(-0.4, 0.5, 0.3, -0.2))
    Xc <- cbind(1, w, z1, z2)
    err[paste(lk, "APE continuous")] <- relerr(
      function(bb) powerape:::ape_cont_est(bb, Xc, 2L, lk), c(-0.4, 0.35, 0.3, -0.2))
    for (tp in types) {
      f1 <- if (tp[1] == "binary") d else w
      f2 <- if (tp[2] == "binary") m else v
      Xi <- cbind(1, f1, f2, f1 * f2, z1, z2)
      err[paste(lk, "AIE", tp[1], "x", tp[2])] <- relerr(
        function(bb) powerape:::aie_est(bb, Xi, lk, tp[1], tp[2]),
        c(-0.3, 0.4, 0.25, 0.3, 0.2, -0.1))
    }
  }
  ## control-function estimators, on fits built by cf_fit (linear and probit
  ## first stages, with and without a moderator)
  zi <- rnorm(n)
  u <- rnorm(n); e <- 0.5 * u + sqrt(0.75) * rnorm(n)
  y2c <- 0.6 * zi + 0.3 * z1 + e
  y2b <- as.integer(0.2 + 0.7 * zi + 0.3 * z1 + e > 0)
  y0c <- as.integer(-0.3 + 0.4 * y2c + 0.3 * z1 + 0.2 * m + 0.2 * y2c * m + u > 0)
  y0b <- as.integer(-0.3 + 0.5 * y2b + 0.3 * z1 + 0.2 * m + 0.2 * y2b * m + u > 0)
  X1 <- cbind(1, z1, zi)
  cf <- powerape:::cf_fit
  fc <- cf(y0c, y2c, Z = cbind(z1), X1 = X1, first = "linear")
  fb <- cf(y0b, y2b, Z = cbind(z1), X1 = X1, first = "probit")
  fcm <- cf(y0c, y2c, m = m, Z = cbind(z1), X1 = X1, first = "linear")
  fbm <- cf(y0b, y2b, m = m, Z = cbind(z1), X1 = X1, first = "probit")
  cfa <- powerape:::cf_ape_est
  cfi <- powerape:::cf_aie_est
  err["CF APE continuous"] <- relerr(function(th) cfa(th, fc$W, "normal", fc$cf_col), fc$theta)
  err["CF APE binary"] <- relerr(function(th) cfa(th, fb$W, "binary", fb$cf_col), fb$theta)
  err["CF APE continuous, moderated"] <- relerr(
    function(th) cfa(th, fcm$W, "normal", fcm$cf_col, fcm$cfm_col), fcm$theta)
  err["CF APE binary, moderated"] <- relerr(
    function(th) cfa(th, fbm$W, "binary", fbm$cf_col, fbm$cfm_col), fbm$theta)
  err["CF AIE continuous x binary"] <- relerr(
    function(th) cfi(th, fcm$W, "normal", fcm$cf_col, fcm$cfm_col), fcm$theta)
  err["CF AIE binary x binary"] <- relerr(
    function(th) cfi(th, fbm$W, "binary", fbm$cf_col, fbm$cfm_col), fbm$theta)
  expect_length(err, 18L)
  expect_lt(max(err), 1e-6, label = sprintf("worst gradient error (%s)", names(which.max(err))))
})

test_that("se = 'robust' is HC0 and the cluster sandwich is vcovCL to machine precision", {
  skip_if_not_installed("sandwich")
  ## a tight convergence removes the one-iteration-stale IRLS weights glm
  ## stores, so the comparison can be exact rather than 1%
  set.seed(5)
  n <- 1500
  df <- data.frame(x1 = rnorm(n), x2 = rbinom(n, 1, 0.4))
  df$y <- rbinom(n, 1, pnorm(-0.3 + 0.5 * df$x1 - 0.4 * df$x2))
  ctl <- glm.control(epsilon = 1e-15, maxit = 200)
  for (lk in c("probit", "logit")) {
    gm <- glm(y ~ x1 + x2, binomial(lk), data = df, control = ctl)
    V <- powerape:::vcov_sandwich(gm, model.matrix(gm), lk)
    Vh <- sandwich::vcovHC(gm, type = "HC0")
    expect_lt(max(abs(V - Vh) / abs(Vh)), 1e-6)
  }
  gm <- glm(y ~ x1 + x2, binomial("probit"), data = df, control = ctl)
  cid <- rep(1:100, each = 15)
  Vc <- powerape:::vcov_cluster(gm, model.matrix(gm), cid, "probit")
  Vch <- sandwich::vcovCL(gm, cluster = cid, type = "HC0", cadjust = TRUE)
  expect_lt(max(abs(Vc - Vch) / abs(Vch)), 1e-6)
})
