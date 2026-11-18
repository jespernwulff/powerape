# Three-way interaction designs. A focal variable, a
# moderator and a second, binary moderator w; the estimand is the two-way
# AIE of (focal, moderator) at w = 1 minus the same at w = 0 (ginteff's
# three-way effect). Anchors: ginteff on fixed data, finite differences,
# the saturated all-binary design, and the design's own integration draw.

se3 <- function(est, V) sqrt(drop(t(est$jac) %*% V %*% est$jac))

v3 <- list(
  fb = pa_var("treat", "binary", p = 0.5),
  fc = pa_var("rd", "normal", mean = 2, sd = 1.5),
  mb = pa_var("female", "binary", p = 0.45),
  mc = pa_var("size", "normal", mean = -1, sd = 2),
  w  = pa_var("public", "binary", p = 0.4))

# The four type pairs of (focal, moderator), with anchors scaled to stay
# feasible for the continuous variables.
pairs3 <- list(
  bb = list(f = v3$fb, m = v3$mb, main = c(0.08, 0.04, 0.03),
            tw = c(0.03, 0.02, -0.01), target = 0.05),
  cb = list(f = v3$fc, m = v3$mb, main = c(0.03, 0.04, 0.03),
            tw = c(0.015, 0.01, -0.005), target = 0.025),
  bc = list(f = v3$fb, m = v3$mc, main = c(0.08, 0.012, 0.03),
            tw = c(0.015, 0.01, -0.005), target = 0.025),
  cc = list(f = v3$fc, m = v3$mc, main = c(0.03, 0.012, 0.03),
            tw = c(0.0075, 0.005, -0.0025), target = 0.0125))

build3 <- function(s, link = "probit", n_int = 4e4) {
  d <- ape_dgp(model = link, focal = s$f, moderator = s$m, moderator2 = v3$w,
               covariates = list(pa_var("age", "normal", mean = 40, sd = 10)),
               correlation = 0.2, baseline = 0.30, signal = 0.15, n_int = n_int)
  set_aie3(d, target = s$target, main = s$main, two_way = s$tw)
}

# Design matrix in the package's layout from the three variables.
X3 <- function(d, m, w, Z = NULL) cbind(1, d, m, d * m, w, d * w, m * w, d * m * w, Z)

bbb <- function(link = "probit", target = 0.05) {
  set_aie3(ape_dgp(model = link, focal = v3$fb, moderator = v3$mb,
                   moderator2 = v3$w, baseline = 0.30),
           target = target, main = c(0.08, 0.04, 0.02), two_way = c(0.03, -0.01, 0.02))
}

test_that("the builder and the pinning function guard their inputs", {
  expect_error(ape_dgp(focal = v3$fb, moderator2 = v3$w, baseline = 0.3),
               "needs a `moderator`")
  expect_error(ape_dgp(focal = v3$fb, moderator = v3$mb, moderator2 = v3$mc,
                       baseline = 0.3), "must be a binary")
  expect_error(ape_dgp(focal = v3$fb, moderator = v3$mb, moderator2 = v3$mb,
                       baseline = 0.3), "three different variables")
  expect_error(ape_dgp(focal = v3$fb, moderator = v3$mb, moderator2 = v3$w,
                       correlation = diag(2), baseline = 0.3),
               "3 x 3 matrix over (focal, moderator, moderator2", fixed = TRUE)
  d0 <- ape_dgp(focal = v3$fb, moderator = v3$mb, moderator2 = v3$w, baseline = 0.3)
  expect_error(set_aie3(d0, 0.05), "needs every anchor")
  expect_error(set_aie3(d0, 0.05, main = c(0.1, 0.1), two_way = c(0, 0, 0)),
               "`main` must hold three numbers")
  expect_error(set_aie3(d0, 0.05, main = c(0.1, 0.1, 0.1), two_way = 0),
               "`two_way` must hold three numbers")
  expect_error(set_aie3(d0, 0.9, main = c(0.08, 0.04, 0.02), two_way = c(0.03, 0, 0)),
               "three-way AIE 0.9000 is infeasible")
  expect_error(set_aie3(d0, 0.05, main = c(0.08, 0.04, 0.02), two_way = c(0.7, 0, 0)),
               "two-way AIE (treat x female) 0.7000 is infeasible", fixed = TRUE)
  expect_error(set_aie3(d0, 0.05, main = c(0.8, 0.04, 0.02), two_way = c(0.03, 0, 0)),
               "main-effect APE (treat)", fixed = TRUE)
  ## a two-way design is not a three-way design, and the other way round
  d2 <- ape_dgp(focal = v3$fb, moderator = v3$mb, baseline = 0.3)
  expect_error(set_aie3(d2, 0.05, main = c(0.1, 0.1, 0.1), two_way = c(0, 0, 0)),
               "no second moderator")
  expect_error(set_aie(d0, 0.05, 0.1, 0.1), "set_aie3()", fixed = TRUE)
  expect_error(set_ape(d0, 0.05), "set_aie3()", fixed = TRUE)
  expect_error(true_aie3(d0), "call set_aie3() first", fixed = TRUE)
  expect_error(ape_power(d0, 500, "detect", nsim = 20), "call set_aie3() first",
               fixed = TRUE)
  expect_output(print(d0), "call set_aie3()", fixed = TRUE)
})

test_that("functions outside the wired set refuse a three-way design", {
  d <- bbb()
  expect_error(true_aie(d), "true_aie3()", fixed = TRUE)
  expect_error(true_ape(d), "true_aie3()", fixed = TRUE)
  expect_error(ape_robust(d, n = 1000, claim = "detect"),
               "not available for three-way interaction designs")
  p <- suppressWarnings(ape_power(d, 800, "detect", nsim = 20, seed = 1))
  expect_error(power_at(p, 0.01), "not available for three-way")
  expect_error(ape_alpha(p), "not available for three-way")
  expect_error(power_statement(p), "not available for three-way")
  expect_error(ape_power(d, 800, "minimum", sesoi = 0.08, nsim = 20),
               "true three-way AIE")
})

test_that("with three binary variables the anchors are the eight cell rates", {
  for (lk in c("probit", "logit")) {
    d <- bbb(lk)
    ## from the anchors by hand: baseline .30; mains .08/.04/.02; two-way
    ## .03/-.01/.02 at the third variable's reference; three-way .05
    want <- array(NA_real_, c(2, 2, 2))
    want[1, 1, 1] <- 0.30
    want[2, 1, 1] <- 0.38; want[1, 2, 1] <- 0.34; want[1, 1, 2] <- 0.32
    want[2, 2, 1] <- 0.38 + 0.34 - 0.30 + 0.03
    want[2, 1, 2] <- 0.38 + 0.32 - 0.30 - 0.01
    want[1, 2, 2] <- 0.34 + 0.32 - 0.30 + 0.02
    want[2, 2, 2] <- 0.05 + (want[2, 2, 1] - want[1, 2, 1] - want[2, 1, 1] + want[1, 1, 1]) +
      want[1, 2, 2] + want[2, 1, 2] - want[1, 1, 2]
    expect_equal(d$cells3, want, tolerance = 1e-7)
    ## the same rates from the raw coefficients on the design-matrix rows
    g <- expand.grid(d = 0:1, m = 0:1, w = 0:1)
    pr <- d$G(drop(X3(g$d, g$m, g$w) %*% powerape:::beta_true(d)))
    expect_equal(pr, as.vector(want), tolerance = 1e-7)
    expect_equal(true_aie3(d), 0.05, tolerance = 1e-7)
    expect_identical(d$estimand, "aie3")
  }
})

test_that("anchors implied by an index without products solve to zero products", {
  ## three binary variables, each shifting the index by 0.5: the implied
  ## main, two-way and three-way effects in points are not zero, and
  ## pinning them must return the additive index
  for (lk in c("probit", "logit")) {
    lf <- powerape:::link_funs(lk)
    r <- function(k) lf$G(lf$Ginv(0.10) + 0.5 * k)   # rate with k variables at 1
    m1 <- r(1) - r(0)
    tw <- r(2) - 2 * r(1) + r(0)
    th <- (r(3) - 2 * r(2) + r(1)) - tw
    d <- set_aie3(ape_dgp(model = lk, focal = v3$fb, moderator = v3$mb,
                          moderator2 = v3$w, baseline = 0.10),
                  th, main = rep(m1, 3), two_way = rep(tw, 3))
    expect_lt(max(abs(c(d$beta_int, d$beta_int_f2, d$beta_int_m2, d$beta_int3))), 1e-6)
    expect_lt(max(abs(c(d$beta_focal, d$beta_mod, d$beta_mod2) - 0.5)), 1e-6)
    ## a zero triple coefficient is not a zero three-way effect in points
    if (lk == "probit") expect_equal(th, -0.0293, tolerance = 1e-2)
  }
})

test_that("every anchor is reproduced from the solved coefficients", {
  for (k in names(pairs3)) for (lk in c("probit", "logit")) {
    s <- pairs3[[k]]
    d <- build3(s, lk)
    a <- powerape:::aie3_anchors(d)
    expect_equal(a$main, s$main, tolerance = 1e-6, info = paste(k, lk))
    expect_equal(a$two_way, s$tw, tolerance = 1e-6, info = paste(k, lk))
    expect_equal(a$target, s$target, tolerance = 1e-6, info = paste(k, lk))
    expect_equal(true_aie3(d), s$target, tolerance = 1e-6, info = paste(k, lk))
    expect_equal(d$main3, s$main)
    expect_equal(d$two_way3, s$tw)
  }
})

test_that("the estimator at the true coefficients returns the pinned effect", {
  ## the check of the centred-to-raw coefficient map: continuous variables
  ## have nonzero means, and everything is correlated
  for (k in names(pairs3)) for (lk in c("probit", "logit")) {
    s <- pairs3[[k]]
    d <- build3(s, lk)
    bt <- powerape:::beta_true(d)
    ## (a) on the design's own integration draw: exact
    rd <- powerape:::reg_draw_3way(d$focal, d$moderator, d$moderator2,
                                   d$covariates, d$R, d$n_int, d$seed_int)
    e <- powerape:::aie3_est(bt, X3(rd$xd, rd$md, rd$wd, rd$Z), lk,
                             s$f$type, s$m$type)
    expect_equal(e$ape, s$target, tolerance = 1e-6, info = paste(k, lk))
    ## (b) on a fresh draw from the simulation's own generator
    xx <- powerape:::with_seed(11, powerape:::draw_x(d, 1e5))
    expect_equal(ncol(xx$X), 9L)
    expect_equal(xx$X[, 8], xx$X[, 2] * xx$X[, 3] * xx$X[, 5])
    expect_equal(xx$X[, 6], xx$X[, 2] * xx$X[, 5])
    expect_equal(xx$X[, 7], xx$X[, 3] * xx$X[, 5])
    e2 <- powerape:::aie3_est(bt, xx$X, lk, s$f$type, s$m$type)
    expect_lt(abs(e2$ape - s$target), 0.0015)
  }
})

test_that("the three-way effect and its SE match ginteff on fixed data", {
  skip_if_not_installed("ginteff")
  set.seed(31)
  n <- 1500
  d <- rbinom(n, 1, 0.5); m <- rbinom(n, 1, 0.45); w <- rbinom(n, 1, 0.5)
  c1 <- rnorm(n); c2 <- rnorm(n, 1, 2); x <- rnorm(n)
  idx <- function(a, b) -0.7 + 0.35 * a + 0.2 * b + 0.2 * w + 0.15 * a * b -
    0.1 * a * w + 0.1 * b * w + 0.2 * a * b * w + 0.4 * x
  run <- function(a, b, ta, tb, link) {
    G <- if (link == "probit") pnorm else plogis
    y <- rbinom(n, 1, G(idx(a, b)))
    df <- data.frame(y, a = if (ta == "binary") factor(a) else a,
                     b = if (tb == "binary") factor(b) else b, w = factor(w), x)
    fit <- glm(y ~ a * b * w + x, binomial(link), data = df)
    g <- ginteff::ginteff(fit, dydxs = c("a", "b", "w"))
    X <- X3(a, b, w, x)
    ft <- powerape:::fit_index_model(X, y, link)
    e <- powerape:::aie3_est(ft$fit$coefficients, X, link, ta, tb)
    list(e = e$ape, se = se3(e, powerape:::vcov_from_glmfit(ft$fit)),
         g = unname(g$aie), gse = unname(g$se))
  }
  cases <- list(list(d, m, "binary", "binary", "probit"),
                list(c1, m, "normal", "binary", "probit"),
                list(d, c1, "binary", "normal", "logit"),
                list(c1, (c2 - 1) / 2, "normal", "normal", "probit"))
  for (cs in cases) {
    r <- do.call(run, cs)
    ## two continuous variables: ginteff takes the cross-partial by
    ## numerical differentiation, which agrees to five digits (1.3e-5
    ## relative here); the other pairs agree to nine digits or more
    tol <- if (cs[[3]] == "normal" && cs[[4]] == "normal") 1e-4 else 1e-7
    expect_equal(r$e, r$g, tolerance = tol, info = paste(cs[[3]], cs[[4]]))
    expect_equal(r$se, r$gse, tolerance = tol, info = paste(cs[[3]], cs[[4]]))
  }
})

test_that("the Jacobian equals finite differences, and the effect is symmetric", {
  set.seed(32)
  n <- 400
  db <- rbinom(n, 1, 0.5); mb <- rbinom(n, 1, 0.4); w <- rbinom(n, 1, 0.5)
  dc <- rnorm(n, 1, 1.2); mc <- rnorm(n, -0.5, 0.8); z <- rnorm(n)
  beta <- c(-0.6, 0.3, 0.2, 0.15, 0.25, -0.1, 0.12, 0.18, 0.4)
  for (lk in c("probit", "logit")) {
    for (cs in list(list(db, mb, "binary", "binary"), list(dc, mb, "normal", "binary"),
                    list(db, mc, "binary", "normal"), list(dc, mc, "normal", "normal"))) {
      X <- X3(cs[[1]], cs[[2]], w, z)
      e <- powerape:::aie3_est(beta, X, lk, cs[[3]], cs[[4]])
      fd <- vapply(seq_along(beta), function(j) {
        h <- 1e-6
        bp <- beta; bp[j] <- bp[j] + h
        bm <- beta; bm[j] <- bm[j] - h
        (powerape:::aie3_est(bp, X, lk, cs[[3]], cs[[4]])$ape -
           powerape:::aie3_est(bm, X, lk, cs[[3]], cs[[4]])$ape) / (2 * h)
      }, numeric(1))
      expect_equal(e$jac, fd, tolerance = 1e-6, info = paste(lk, cs[[3]], cs[[4]]))
      ## exchange the roles of the focal variable and the moderator
      Xs <- X3(cs[[2]], cs[[1]], w, z)
      bs <- beta[c(1, 3, 2, 4, 5, 7, 6, 8, 9)]
      es <- powerape:::aie3_est(bs, Xs, lk, cs[[4]], cs[[3]])
      expect_equal(es$ape, e$ape, tolerance = 1e-12)
    }
    ## all binary: exchange the moderator and the second moderator
    X <- X3(db, mb, w, z)
    e <- powerape:::aie3_est(beta, X, lk, "binary", "binary")
    Xs <- X3(db, w, mb, z)
    bs <- beta[c(1, 2, 5, 6, 3, 4, 7, 8, 9)]
    expect_equal(powerape:::aie3_est(bs, Xs, lk, "binary", "binary")$ape, e$ape,
                 tolerance = 1e-12)
  }
})

test_that("saturated design: triple difference of cell means, binomial SE", {
  set.seed(33)
  n <- 3000
  d <- rbinom(n, 1, 0.5); m <- rbinom(n, 1, 0.45); w <- rbinom(n, 1, 0.4)
  y <- rbinom(n, 1, plogis(-0.8 + 0.4 * d + 0.3 * m + 0.2 * w + 0.3 * d * m * w))
  X <- X3(d, m, w)
  cellmean <- function(a, b, c) mean(y[d == a & m == b & w == c])
  celln <- function(a, b, c) sum(d == a & m == b & w == c)
  g <- expand.grid(a = 0:1, b = 0:1, c = 0:1)
  pm <- mapply(cellmean, g$a, g$b, g$c)
  nn <- mapply(celln, g$a, g$b, g$c)
  sgn <- ifelse((g$a + g$b + g$c) %% 2 == 1, 1, -1)   # +1 for 111, 100, 010, 001
  triple <- sum(sgn * pm)
  se_bin <- sqrt(sum(pm * (1 - pm) / nn))
  for (lk in c("probit", "logit")) {
    ft <- powerape:::fit_index_model(X, y, lk)
    e <- powerape:::aie3_est(ft$fit$coefficients, X, lk, "binary", "binary")
    expect_equal(e$ape, triple, tolerance = 1e-6)
    expect_equal(se3(e, powerape:::vcov_from_glmfit(ft$fit)), se_bin, tolerance = 1e-5)
  }
  ## the linear analysis: the triple product's coefficient is the same
  ## number, and its HC1 variance the same sum times n / (n - 8)
  fo <- powerape:::fit_ols(X, y)
  expect_equal(fo$coefficients[[8]], triple, tolerance = 1e-10)
  expect_equal(sqrt(powerape:::vcov_ols_hc1(fo, X)[8, 8]),
               se_bin * sqrt(n / (n - 8)), tolerance = 1e-8)
})

test_that("the cell check sees every cell of the binary layout", {
  g <- expand.grid(d = 0:1, m = 0:1, w = 0:1)
  g <- g[rep(seq_len(8), each = 20), ]
  X <- X3(g$d, g$m, g$w)
  y <- rep(c(rep(1, 6), rep(0, 14)), 8)
  cc <- powerape:::aie3_cell_check(y, X, TRUE, TRUE)
  expect_false(cc$sep)
  expect_equal(cc$min_cell, 6)
  ## no events in the (1, 1, 1) cell only
  y2 <- y
  y2[g$d == 1 & g$m == 1 & g$w == 1] <- 0
  expect_true(powerape:::aie3_cell_check(y2, X, TRUE, TRUE)$sep)
  ## with a continuous focal variable the cells are moderator x second
  ## moderator, and that pooled cell still holds events
  cc2 <- powerape:::aie3_cell_check(y2, X, FALSE, TRUE)
  expect_false(cc2$sep)
  expect_equal(cc2$min_cell, 6)
  ## both continuous: the two levels of the second moderator
  expect_equal(powerape:::aie3_cell_check(y, X, FALSE, FALSE)$min_cell, 24)
})

test_that("the engine runs the three-way estimator (replay of one study)", {
  d <- build3(pairs3$cb)
  bt <- powerape:::beta_true(d)
  n <- 900
  s1 <- powerape:::sim_ci(d, n, 1L, 0.95, seed = 5)
  rep1 <- powerape:::with_seed(5, {
    xx <- powerape:::draw_x(d, n)
    y <- rbinom(n, 1L, d$G(drop(xx$X %*% bt)))
    list(X = xx$X, y = y)
  })
  ft <- powerape:::fit_index_model(rep1$X, rep1$y, "probit", start = bt)
  e <- powerape:::aie3_est(ft$fit$coefficients, rep1$X, "probit", "normal", "binary")
  expect_equal(s1$est, e$ape, tolerance = 1e-12)
  expect_equal(s1$se, se3(e, powerape:::vcov_from_glmfit(ft$fit)), tolerance = 1e-12)
  ## linear analysis: the coefficient of the triple product, HC1
  sl <- powerape:::sim_ci(set_analysis(d, "linear"), n, 1L, 0.95, seed = 5)
  fo <- powerape:::fit_ols(rep1$X, rep1$y)
  expect_equal(sl$est, fo$coefficients[[8]], tolerance = 1e-12)
  expect_equal(sl$se, sqrt(powerape:::vcov_ols_hc1(fo, rep1$X)[8, 8]), tolerance = 1e-12)
})

test_that("simulated estimates and SEs agree with the closed form", {
  d <- bbb(target = 0.10)
  n <- 3000
  ## independent binary variables: cell shares are products of the margins
  g <- expand.grid(d = 0:1, m = 0:1, w = 0:1)
  share <- ifelse(g$d == 1, 0.5, 0.5) * ifelse(g$m == 1, 0.45, 0.55) *
    ifelse(g$w == 1, 0.4, 0.6)
  pc <- as.vector(d$cells3)
  se_cf <- sqrt(sum(pc * (1 - pc) / (n * share)))
  sim <- powerape:::sim_ci(d, n, 400L, 0.95, seed = 7)
  expect_true(all(sim$ok))
  expect_lt(abs(mean(sim$est) - 0.10), 4 * sd(sim$est) / sqrt(400))
  expect_equal(mean(sim$se), se_cf, tolerance = 0.02)
  expect_equal(sd(sim$est), se_cf, tolerance = 0.12)
  p <- suppressWarnings(ape_power(d, n, "detect", nsim = 400, seed = 7))
  p_cf <- pnorm(0.10 / se_cf - qnorm(0.975))
  expect_lt(abs(p$power - p_cf), 4 * sqrt(p_cf * (1 - p_cf) / 400) + 0.01)
  expect_identical(p$estimand, "aie3")
  expect_output(print(p), "assumed true three-way AIE +0.1000", fixed = TRUE)
  ## saturated design: the linear analysis gives the same estimates
  dl <- set_analysis(d, "linear")
  siml <- powerape:::sim_ci(dl, n, 50L, 0.95, seed = 7)
  expect_equal(siml$est, sim$est[1:50], tolerance = 1e-6)
  expect_equal(powerape:::linear_centre(dl), 0.10, tolerance = 1e-8)
  expect_output(print(dl), "coefficient of treat x female x public", fixed = TRUE)
})

test_that("the searches work on a three-way design", {
  d <- bbb(target = 0.20)
  g <- expand.grid(d = 0:1, m = 0:1, w = 0:1)
  share <- 0.5 * ifelse(g$m == 1, 0.45, 0.55) * ifelse(g$w == 1, 0.4, 0.6)
  pc <- as.vector(d$cells3)
  n_cf <- (qnorm(0.975) + qnorm(0.80))^2 * sum(pc * (1 - pc) / share) / 0.20^2
  r <- suppressWarnings(ape_n(d, power = 0.80, claim = "detect", nsim = 300, seed = 8))
  expect_equal(r$n, n_cf, tolerance = 0.10)
  expect_gt(r$power, 0.76)
  expect_output(print(r), "assumed true three-way AIE", fixed = TRUE)
  md <- suppressWarnings(ape_mde(d, n = round(n_cf), power = 0.80, claim = "detect",
                                 nsim = 300, seed = 9, confirm = FALSE))
  expect_equal(md$mde, 0.20, tolerance = 0.10)
  expect_identical(md$estimand, "aie3")
  ## the six other anchors are held, the design is re-pinned at the answer
  expect_equal(md$dgp$main3, d$main3)
  expect_equal(md$dgp$two_way3, d$two_way3)
  expect_equal(true_aie3(md$dgp), md$mde, tolerance = 1e-6)
  ## the short re-pin (only the triple product is solved again) equals a
  ## full pin at the same anchors
  full <- set_aie3(d, md$mde, main = d$main3, two_way = d$two_way3)
  expect_equal(powerape:::beta_true(md$dgp), powerape:::beta_true(full), tolerance = 1e-8)
  expect_equal(md$dgp$cells3, full$cells3, tolerance = 1e-8)
  dm <- build3(pairs3$cb)
  expect_equal(powerape:::beta_true(powerape:::aie3_repin(dm, 0.04)),
               powerape:::beta_true(set_aie3(dm, 0.04, main = dm$main3, two_way = dm$two_way3)),
               tolerance = 1e-8)
  expect_output(print(md), "minimum detectable three-way AIE", fixed = TRUE)
  ## an unpinned design has no anchors to hold
  d0 <- ape_dgp(focal = v3$fb, moderator = v3$mb, moderator2 = v3$w, baseline = 0.3)
  expect_error(ape_mde(d0, n = 2000, claim = "detect", nsim = 50), "pin the design first")
  cv <- suppressWarnings(ape_curve(d, n = c(1500, 3000), claim = "detect",
                                   nsim = 100, seed = 10))
  expect_true(cv$results$power[2] > cv$results$power[1])
  expect_output(print(cv), "three-way AIE", fixed = TRUE)
})

test_that("the design print lists anchors, and cell rates when all are binary", {
  out <- capture.output(print(bbb()))
  expect_true(any(grepl("second moderator: public", out, fixed = TRUE)))
  expect_true(any(grepl("true three-way AIE (treat x female x public) pinned at +0.0500",
                        out, fixed = TRUE)))
  expect_true(any(grepl("cell rates P(Y=1 | treat, female, public): 000 0.300 | 100 0.380",
                        out, fixed = TRUE)))
  out2 <- capture.output(print(build3(pairs3$cb)))
  expect_false(any(grepl("cell rates", out2, fixed = TRUE)))
  expect_true(any(grepl("two-way AIEs with the third variable at reference", out2,
                        fixed = TRUE)))
})
