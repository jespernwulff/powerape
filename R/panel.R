# Panel designs: correlated random effects probit/logit ------------------------
# DESIGN.md section 11. Model: y_it = 1[b0 + b_d d_it + z'gamma + c_i + u_it > 0],
# c_i = xi'(centered observed unit means of time-varying regressors) + a_i,
# a_i ~ N(0, var_a). Estimation: pooled probit/logit with Mundlak means and
# unit-clustered SEs. Estimand: ASF-based APE (means held at observed values).
# True values integrate a_i in closed form via scale_a = sqrt(1 + var_a).

# One draw of all regressor paths for n_units x T rows (unit-major order).
# Cross-regressor dependence R is applied at both the unit and within level.
panel_regressor_draw <- function(vars, iccs, R, n_units, T_p, seed = NULL) {
  k <- length(vars)
  with_seed(seed, {
    Rc <- chol(R)
    V <- matrix(rnorm(n_units * k), n_units) %*% Rc
    E <- matrix(rnorm(n_units * T_p * k), n_units * T_p) %*% Rc
    ui <- rep(seq_len(n_units), each = T_p)
    x <- vector("list", k)
    xbar <- vector("list", k)
    for (j in seq_len(k)) {
      u <- sqrt(iccs[j]) * V[ui, j] + sqrt(1 - iccs[j]) * E[, j]
      xj <- transform_u(u, vars[[j]])
      x[[j]] <- xj
      xbar[[j]] <- rowsum(xj, ui)[, 1] / T_p
    }
    list(x = x, xbar = xbar, unit = ui)
  })
}

# Deterministic integration pieces for calibration/inversion/truth:
# centered focal draws, q_it = z'gamma + m_i (heterogeneity's Mundlak part).
panel_integration <- function(dgp) {
  vars <- c(list(dgp$focal), dgp$covariates)
  rd <- panel_regressor_draw(vars, dgp$iccs, dgp$R, dgp$n_int, dgp$n_periods,
                             dgp$seed_int)
  idxz <- if (dgp$k > 0L) {
    Z <- do.call(cbind, rd$x[-1L])
    drop(Z %*% dgp$gamma)
  } else {
    rep(0, dgp$n_int * dgp$n_periods)
  }
  m_i <- if (length(dgp$tv_idx) && any(dgp$xi != 0)) {
    Mbar <- do.call(cbind, rd$xbar[dgp$tv_idx])
    drop(sweep(Mbar, 2, dgp$xbar_center) %*% dgp$xi)
  } else {
    rep(0, dgp$n_int)
  }
  list(xd_c = rd$x[[1L]] - dgp$focal_ref, q = idxz + m_i[rd$unit])
}

#' Specify a panel DGP: correlated random effects probit/logit
#'
#' Panel power analysis under the most general applied panel binary-response
#' model: a correlated random effects (CRE) specification estimated by
#' pooled probit/logit with Mundlak means and unit-clustered standard
#' errors (Wooldridge, 2010, ch. 15). Unobserved unit heterogeneity
#' `c_i = xi'(centered observed unit means of time-varying regressors) + a_i`
#' may correlate with the regressors through the Mundlak part; the target
#' effect remains the ASF-based average partial effect, which the pooled
#' CRE estimator recovers even though its coefficients are attenuated.
#'
#' Sample-size arguments of [ape_power()], [ape_curve()], [ape_n()], and
#' [ape_robust()] refer to the **number of units (clusters)** for panel
#' DGPs; each unit contributes `n_periods` observations.
#'
#' @inheritParams ape_dgp
#' @param focal,covariates [pa_var()] objects; use their `icc` argument to
#'   set within-unit persistence (default 0; `icc = 1` = time-constant).
#' @param n_periods Number of periods per unit (balanced; >= 2).
#' @param rho Latent unit-effect share `Var(c) / (Var(c) + 1)` — the
#'   xtprobit-style rho. 0 = no unobserved heterogeneity.
#' @param cre_share Fraction of `Var(c)` loaded (with equal standardized
#'   weights) on the centered observed unit means of the time-varying
#'   regressors; the remainder is the independent component `a_i`.
#'   `cre_share = 0` gives a pure random-effects world; the estimator
#'   includes Mundlak means either way.
#' @param n_int,seed_int Size (in **units**, each contributing `n_periods`
#'   rows) and seed of the deterministic integration draw used for
#'   calibration and inversion.
#'
#' @return An object of class `powerape_dgp` (route `"panel"`). Pin the
#'   effect with [set_ape()]; moderators (AIE) and the empirical/pilot
#'   routes are not yet available for panel DGPs.
#' @examples
#' \donttest{
#' d <- ape_dgp_panel(
#'   focal = pa_var("treat", "binary", p = 0.5, icc = 1),  # unit-level treatment
#'   covariates = list(pa_var("size", "normal", icc = 0.6)),
#'   n_periods = 4, rho = 0.3, cre_share = 0.5,
#'   baseline = 0.30, signal = 0.10
#' )
#' d <- set_ape(d, target = 0.08)
#' ape_power(d, n = 150, claim = "detect", nsim = 300, seed = 1)  # n = units
#' }
#' @export
ape_dgp_panel <- function(model = c("probit", "logit"), focal, covariates = list(),
                          n_periods, rho = 0, cre_share = 0, correlation = NULL,
                          baseline, signal = 0,
                          n_int = 1e5, seed_int = 20260814L) {
  model <- match.arg(model)
  if (inherits(covariates, "pa_var")) covariates <- list(covariates)
  stopifnot(inherits(focal, "pa_var"),
            all(vapply(covariates, inherits, logical(1), "pa_var")))
  T_p <- as.integer(n_periods)
  stopifnot(length(T_p) == 1L, T_p >= 2L)
  stopifnot(is.numeric(rho), length(rho) == 1L, rho >= 0, rho < 1)
  stopifnot(is.numeric(cre_share), length(cre_share) == 1L,
            cre_share >= 0, cre_share <= 1)
  stopifnot(is.numeric(baseline), length(baseline) == 1L,
            baseline > 0, baseline < 1)
  stopifnot(is.numeric(signal), length(signal) == 1L, signal >= 0, signal < 1)
  stopifnot(is.numeric(n_int), n_int >= 2e4)

  vars <- c(list(focal), covariates)
  k <- length(covariates)
  iccs <- vapply(vars, function(v) v$icc %||% 0, numeric(1))
  tv_idx <- which(iccs < 1)                # Mundlak means only for time-varying
  if (k == 0L && signal > 0)
    stop("`signal` > 0 requires at least one covariate.", call. = FALSE)
  if (cre_share > 0 && rho > 0 && length(tv_idx) == 0L)
    stop(paste("cre_share > 0 needs at least one time-varying regressor",
               "(icc < 1) for the heterogeneity to load on."), call. = FALSE)
  R <- build_corr(k + 1L, correlation)
  lf <- link_funs(model)

  rd <- panel_regressor_draw(vars, iccs, R, n_int, T_p, seed_int)

  ## nuisance coefficients from per-period signal (equal standardized weights)
  if (k > 0L) {
    Z <- do.call(cbind, rd$x[-1L])
    if (signal == 0) {
      gamma <- rep(0, k)
    } else {
      g_unit <- 1 / vapply(covariates, pa_var_sd, numeric(1))
      v_unit <- var(drop(Z %*% g_unit))
      gamma <- sqrt((signal / (1 - signal)) / v_unit) * g_unit
    }
    idxz <- drop(Z %*% gamma)
  } else {
    gamma <- numeric(0)
    idxz <- rep(0, n_int * T_p)
  }
  s2 <- if (k > 0L) var(idxz) else 0

  ## heterogeneity: Var(c) from rho; Mundlak loadings from cre_share
  var_c <- if (rho > 0) rho / (1 - rho) else 0
  if (var_c > 0 && cre_share > 0) {
    Mbar <- do.call(cbind, rd$xbar[tv_idx])
    xbar_center <- colMeans(Mbar)
    sds <- apply(Mbar, 2, sd)
    if (any(sds < 1e-10))
      stop("A time-varying regressor has (numerically) constant unit means.",
           call. = FALSE)
    w_unit <- 1 / sds
    m_unit <- drop(sweep(Mbar, 2, xbar_center) %*% w_unit)
    xi <- sqrt(cre_share * var_c / var(m_unit)) * w_unit
    m_i <- drop(sweep(Mbar, 2, xbar_center) %*% xi)
  } else {
    xbar_center <- if (length(tv_idx)) {
      colMeans(do.call(cbind, rd$xbar[tv_idx]))
    } else {
      numeric(0)
    }
    xi <- rep(0, length(tv_idx))
    m_i <- rep(0, n_int)
  }
  var_a <- (1 - cre_share) * var_c
  scale_a <- sqrt(1 + var_a)

  ## intercept from the counterfactual baseline (focal at reference),
  ## integrating a_i in closed form
  q <- idxz + m_i[rd$unit]
  f0 <- function(b0) mean(lf$G((b0 + q) / scale_a)) - baseline
  beta0 <- uniroot(f0, c(-15, 15) * scale_a, extendInt = "upX", tol = 1e-10)$root

  structure(list(
    model = model, link = model, G = lf$G, g = lf$g,
    focal = focal, moderator = NULL, covariates = covariates, k = k,
    R = R, R_chol = chol(R),
    baseline = baseline, signal = signal,
    signal_implied = s2 / (1 + s2), s2 = s2,
    gamma = gamma, beta0 = beta0,
    beta_focal = NULL, beta_mod = NULL, beta_int = NULL,
    estimand = NULL, target_est = NULL, p1 = NULL,
    main_focal = NULL, main_moderator = NULL,
    route = "panel", builder = "panel",
    n_periods = T_p, rho = rho, cre_share = cre_share,
    var_a = var_a, scale_a = scale_a,
    iccs = iccs, tv_idx = tv_idx, xi = xi, xbar_center = xbar_center,
    focal_ref = if (focal$type == "binary") 0 else focal$mean,
    mod_ref = NULL,
    n_int = as.integer(n_int), seed_int = seed_int,
    spec = list(model = model, focal = focal, covariates = covariates,
                n_periods = n_periods, rho = rho, cre_share = cre_share,
                correlation = correlation, baseline = baseline,
                signal = signal, n_int = n_int, seed_int = seed_int)
  ), class = "powerape_dgp")
}

## ---- inversion and truth (panel) --------------------------------------------

# Binary focal: solve mean(G((b0 + b + q)/s)) - p0 = target.
solve_binary_panel <- function(dgp, target, q) {
  s <- dgp$scale_a
  G <- dgp$G
  p0 <- mean(G((dgp$beta0 + q) / s))
  eps <- 1e-4
  if (target <= -p0 + eps || target >= (1 - p0) - eps)
    stop(sprintf(paste(
      "Target APE %.4f is infeasible for this panel DGP: with baseline",
      "P(Y=1 | reference) = %.3f the attainable APE range is (%.3f, %.3f)."),
      target, p0, -p0, 1 - p0), call. = FALSE)
  uniroot(function(b) mean(G((dgp$beta0 + b + q) / s)) - p0 - target,
          c(-8, 8) * s, extendInt = "upX", tol = 1e-9)$root
}

# Continuous focal: APE(b) = (b/s) mean(g((b0 + b*xd_c + q)/s)); plateaus.
solve_cont_panel <- function(dgp, target, xd_c, q) {
  s <- dgp$scale_a
  g <- dgp$g
  f <- function(b) (b / s) * mean(g((dgp$beta0 + b * xd_c + q) / s))
  sgn <- if (target < 0) -1 else 1
  sd_d <- pa_var_sd(dgp$focal)
  bgrid <- sgn * exp(seq(log(1e-3 / sd_d), log(200 * s / sd_d), length.out = 100))
  vals <- vapply(bgrid, f, numeric(1))
  amax <- max(abs(vals))
  if (abs(target) >= amax * (1 - 1e-3))
    stop(sprintf(paste(
      "Target APE %.4f is infeasible for this panel DGP: the attainable",
      "range is about (%.4f, %.4f). The average derivative plateaus as the",
      "coefficient grows."), target, -amax, amax), call. = FALSE)
  i <- which(abs(vals) >= abs(target))[1L]
  lo <- if (i == 1L) 0 else bgrid[i - 1L]
  uniroot(function(b) f(b) - target, sort(c(lo, bgrid[i])), tol = 1e-9)$root
}

true_ape_panel <- function(dgp) {
  id <- panel_integration(dgp)
  s <- dgp$scale_a
  if (dgp$focal$type == "binary") {
    mean(dgp$G((dgp$beta0 + dgp$beta_focal + id$q) / s) -
           dgp$G((dgp$beta0 + id$q) / s))
  } else {
    (dgp$beta_focal / s) *
      mean(dgp$g((dgp$beta0 + dgp$beta_focal * id$xd_c + id$q) / s))
  }
}

## ---- engine draw (panel) -----------------------------------------------------

# One simulated panel of n units x T rows. Returns the design matrix with
# Mundlak means (raw, uncentered) for time-varying regressors, the cluster
# id, the true outcome probabilities (structural index + drawn a_i), and
# pseudo-true start values (attenuated coefficients) for the pooled fit.
draw_x_panel <- function(dgp, n_units) {
  if (is.null(dgp$beta_focal))
    stop("This DGP has no target effect yet - call set_ape() first.", call. = FALSE)
  vars <- c(list(dgp$focal), dgp$covariates)
  T_p <- dgp$n_periods
  rd <- panel_regressor_draw(vars, dgp$iccs, dgp$R, n_units, T_p)
  d <- rd$x[[1L]]
  Z <- if (dgp$k > 0L) do.call(cbind, rd$x[-1L]) else NULL
  Mbar_cols <- if (length(dgp$tv_idx)) {
    Mb <- do.call(cbind, rd$xbar[dgp$tv_idx])
    Mb[rd$unit, , drop = FALSE]
  } else NULL
  X <- cbind(1, d, Z, Mbar_cols)

  ## structural probabilities
  idxz <- if (dgp$k > 0L) drop(Z %*% dgp$gamma) else 0
  m_i <- if (length(dgp$tv_idx) && any(dgp$xi != 0)) {
    drop(sweep(do.call(cbind, rd$xbar[dgp$tv_idx]), 2, dgp$xbar_center) %*% dgp$xi)
  } else {
    rep(0, n_units)
  }
  a_i <- if (dgp$var_a > 0) rnorm(n_units, 0, sqrt(dgp$var_a)) else rep(0, n_units)
  eta <- dgp$beta0 + dgp$beta_focal * d + idxz + (m_i + a_i)[rd$unit]
  pr <- dgp$G(eta)

  ## pseudo-true (attenuated) coefficients as IRLS start values
  s <- dgp$scale_a
  start <- c(dgp$beta0 - sum(dgp$xi * dgp$xbar_center),
             dgp$beta_focal, dgp$gamma, dgp$xi) / s

  list(X = X, focal_col = 2L, id = rd$unit, pr = pr, start = start)
}
