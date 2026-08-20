# Panel designs: correlated random effects probit -----------------------------
# DESIGN.md section 11. Model: y_it = 1[b0 + b_d d_it + z'gamma + c_i + u_it > 0],
# c_i = xi'(centered observed unit means of time-varying regressors) + a_i,
# a_i ~ N(0, var_a). Estimation: pooled probit with Mundlak means and
# unit-clustered SEs. Estimand: ASF-based APE (means held at observed values).
# True values integrate a_i exactly: probit in closed form via
# scale_a = sqrt(1 + var_a).
# Probit-only by design (JW, 2026-08-20): the normal-normal convolution that
# makes the pooled CRE estimator recover the ASF has no logistic counterpart,
# so a pooled CRE logit is a QMLE approximation; ape_dgp_panel refuses
# model = "logit" with a teaching error. The Gauss-Hermite machinery in
# panel_mix_funs is retained for the general-link mixture (and gauss_hermite
# remains unit-tested), but no public route reaches its logit branch.

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
# centered focal (and moderator) draws, q_it = z'gamma + m_i (the
# heterogeneity's Mundlak part). The design variables' structural terms are
# excluded -- their coefficients are what inversion solves for.
panel_integration <- function(dgp) {
  has_m <- !is.null(dgp$moderator)
  vars <- c(list(dgp$focal), if (has_m) list(dgp$moderator), dgp$covariates)
  rd <- panel_regressor_draw(vars, dgp$iccs, dgp$R, dgp$n_int, dgp$n_periods,
                             dgp$seed_int)
  idxz <- if (dgp$k > 0L) {
    Z <- do.call(cbind, rd$x[-seq_len(1L + has_m)])
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
  list(xd_c = rd$x[[1L]] - dgp$focal_ref,
       xm_c = if (has_m) rd$x[[2L]] - dgp$mod_ref else NULL,
       q = idxz + m_i[rd$unit])
}

# Gauss-Hermite nodes/weights (physicists' convention) via Golub-Welsch:
# integral f(t) exp(-t^2) dt = sum w_j f(t_j).
gauss_hermite <- function(n) {
  J <- matrix(0, n, n)
  b <- sqrt(seq_len(n - 1L) / 2)
  J[cbind(seq_len(n - 1L), seq_len(n - 1L) + 1L)] <- b
  J[cbind(seq_len(n - 1L) + 1L, seq_len(n - 1L))] <- b
  e <- eigen(J, symmetric = TRUE)
  list(nodes = e$values, weights = sqrt(pi) * e$vectors[1L, ]^2)
}

# E_a[G(x + a)] with a ~ N(0, var_a), and its derivative in x. Probit:
# exact closed form (normal convolution). Logit: 20-node Gauss-Hermite --
# the probit-style index/scale division has no logistic analog.
panel_mix_funs <- function(model, var_a) {
  lf <- link_funs(model)
  if (var_a == 0) return(list(P = lf$G, p = lf$g))
  if (model == "probit") {
    s <- sqrt(1 + var_a)
    return(list(P = function(x) pnorm(x / s),
                p = function(x) dnorm(x / s) / s))
  }
  gh <- gauss_hermite(20L)
  aj <- sqrt(2 * var_a) * gh$nodes
  wj <- gh$weights / sqrt(pi)
  mix <- function(f) function(x) {
    out <- 0
    for (j in seq_along(aj)) out <- out + wj[j] * f(x + aj[j])
    out
  }
  list(P = mix(lf$G), p = mix(lf$g))
}

#' Specify a panel DGP: correlated random effects probit
#'
#' Panel power analysis under the most general applied panel binary-response
#' model: a correlated random effects (CRE) specification estimated by
#' pooled probit with Mundlak means and unit-clustered standard
#' errors (Wooldridge, 2010, ch. 15). Unobserved unit heterogeneity
#' `c_i = xi'(centered observed unit means of time-varying regressors) + a_i`
#' may correlate with the regressors through the Mundlak part; the target
#' effect remains the ASF-based average partial effect, which the pooled
#' CRE estimator recovers even though its coefficients are attenuated.
#'
#' The panel route is **probit-only**. The exactness above rests on the
#' normal-normal convolution (a normal unit effect rescales the probit
#' index by `sqrt(1 + var_a)`); the logistic link has no such stability,
#' so a pooled CRE logit estimates only a quasi-maximum-likelihood
#' approximation of the ASF. Rather than price designs with an
#' approximate estimand, `ape_dgp_panel(model = "logit")` is refused with
#' an error explaining this. On the probability scale a probit world
#' calibrated to the same baseline and effect is practically
#' indistinguishable from the logit world it replaces.
#'
#' Sample-size arguments of [ape_power()], [ape_curve()], [ape_n()], and
#' [ape_robust()] refer to the **number of units (clusters)** for panel
#' DGPs; each unit contributes `n_periods` observations.
#'
#' @inheritParams ape_dgp
#' @param model `"probit"` (the only panel link; see Details for why a
#'   `"logit"` request is refused with a teaching error).
#' @param focal,covariates [pa_var()] objects; use their `icc` argument to
#'   set within-unit persistence (default 0; `icc = 1` = time-constant).
#' @param moderator Optional binary [pa_var()] for panel AIE designs: the
#'   index gains the moderator and a focal-by-moderator interaction, the
#'   estimand becomes the AIE (pin it with [set_aie()]), and the estimating
#'   model gains the interaction's own Mundlak mean whenever the product is
#'   time-varying. Currently binary focal x binary moderator only.
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
#'   effect with [set_ape()] (no moderator) or [set_aie()] (with
#'   moderator); the empirical/pilot routes are not yet available for
#'   panel DGPs.
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
ape_dgp_panel <- function(model = "probit", focal, moderator = NULL,
                          covariates = list(),
                          n_periods, rho = 0, cre_share = 0, correlation = NULL,
                          baseline, signal = 0,
                          n_int = 1e5, seed_int = 20260814L) {
  model <- as.character(model)[1L]
  if (identical(model, "logit"))
    stop(paste0(
      "Panel designs are probit-only.\n",
      "The device that lets a pooled CRE estimator recover the average\n",
      "structural function is specific to the probit link: a normal unit\n",
      "effect convolved with a normal error just rescales the index by\n",
      "sqrt(1 + var_a). The logistic link has no such stability, so a pooled\n",
      "CRE logit estimates only a quasi-ML approximation of the ASF, and\n",
      "powerape does not price designs with an approximate estimand.\n",
      "Use ape_dgp_panel(model = \"probit\", ...): the probit panel route is\n",
      "validated exactly (battery V8, V12, V16), and on the probability\n",
      "scale the two links describe nearly identical worlds anyway."),
      call. = FALSE)
  if (!identical(model, "probit"))
    stop("`model` must be \"probit\" for panel DGPs.", call. = FALSE)
  if (inherits(covariates, "pa_var")) covariates <- list(covariates)
  stopifnot(inherits(focal, "pa_var"),
            all(vapply(covariates, inherits, logical(1), "pa_var")))
  has_m <- !is.null(moderator)
  if (has_m) {
    stopifnot(inherits(moderator, "pa_var"))
    if (identical(moderator$name, focal$name))
      stop("`moderator` must be a different variable than `focal`.", call. = FALSE)
    if (focal$type != "binary" || moderator$type != "binary")
      stop(paste("Panel AIE designs currently support a binary focal x binary",
                 "moderator; continuous pairs are on the roadmap."),
           call. = FALSE)
  }
  T_p <- as.integer(n_periods)
  stopifnot(length(T_p) == 1L, T_p >= 2L)
  stopifnot(is.numeric(rho), length(rho) == 1L, rho >= 0, rho < 1)
  stopifnot(is.numeric(cre_share), length(cre_share) == 1L,
            cre_share >= 0, cre_share <= 1)
  stopifnot(is.numeric(baseline), length(baseline) == 1L,
            baseline > 0, baseline < 1)
  stopifnot(is.numeric(signal), length(signal) == 1L, signal >= 0, signal < 1)
  stopifnot(is.numeric(n_int), n_int >= 2e4)

  vars <- c(list(focal), if (has_m) list(moderator), covariates)
  k <- length(covariates)
  iccs <- vapply(vars, function(v) v$icc %||% 0, numeric(1))
  tv_idx <- which(iccs < 1)                # Mundlak means only for time-varying
  if (k == 0L && signal > 0)
    stop("`signal` > 0 requires at least one covariate.", call. = FALSE)
  if (cre_share > 0 && rho > 0 && length(tv_idx) == 0L)
    stop(paste("cre_share > 0 needs at least one time-varying regressor",
               "(icc < 1) for the heterogeneity to load on."), call. = FALSE)
  R <- build_corr(k + 1L + has_m, correlation)
  lf <- link_funs(model)

  rd <- panel_regressor_draw(vars, iccs, R, n_int, T_p, seed_int)

  ## nuisance coefficients from per-period signal (equal standardized weights)
  if (k > 0L) {
    Z <- do.call(cbind, rd$x[-seq_len(1L + has_m)])
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
  mx <- panel_mix_funs(model, var_a)

  ## intercept from the counterfactual baseline (focal at reference),
  ## integrating a_i exactly (closed form / quadrature)
  q <- idxz + m_i[rd$unit]
  f0 <- function(b0) mean(mx$P(b0 + q)) - baseline
  beta0 <- uniroot(f0, c(-15, 15) * scale_a, extendInt = "upX", tol = 1e-10)$root

  structure(list(
    model = model, link = model, G = lf$G, g = lf$g,
    focal = focal, moderator = moderator, covariates = covariates, k = k,
    R = R, R_chol = chol(R),
    baseline = baseline, signal = signal,
    signal_implied = s2 / (1 + s2), s2 = s2,
    gamma = gamma, beta0 = beta0,
    beta_focal = NULL, beta_mod = NULL, beta_int = NULL,
    estimand = NULL, target_est = NULL, p1 = NULL,
    main_focal = NULL, main_moderator = NULL,
    route = "panel", builder = "panel",
    n_periods = T_p, rho = rho, cre_share = cre_share,
    var_a = var_a, scale_a = scale_a, mix_P = mx$P, mix_p = mx$p,
    iccs = iccs, tv_idx = tv_idx, xi = xi, xbar_center = xbar_center,
    focal_ref = if (focal$type == "binary") 0 else focal$mean,
    mod_ref = if (has_m) 0 else NULL,
    n_int = as.integer(n_int), seed_int = seed_int,
    spec = list(model = model, focal = focal, moderator = moderator,
                covariates = covariates,
                n_periods = n_periods, rho = rho, cre_share = cre_share,
                correlation = correlation, baseline = baseline,
                signal = signal, n_int = n_int, seed_int = seed_int)
  ), class = "powerape_dgp")
}

## ---- inversion and truth (panel) --------------------------------------------

# Binary focal: solve E_a mean(G(b0 + b + q + a)) - p0 = target.
solve_binary_panel <- function(dgp, target, q, label = "APE") {
  P <- dgp$mix_P
  p0 <- mean(P(dgp$beta0 + q))
  eps <- 1e-4
  if (target <= -p0 + eps || target >= (1 - p0) - eps)
    stop(sprintf(paste(
      "Target %s %.4f is infeasible for this panel DGP: with baseline",
      "P(Y=1 | reference) = %.3f the attainable range is (%.3f, %.3f)."),
      label, target, p0, -p0, 1 - p0), call. = FALSE)
  uniroot(function(b) mean(P(dgp$beta0 + b + q)) - p0 - target,
          c(-8, 8) * dgp$scale_a, extendInt = "upX", tol = 1e-9)$root
}

# Continuous focal: APE(b) = b E_a mean(g(b0 + b*xd_c + q + a)); plateaus.
solve_cont_panel <- function(dgp, target, xd_c, q) {
  s <- dgp$scale_a
  f <- function(b) b * mean(dgp$mix_p(dgp$beta0 + b * xd_c + q))
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
  if (dgp$focal$type == "binary") {
    mean(dgp$mix_P(dgp$beta0 + dgp$beta_focal + id$q) -
           dgp$mix_P(dgp$beta0 + id$q))
  } else {
    dgp$beta_focal *
      mean(dgp$mix_p(dgp$beta0 + dgp$beta_focal * id$xd_c + id$q))
  }
}

# True panel AIE (binary x binary): ASF double difference over the four
# counterfactual cells, Mundlak means held fixed, a_i integrated exactly.
true_aie_panel <- function(dgp) {
  id <- panel_integration(dgp)
  P <- dgp$mix_P
  b0 <- dgp$beta0
  bf <- dgp$beta_focal
  bm <- dgp$beta_mod
  bi <- dgp$beta_int
  mean(P(b0 + bf + bm + bi + id$q)) - mean(P(b0 + bm + id$q)) -
    mean(P(b0 + bf + id$q)) + mean(P(b0 + id$q))
}

## ---- engine draw (panel) -----------------------------------------------------

# One simulated panel of n units x T rows. Returns the design matrix with
# Mundlak means (raw, uncentered) for time-varying regressors -- plus, for
# moderated designs, the interaction's own mean whenever the product is
# time-varying (the mean of a product is not the product of means) -- the
# cluster id, the true outcome probabilities (structural index + drawn
# a_i), and pseudo-true start values (attenuated coefficients).
draw_x_panel <- function(dgp, n_units) {
  if (is.null(dgp$beta_focal))
    stop("This DGP has no target effect yet - call set_ape() (or set_aie()) first.",
         call. = FALSE)
  has_m <- !is.null(dgp$moderator)
  vars <- c(list(dgp$focal), if (has_m) list(dgp$moderator), dgp$covariates)
  T_p <- dgp$n_periods
  rd <- panel_regressor_draw(vars, dgp$iccs, dgp$R, n_units, T_p)
  d <- rd$x[[1L]]
  m <- if (has_m) rd$x[[2L]] else NULL
  Z <- if (dgp$k > 0L) do.call(cbind, rd$x[-seq_len(1L + has_m)]) else NULL
  Mbar_cols <- if (length(dgp$tv_idx)) {
    Mb <- do.call(cbind, rd$xbar[dgp$tv_idx])
    Mb[rd$unit, , drop = FALSE]
  } else NULL
  if (has_m) {
    dm <- d * m
    dmbar_col <- if (dgp$iccs[1L] < 1 || dgp$iccs[2L] < 1) {
      (rowsum(dm, rd$unit)[, 1L] / T_p)[rd$unit]
    } else NULL
    X <- cbind(1, d, m, dm, Z, Mbar_cols, dmbar_col)
  } else {
    X <- cbind(1, d, Z, Mbar_cols)
  }

  ## structural probabilities
  idxz <- if (dgp$k > 0L) drop(Z %*% dgp$gamma) else 0
  m_i <- if (length(dgp$tv_idx) && any(dgp$xi != 0)) {
    drop(sweep(do.call(cbind, rd$xbar[dgp$tv_idx]), 2, dgp$xbar_center) %*% dgp$xi)
  } else {
    rep(0, n_units)
  }
  a_i <- if (dgp$var_a > 0) rnorm(n_units, 0, sqrt(dgp$var_a)) else rep(0, n_units)
  eta <- dgp$beta0 + dgp$beta_focal * d + idxz + (m_i + a_i)[rd$unit]
  if (has_m) eta <- eta + dgp$beta_mod * m + dgp$beta_int * dm
  pr <- dgp$G(eta)

  ## pseudo-true (attenuated) coefficients as IRLS start values; the
  ## interaction-mean column's structural coefficient is 0
  s <- dgp$scale_a
  start <- c(dgp$beta0 - sum(dgp$xi * dgp$xbar_center),
             dgp$beta_focal,
             if (has_m) c(dgp$beta_mod, dgp$beta_int),
             dgp$gamma, dgp$xi,
             if (has_m && (dgp$iccs[1L] < 1 || dgp$iccs[2L] < 1)) 0) / s

  list(X = X, focal_col = 2L, id = rd$unit, pr = pr, start = start)
}
