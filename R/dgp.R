# DGP construction and calibration --------------------------------------------

link_funs <- function(model) {
  switch(model,
         probit = list(G = pnorm, g = dnorm, Ginv = qnorm),
         logit  = list(G = plogis, g = dlogis, Ginv = qlogis),
         stop("Unknown model: ", model, call. = FALSE))
}

# Correlation matrix over (focal[, moderator], covariates).
build_corr <- function(K, correlation) {
  R <- if (is.null(correlation)) {
    diag(K)
  } else if (is.matrix(correlation)) {
    correlation
  } else {
    stopifnot(is.numeric(correlation), length(correlation) == 1L)
    m <- matrix(correlation, K, K)
    diag(m) <- 1
    m
  }
  if (!all(dim(R) == c(K, K)))
    stop(sprintf("`correlation` must be a %d x %d matrix over (focal%s, covariates).",
                 K, K, if (K > 1) ", moderator" else ""), call. = FALSE)
  if (max(abs(R - t(R))) > 1e-8 || max(abs(diag(R) - 1)) > 1e-8)
    stop("`correlation` must be symmetric with unit diagonal.", call. = FALSE)
  ok <- tryCatch({ chol(R); TRUE }, error = function(e) FALSE)
  if (!ok) stop("`correlation` is not positive definite.", call. = FALSE)
  R
}

# n x K matrix of standard-normal rows with correlation R (R_chol = chol(R)).
draw_std_mvn <- function(n, R_chol) {
  matrix(rnorm(n * ncol(R_chol)), n) %*% R_chol
}

# Deterministic draw of ALL regressors (focal, then moderator if any, then
# covariates) from the Gaussian copula; reproduced exactly from seed_int.
reg_draw_parametric <- function(focal, moderator, covariates, R, n_int, seed_int) {
  k <- length(covariates)
  has_m <- !is.null(moderator)
  with_seed(seed_int, {
    U <- draw_std_mvn(n_int, chol(R))
    xd <- transform_u(U[, 1], focal)
    md <- if (has_m) transform_u(U[, 2], moderator) else NULL
    z_off <- 1L + has_m
    Z <- if (k > 0L) {
      Zm <- U[, z_off + seq_len(k), drop = FALSE]
      for (j in seq_len(k)) Zm[, j] <- transform_u(Zm[, j], covariates[[j]])
      Zm
    } else {
      matrix(numeric(0), n_int, 0L)
    }
    list(xd = xd, md = md, Z = Z)
  })
}

# Integration set for calibration, inversion, and truth evaluation: centered
# focal draws, centered moderator draws (if any), and the nuisance index
# z'gamma, under the route's covariate distribution. Deterministic.
integration_draw <- function(dgp) {
  if (identical(dgp$route, "iv")) return(iv_integration(dgp))
  if (dgp$route == "empirical") {
    idxz <- if (dgp$k > 0L) drop(dgp$emp_x %*% dgp$gamma) else rep(0, dgp$n_emp)
    if (!is.null(dgp$emp_xd)) {
      xd_c <- dgp$emp_xd - dgp$focal_ref
    } else if (dgp$focal$type != "binary") {
      # Independent parametric focal over empirical covariates: Monte Carlo
      # over the product measure, tiling the pilot rows for precision.
      m <- max(1L, ceiling(dgp$n_int / length(idxz)))
      idxz <- rep(idxz, m)
      xd_c <- with_seed(dgp$seed_int,
                        transform_u(rnorm(length(idxz)), dgp$focal)) - dgp$focal_ref
    } else {
      xd_c <- NULL
    }
    list(xd_c = xd_c, md_c = NULL, idxz = idxz)
  } else {
    rd <- reg_draw_parametric(dgp$focal, dgp$moderator, dgp$covariates,
                              dgp$R, dgp$n_int, dgp$seed_int)
    idxz <- if (dgp$k > 0L) drop(rd$Z %*% dgp$gamma) else rep(0, dgp$n_int)
    list(xd_c = rd$xd - dgp$focal_ref,
         md_c = if (is.null(rd$md)) NULL else rd$md - dgp$mod_ref,
         idxz = idxz)
  }
}

#' Specify the data-generating process for an APE/AIE power analysis
#'
#' Describes the world the study will sample from: an index model
#' (probit/logit), a binary or continuous (normal) focal variable, an
#' optional moderator (for average-interaction-effect designs), covariates
#' with a Gaussian-copula dependence structure, and an outcome anchored on
#' the friendly scale -- the baseline outcome rate plus the strength of the
#' nuisance covariates. Effect coefficients are set later, in effect-size
#' units, by [set_ape()] (no moderator) or [set_aie()] (with moderator).
#'
#' @param model `"probit"` (default) or `"logit"`.
#' @param focal A [pa_var()], `"binary"` or `"normal"`. For a continuous
#'   focal the APE is the average derivative, and reference = its mean.
#' @param moderator Optional [pa_var()]. When present, the index includes
#'   the moderator and a focal-by-moderator interaction, the estimand is
#'   the average interaction effect (AIE, as in the `ginteff` package), and
#'   the DGP must be pinned with [set_aie()].
#' @param covariates List of [pa_var()] objects (may be empty).
#' @param correlation `NULL` (independence), a scalar exchangeable
#'   correlation, or a full correlation matrix over
#'   `(focal, moderator, covariates)` -- Gaussian-copula dependence.
#' @param baseline Baseline outcome rate in (0, 1): the average probability
#'   `P(Y = 1)` with the focal (and moderator, if any) set to their
#'   reference values (0 for binary, the mean for continuous), covariates
#'   at their population distribution. The intercept is calibrated to
#'   reproduce it.
#' @param signal Nuisance-signal strength: the McKelvey-Zavoina-style latent
#'   pseudo-R-squared of the nuisance index, i.e. `Var(z'gamma) /
#'   (Var(z'gamma) + 1)`. Equal standardized weights, scaled to hit
#'   `signal`. Ignored when `gamma` is supplied.
#' @param gamma Optional explicit nuisance coefficients (length =
#'   `length(covariates)`); overrides `signal`.
#' @param n_int,seed_int Size and seed of the deterministic Monte Carlo
#'   integration draw used for calibration and inversion.
#'
#' @return An object of class `powerape_dgp`. Pin the effect with
#'   [set_ape()] / [set_aie()], then use [ape_power()], [ape_curve()],
#'   [ape_n()], or [ape_robust()].
#' @examples
#' d <- ape_dgp(
#'   model = "probit",
#'   focal = pa_var("treat", "binary", p = 0.5),
#'   covariates = list(pa_var("age", "normal", mean = 45, sd = 12)),
#'   baseline = 0.30, signal = 0.15
#' )
#' d
#' @export
ape_dgp <- function(model = c("probit", "logit"), focal, moderator = NULL,
                    covariates = list(), correlation = NULL, baseline,
                    signal = 0, gamma = NULL,
                    n_int = 4e5, seed_int = 20260812L) {
  model <- match.arg(model)
  if (inherits(covariates, "pa_var")) covariates <- list(covariates)
  stopifnot(inherits(focal, "pa_var"),
            all(vapply(covariates, inherits, logical(1), "pa_var")))
  if (!is.null(moderator)) {
    stopifnot(inherits(moderator, "pa_var"))
    if (identical(moderator$name, focal$name))
      stop("`moderator` must be a different variable than `focal`.", call. = FALSE)
  }
  stopifnot(is.numeric(baseline), length(baseline) == 1L,
            baseline > 0, baseline < 1)
  stopifnot(is.numeric(n_int), n_int >= 1e4)

  k <- length(covariates)
  has_m <- !is.null(moderator)
  K <- 1L + has_m + k
  R <- build_corr(K, correlation)
  lf <- link_funs(model)

  if (is.null(gamma)) {
    stopifnot(is.numeric(signal), length(signal) == 1L, signal >= 0, signal < 1)
    if (k == 0L && signal > 0)
      stop("`signal` > 0 requires at least one covariate.", call. = FALSE)
  } else {
    if (!is.numeric(gamma) || length(gamma) != k)
      stop(sprintf("`gamma` must be a numeric vector of length %d.", k), call. = FALSE)
  }

  if (k == 0L) {
    gamma_out <- numeric(0)
    s2 <- 0
    signal_implied <- 0
    beta0 <- lf$Ginv(baseline)
  } else {
    rd <- reg_draw_parametric(focal, moderator, covariates, R, n_int, seed_int)
    Z <- rd$Z
    if (is.null(gamma)) {
      if (signal == 0) {
        gamma_out <- rep(0, k)
      } else {
        g_unit <- 1 / vapply(covariates, pa_var_sd, numeric(1))
        v_unit <- var(drop(Z %*% g_unit))
        gamma_out <- sqrt((signal / (1 - signal)) / v_unit) * g_unit
      }
    } else {
      gamma_out <- gamma
    }
    idx <- drop(Z %*% gamma_out)
    s2 <- var(idx)
    signal_implied <- s2 / (1 + s2)
    f0 <- function(b0) mean(lf$G(b0 + idx)) - baseline
    beta0 <- uniroot(f0, c(-15, 15), extendInt = "upX", tol = 1e-10)$root
  }

  structure(list(
    model = model, link = model, G = lf$G, g = lf$g,
    focal = focal, moderator = moderator, covariates = covariates, k = k,
    R = R, R_chol = chol(R),
    baseline = baseline,
    signal = if (is.null(gamma)) signal else NA_real_,
    signal_implied = signal_implied, s2 = s2,
    gamma = gamma_out, beta0 = beta0,
    beta_focal = NULL, beta_mod = NULL, beta_int = NULL,
    estimand = NULL, target_est = NULL, p1 = NULL,
    main_focal = NULL, main_moderator = NULL,
    route = "parametric", builder = "parametric",
    focal_ref = if (focal$type == "binary") 0 else focal$mean,
    mod_ref = if (!has_m) NULL else if (moderator$type == "binary") 0 else moderator$mean,
    n_int = n_int, seed_int = seed_int,
    spec = list(model = model, focal = focal, moderator = moderator,
                covariates = covariates, correlation = correlation,
                baseline = baseline, signal = signal, gamma = gamma,
                n_int = n_int, seed_int = seed_int)
  ), class = "powerape_dgp")
}

# One simulated design matrix of size n. Layout:
#   no moderator:  [1, focal, Z]           focal_col = 2
#   moderator:     [1, focal, mod, focal*mod, Z]   cols 2, 3, 4
draw_x <- function(dgp, n) {
  if (identical(dgp$route, "panel")) return(draw_x_panel(dgp, n))
  if (identical(dgp$route, "iv")) return(draw_x_iv(dgp, n))
  if (dgp$route == "empirical") return(draw_x_empirical(dgp, n))
  U <- draw_std_mvn(n, dgp$R_chol)
  d <- transform_u(U[, 1], dgp$focal)
  has_m <- !is.null(dgp$moderator)
  m <- if (has_m) transform_u(U[, 2], dgp$moderator) else NULL
  z_off <- 1L + has_m
  Z <- if (dgp$k > 0L) {
    Zm <- U[, z_off + seq_len(dgp$k), drop = FALSE]
    for (j in seq_len(dgp$k)) Zm[, j] <- transform_u(Zm[, j], dgp$covariates[[j]])
    Zm
  } else NULL
  X <- if (has_m) cbind(1, d, m, d * m, Z) else cbind(1, d, Z)
  list(X = X, focal_col = 2L)
}

# True coefficient vector in design-matrix (raw, uncentered) order. The DGP
# is calibrated with focal (and moderator) centered at their references, so
# the raw coefficients absorb the centering.
beta_true <- function(dgp) {
  if (is.null(dgp$beta_focal))
    stop("This DGP has no target effect yet - call set_ape() or set_aie() first.",
         call. = FALSE)
  if (is.null(dgp$moderator)) {
    c(dgp$beta0 - dgp$beta_focal * dgp$focal_ref, dgp$beta_focal, dgp$gamma)
  } else {
    bf <- dgp$beta_focal; bm <- dgp$beta_mod; bi <- dgp$beta_int
    rd <- dgp$focal_ref; rm_ <- dgp$mod_ref
    c(dgp$beta0 - bf * rd - bm * rm_ + bi * rd * rm_,
      bf - bi * rm_,
      bm - bi * rd,
      bi,
      dgp$gamma)
  }
}
