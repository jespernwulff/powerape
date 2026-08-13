# Empirical-covariate route: resample pilot data ------------------------------

#' Specify a DGP with covariates resampled from pilot data
#'
#' Instead of parametric marginals and a copula, the covariate distribution
#' is the empirical joint distribution of a pilot (or register) dataset:
#' simulation resamples its rows with replacement, and all calibration
#' integrals are exact averages over the rows. This preserves the covariate
#' dependence structure for free and is the recommended route when data
#' exist (DESIGN.md section 3.1).
#'
#' The focal variable comes in two flavors:
#' * a **column name** in `data`: the focal variable is resampled jointly
#'   with the covariates, preserving its empirical dependence with them
#'   (observational mental model). Type is inferred (binary if the column is
#'   0/1, continuous otherwise; a continuous focal's reference value is its
#'   mean).
#' * a **[pa_var()]**: the focal variable is drawn independently of the
#'   resampled rows (randomized-treatment mental model); all of `data` is
#'   treated as covariates.
#'
#' @param model `"probit"` (default) or `"logit"`.
#' @param data Data frame or matrix of numeric columns, no missing values.
#'   Dummy-encode factors first.
#' @param focal Column name in `data`, or a [pa_var()].
#' @param baseline Baseline outcome rate in (0, 1): the average `P(Y = 1)`
#'   with the focal at its reference value (0 for binary, the mean for
#'   continuous), covariates at their empirical distribution.
#' @param signal Nuisance-signal strength (latent pseudo-R-squared of the
#'   nuisance index), as in [ape_dgp()]; equal weights on the standardized
#'   scale, scaled on the empirical rows. Ignored when `gamma` is supplied.
#' @param gamma Optional explicit nuisance coefficients, one per covariate
#'   column.
#' @param seed_int Seed for the deterministic draw used only when a
#'   continuous `pa_var` focal is combined with empirical covariates.
#'
#' @return An object of class `powerape_dgp` (route `"empirical"`).
#' @examples
#' pilot <- data.frame(treat = rbinom(600, 1, 0.5),
#'                     age = rnorm(600, 45, 12))
#' d <- ape_dgp_empirical(data = pilot, focal = "treat",
#'                        baseline = 0.30, signal = 0.15)
#' d
#' @export
ape_dgp_empirical <- function(model = c("probit", "logit"), data, focal,
                              baseline, signal = 0, gamma = NULL,
                              seed_int = 20260812L) {
  model <- match.arg(model)
  lf <- link_funs(model)
  stopifnot(is.numeric(baseline), length(baseline) == 1L,
            baseline > 0, baseline < 1)

  dat <- as.data.frame(data)
  if (!nrow(dat)) stop("`data` has no rows.", call. = FALSE)
  num_ok <- vapply(dat, is.numeric, logical(1))
  if (!all(num_ok))
    stop("All columns of `data` must be numeric (dummy-encode factors first). Offending: ",
         paste(names(dat)[!num_ok], collapse = ", "), call. = FALSE)
  M <- as.matrix(dat)
  if (!all(is.finite(M)))
    stop("`data` contains NA/NaN/Inf values; supply complete cases.", call. = FALSE)
  if (nrow(M) < 30L)
    stop("The empirical route needs at least 30 pilot rows.", call. = FALSE)
  if (nrow(M) < 500L)
    warning("Fewer than 500 pilot rows: pilot sampling noise becomes part of the design (cf. Hancock & Feng, 2025).")

  if (is.character(focal)) {
    stopifnot(length(focal) == 1L)
    if (!focal %in% colnames(M))
      stop(sprintf("Focal column \"%s\" not found in `data`.", focal), call. = FALSE)
    xd <- M[, focal]
    covs <- M[, setdiff(colnames(M), focal), drop = FALSE]
    if (length(unique(xd)) < 2L)
      stop("The focal column is constant.", call. = FALSE)
    focal_var <- if (all(xd %in% c(0, 1))) {
      pa_var(focal, "binary", p = mean(xd))
    } else {
      pa_var(focal, "normal", mean = mean(xd), sd = sd(xd))
    }
    emp_xd <- xd
  } else if (inherits(focal, "pa_var")) {
    focal_var <- focal
    covs <- M
    emp_xd <- NULL
  } else {
    stop("`focal` must be a column name in `data` or a pa_var().", call. = FALSE)
  }

  k <- ncol(covs)
  if (is.null(gamma)) {
    stopifnot(is.numeric(signal), length(signal) == 1L, signal >= 0, signal < 1)
    if (k == 0L && signal > 0)
      stop("`signal` > 0 requires at least one covariate column.", call. = FALSE)
    if (k == 0L || signal == 0) {
      gamma_out <- rep(0, k)
    } else {
      sds <- apply(covs, 2, sd)
      if (any(sds < 1e-12))
        stop("Constant covariate column(s): ",
             paste(colnames(covs)[sds < 1e-12], collapse = ", "), call. = FALSE)
      g_unit <- 1 / sds
      v_unit <- var(drop(covs %*% g_unit))
      gamma_out <- sqrt((signal / (1 - signal)) / v_unit) * g_unit
    }
  } else {
    if (!is.numeric(gamma) || length(gamma) != k)
      stop(sprintf("`gamma` must be a numeric vector of length %d.", k), call. = FALSE)
    gamma_out <- gamma
  }

  idxz <- if (k > 0L) drop(covs %*% gamma_out) else rep(0, nrow(M))
  s2 <- if (k > 0L) var(idxz) else 0
  f0 <- function(b0) mean(lf$G(b0 + idxz)) - baseline
  beta0 <- uniroot(f0, c(-15, 15), extendInt = "upX", tol = 1e-10)$root

  focal_ref <- if (focal_var$type == "binary") 0
               else if (!is.null(emp_xd)) mean(emp_xd)
               else focal_var$mean

  structure(list(
    model = model, link = model, G = lf$G, g = lf$g,
    focal = focal_var, covariates = colnames(covs), k = k,
    baseline = baseline,
    signal = if (is.null(gamma)) signal else NA_real_,
    signal_implied = s2 / (1 + s2), s2 = s2,
    gamma = gamma_out, beta0 = beta0,
    beta_focal = NULL, beta_mod = NULL, beta_int = NULL,
    moderator = NULL, mod_ref = NULL,
    estimand = NULL, target_est = NULL, p1 = NULL,
    main_focal = NULL, main_moderator = NULL,
    route = "empirical", builder = "empirical",
    emp_x = covs, emp_xd = emp_xd,
    n_emp = nrow(M),
    focal_ref = focal_ref,
    n_int = 4e5, seed_int = seed_int,
    spec = list(model = model, data = dat, focal = focal,
                baseline = baseline, signal = signal, gamma = gamma,
                seed_int = seed_int)
  ), class = "powerape_dgp")
}

# Engine draw for the empirical route: resample pilot rows; an independent
# pa_var focal is drawn from its parametric marginal on top.
draw_x_empirical <- function(dgp, n) {
  idx <- sample.int(dgp$n_emp, n, replace = TRUE)
  d <- if (!is.null(dgp$emp_xd)) dgp$emp_xd[idx] else transform_u(rnorm(n), dgp$focal)
  X <- if (dgp$k > 0L) cbind(1, d, dgp$emp_x[idx, , drop = FALSE]) else cbind(1, d)
  list(X = X, focal_col = 2L)
}
