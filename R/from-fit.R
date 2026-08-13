# Pilot-model route: lift the DGP from a fitted glm ---------------------------

#' Specify a DGP from a fitted pilot model
#'
#' The third covariate route (DESIGN.md section 3.1): everything except the
#' focal coefficient is lifted from a fitted pilot `glm` — the covariate
#' rows come from the pilot's model matrix (resampled with replacement, so
#' their joint distribution is preserved), the nuisance coefficients are
#' the pilot estimates, and the intercept is the pilot's unless `baseline`
#' overrides it. The focal coefficient itself is *not* taken from the
#' pilot: it is re-solved in APE units by [set_ape()], keeping the "power
#' for the effect size you care about" logic intact.
#'
#' The pilot fit must be a binomial `glm` with a probit or logit link and a
#' main-effects role for the focal variable (terms interacting with the
#' focal variable are not supported here — use [ape_dgp()] + [set_aie()]
#' for interaction designs). Nuisance-side interactions and transformed
#' covariates are fine: they are just columns of the model matrix.
#'
#' @param fit A fitted `glm` with `family = binomial("probit")` or
#'   `binomial("logit")`.
#' @param focal Name of the focal regressor: an untransformed numeric
#'   variable in the model (binary 0/1 or continuous).
#' @param baseline Optional override, in (0, 1): recalibrate the intercept
#'   so the average `P(Y = 1)` with the focal at reference equals this
#'   value. Default `NULL` keeps the pilot-implied baseline (reported in
#'   the object).
#'
#' @return An object of class `powerape_dgp` (route `"empirical"`, builder
#'   `"from_fit"`). Pin the effect with [set_ape()], then use the usual
#'   power functions; [ape_robust()] can vary `baseline`.
#' @examples
#' set.seed(1)
#' pilot <- data.frame(treat = rbinom(600, 1, 0.5),
#'                     age = rnorm(600, 45, 12))
#' pilot$y <- rbinom(600, 1, pnorm(-0.6 + 0.3 * pilot$treat + 0.01 * pilot$age))
#' fit <- glm(y ~ treat + age, binomial("probit"), pilot)
#' d <- ape_dgp_from_fit(fit, focal = "treat")
#' d <- set_ape(d, target = 0.10)
#' d
#' @export
ape_dgp_from_fit <- function(fit, focal, baseline = NULL) {
  if (!inherits(fit, "glm") || !identical(fit$family$family, "binomial"))
    stop("`fit` must be a glm with a binomial family.", call. = FALSE)
  link <- fit$family$link
  if (!link %in% c("probit", "logit"))
    stop("`fit` must use a probit or logit link.", call. = FALSE)
  stopifnot(is.character(focal), length(focal) == 1L)
  if (!is.null(baseline))
    stopifnot(is.numeric(baseline), length(baseline) == 1L,
              baseline > 0, baseline < 1)

  X <- stats::model.matrix(fit)
  if (!"(Intercept)" %in% colnames(X))
    stop("The pilot model must include an intercept.", call. = FALSE)
  if (!focal %in% colnames(X))
    stop(sprintf(paste("Focal variable \"%s\" is not a column of the pilot",
                       "model matrix. It must enter the model as an",
                       "untransformed numeric main effect."), focal),
         call. = FALSE)
  trm <- stats::terms(fit)
  fac <- attr(trm, "factors")
  if (focal %in% rownames(fac)) {
    inv <- colnames(fac)[fac[focal, ] > 0]
    if (any(attr(trm, "order")[match(inv, colnames(fac))] > 1L))
      stop(paste("The pilot-model route supports the focal variable as a",
                 "main effect only; interactions with the focal variable",
                 "are not supported (use ape_dgp() + set_aie() for",
                 "interaction designs)."), call. = FALSE)
  }

  co <- stats::coef(fit)
  if (anyNA(co))
    stop("The pilot fit has NA coefficients (rank deficiency); fix the pilot model first.",
         call. = FALSE)
  w <- fit$prior.weights
  if (!is.null(w) && any(w != 1))
    warning("Prior weights in the pilot fit are ignored; rows are resampled unweighted.")

  xd <- X[, focal]
  if (length(unique(xd)) < 2L)
    stop("The focal column is constant in the pilot data.", call. = FALSE)
  covs <- X[, setdiff(colnames(X), c("(Intercept)", focal)), drop = FALSE]
  n_emp <- nrow(X)
  if (n_emp < 30L)
    stop("The pilot-model route needs at least 30 pilot rows.", call. = FALSE)
  if (n_emp < 500L)
    warning("Fewer than 500 pilot rows: pilot sampling noise becomes part of the design (cf. Hancock & Feng, 2025).")

  focal_var <- if (all(xd %in% c(0, 1))) {
    pa_var(focal, "binary", p = mean(xd))
  } else {
    pa_var(focal, "normal", mean = mean(xd), sd = sd(xd))
  }
  focal_ref <- if (focal_var$type == "binary") 0 else mean(xd)

  lf <- link_funs(link)
  k <- ncol(covs)
  gamma <- if (k > 0L) unname(co[colnames(covs)]) else numeric(0)
  idxz <- if (k > 0L) drop(covs %*% gamma) else rep(0, n_emp)
  s2 <- if (k > 0L) var(idxz) else 0

  ## Centered intercept: pilot intercept plus the fitted focal contribution
  ## at the reference value; or recalibrated to a baseline override.
  if (is.null(baseline)) {
    beta0 <- unname(co[["(Intercept)"]]) + unname(co[[focal]]) * focal_ref
    baseline_out <- mean(lf$G(beta0 + idxz))
  } else {
    f0 <- function(b0) mean(lf$G(b0 + idxz)) - baseline
    beta0 <- uniroot(f0, c(-15, 15), extendInt = "upX", tol = 1e-10)$root
    baseline_out <- baseline
  }

  structure(list(
    model = link, link = link, G = lf$G, g = lf$g,
    focal = focal_var, moderator = NULL, covariates = colnames(covs), k = k,
    baseline = baseline_out,
    signal = NA_real_,
    signal_implied = s2 / (1 + s2), s2 = s2,
    gamma = gamma, beta0 = beta0,
    beta_focal = NULL, beta_mod = NULL, beta_int = NULL,
    estimand = NULL, target_est = NULL, p1 = NULL,
    main_focal = NULL, main_moderator = NULL,
    route = "empirical", builder = "from_fit",
    emp_x = covs, emp_xd = xd, n_emp = n_emp,
    focal_ref = focal_ref, mod_ref = NULL,
    formula_str = paste(deparse(stats::formula(fit)), collapse = " "),
    n_int = 4e5, seed_int = 20260812L,
    spec = list(fit = fit, focal = focal, baseline = baseline)
  ), class = "powerape_dgp")
}
