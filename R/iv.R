# IV designs: endogenous focal via control-function probit ---------------------
# DESIGN.md section 12. Structural model: y = 1[b0 + a*d + (bm*m + bdm*d*m)
# + z'gamma + u > 0], with the focal d endogenous: a linear first stage
# (continuous d) or probit first stage (binary d) on the instruments, and
# (u, v) bivariate normal with correlation `endogeneity`. u is standard
# normal marginally, so the true ASF -- and therefore inversion and truth
# evaluation -- reuses the standard cross-sectional machinery over the
# joint (d, m, x) integration draw. Estimation in the engine replicates
# Stata's cfprobit exactly (see R/cf.R): two-step control function with
# stacked no-bootstrap robust standard errors.

# One draw of the exogenous block and the first stage. Unit-major returns.
iv_draw <- function(dgp, n, seed = NULL) {
  with_seed(seed, {
    U <- draw_std_mvn(n, dgp$R_chol)
    has_m <- !is.null(dgp$moderator)
    m <- if (has_m) transform_u(U[, 1L], dgp$moderator) else NULL
    k <- dgp$k
    X <- if (k > 0L) {
      Xm <- U[, has_m + seq_len(k), drop = FALSE]
      for (j in seq_len(k)) Xm[, j] <- transform_u(Xm[, j], dgp$covariates[[j]])
      Xm
    } else {
      matrix(numeric(0), n, 0L)
    }
    nz <- length(dgp$instruments)
    Zi <- U[, has_m + k + seq_len(nz), drop = FALSE]
    for (j in seq_len(nz)) Zi[, j] <- transform_u(Zi[, j], dgp$instruments[[j]])
    qz <- drop(Zi %*% dgp$pi_z)
    if (dgp$focal$type == "binary") {
      v <- rnorm(n)
      d <- as.integer(dgp$kappa0 + qz + v > 0)
    } else {
      v <- rnorm(n, 0, dgp$sigma_v)
      d <- dgp$kappa0 + qz + v
    }
    list(d = d, m = m, X = X, Zi = Zi, v = v)
  })
}

# Deterministic integration pieces in the standard (xd_c, md_c, idxz) shape,
# so set_ape()/set_aie()/true_ape()/true_aie() work unchanged.
iv_integration <- function(dgp) {
  rd <- iv_draw(dgp, dgp$n_int, dgp$seed_int)
  idxz <- if (dgp$k > 0L) drop(rd$X %*% dgp$gamma) else rep(0, dgp$n_int)
  list(xd_c = rd$d - dgp$focal_ref,
       md_c = if (is.null(rd$m)) NULL else rd$m - dgp$mod_ref,
       idxz = idxz)
}

#' Specify an IV DGP: endogenous focal, control-function probit
#'
#' Power analysis when the focal variable is **endogenous** and the study
#' will use instrumental variables via the control-function approach --
#' the estimator of Stata's `cfprobit` (replicated exactly, including its
#' stacked method-of-moments standard errors; no bootstrap). A continuous
#' focal gets a linear first stage, a binary focal a probit first stage;
#' the control function (first-stage residual or generalized residual)
#' joins the second-stage probit, whose coefficients margins-style ASF
#' averaging turns back into the APE/AIE.
#'
#' The estimand is unchanged: the ASF-based APE of the focal (or the AIE
#' with a binary exogenous moderator). Because the structural error is
#' standard normal marginally, targets, feasibility, and truth evaluation
#' work exactly as in [ape_dgp()].
#'
#' Two documented properties of the control-function estimator (validation
#' battery V13-V15): (1) for a **binary** endogenous focal the generalized
#' residual is the standard Wooldridge approximation, with a small
#' negative bias (about -0.01 at `endogeneity = 0.4`,
#' `iv_strength = 0.2`) that makes the directional detect claim
#' *conservative* at the null (empirical size ~1.6% against the nominal
#' 2.5%; the standard errors themselves are correctly calibrated);
#' (2) for a **continuous** endogenous focal, the margins convention
#' (control functions held at observed values) evaluates the ASF over
#' the observed (focal, residual) pairs, which understates the
#' product-measure ASF slightly (about -0.005 on a 0.08 effect at
#' `endogeneity = 0.5`, `iv_strength = 0.3`) because the focal and its
#' residual are correlated by construction. powerape replicates the
#' field-standard (cfprobit/margins) convention exactly; a
#' product-measure ASF option is on the roadmap.
#'
#' For AIE designs the estimating model includes the focal-by-moderator
#' term, the control function, and the control-function-by-moderator
#' interaction (Stata: `interact()` + `mainonly()`); the first stage
#' includes the moderator **and the instrument-by-moderator interactions**
#' -- the general IV requirement for interaction models. Inference is
#' robust (stacked sandwich) throughout.
#'
#' @inheritParams ape_dgp
#' @param focal [pa_var()], `"normal"` (linear first stage) or `"binary"`
#'   (probit first stage). The endogenous variable.
#' @param moderator Optional **binary** [pa_var()] for AIE designs.
#' @param instruments A [pa_var()] or list of them: the excluded
#'   instruments. Equal standardized first-stage weights.
#' @param endogeneity Correlation of the structural and first-stage errors
#'   (u, v) -- the reason plain probit would be biased. In (-1, 1).
#' @param iv_strength First-stage strength: the share of the focal's
#'   variance (latent variance for a binary focal) explained by the
#'   instruments. In (0, 1); small values mean weak instruments.
#' @param correlation Dependence among the exogenous block (moderator,
#'   covariates, instruments): scalar exchangeable correlation or a matrix.
#'   The focal is generated by the first stage, not drawn.
#'
#' @return An object of class `powerape_dgp` (route `"iv"`). Pin the
#'   effect with [set_ape()] or, with a moderator, [set_aie()].
#' @examples
#' \donttest{
#' d <- ape_dgp_iv(
#'   focal = pa_var("training", "binary", p = 0.4),
#'   covariates = list(pa_var("size", "normal")),
#'   instruments = pa_var("subsidy", "binary", p = 0.5),
#'   endogeneity = 0.4, iv_strength = 0.2,
#'   baseline = 0.30, signal = 0.10
#' )
#' d <- set_ape(d, target = 0.10)
#' ape_power(d, n = 1500, claim = "detect", nsim = 300, seed = 1)
#' }
#' @export
ape_dgp_iv <- function(focal, moderator = NULL, covariates = list(),
                       instruments, endogeneity, iv_strength,
                       correlation = NULL, baseline, signal = 0,
                       n_int = 1e5, seed_int = 20260814L) {
  if (inherits(covariates, "pa_var")) covariates <- list(covariates)
  if (inherits(instruments, "pa_var")) instruments <- list(instruments)
  stopifnot(inherits(focal, "pa_var"),
            all(vapply(covariates, inherits, logical(1), "pa_var")),
            length(instruments) >= 1L,
            all(vapply(instruments, inherits, logical(1), "pa_var")))
  has_m <- !is.null(moderator)
  if (has_m) {
    stopifnot(inherits(moderator, "pa_var"))
    if (moderator$type != "binary")
      stop("IV AIE designs currently require a binary moderator.", call. = FALSE)
    if (identical(moderator$name, focal$name))
      stop("`moderator` must be a different variable than `focal`.", call. = FALSE)
  }
  stopifnot(is.numeric(endogeneity), length(endogeneity) == 1L,
            abs(endogeneity) < 1)
  stopifnot(is.numeric(iv_strength), length(iv_strength) == 1L,
            iv_strength > 0, iv_strength < 1)
  stopifnot(is.numeric(baseline), length(baseline) == 1L,
            baseline > 0, baseline < 1)
  stopifnot(is.numeric(signal), length(signal) == 1L, signal >= 0, signal < 1)
  k <- length(covariates)
  if (k == 0L && signal > 0)
    stop("`signal` > 0 requires at least one covariate.", call. = FALSE)
  nz <- length(instruments)
  R <- build_corr(has_m + k + nz, correlation)
  R_chol <- chol(R)
  lf <- link_funs("probit")

  ## deterministic exogenous draw for calibration
  seed_use <- seed_int
  U <- with_seed(seed_use, draw_std_mvn(n_int, R_chol))
  Xc <- if (k > 0L) {
    Xm <- U[, has_m + seq_len(k), drop = FALSE]
    for (j in seq_len(k)) Xm[, j] <- transform_u(Xm[, j], covariates[[j]])
    Xm
  } else {
    matrix(numeric(0), n_int, 0L)
  }
  Zi <- U[, has_m + k + seq_len(nz), drop = FALSE]
  for (j in seq_len(nz)) Zi[, j] <- transform_u(Zi[, j], instruments[[j]])

  ## first-stage loadings from iv_strength (equal standardized weights)
  w_unit <- 1 / vapply(instruments, pa_var_sd, numeric(1))
  v_unit <- var(drop(Zi %*% w_unit))
  if (focal$type == "binary") {
    var_qz <- iv_strength / (1 - iv_strength)      # latent share
    pi_z <- sqrt(var_qz / v_unit) * w_unit
    sigma_v <- 1
    qz <- drop(Zi %*% pi_z)
    kappa0 <- uniroot(function(b) mean(pnorm(b + qz)) - focal$p,
                      c(-10, 10), extendInt = "upX", tol = 1e-10)$root
  } else {
    sd_d <- pa_var_sd(focal)
    pi_z <- sqrt(iv_strength * sd_d^2 / v_unit) * w_unit
    sigma_v <- sqrt((1 - iv_strength) * sd_d^2)
    qz <- drop(Zi %*% pi_z)
    kappa0 <- focal$mean - mean(qz)
  }

  ## nuisance coefficients from signal (covariates only)
  if (k > 0L) {
    if (signal == 0) {
      gamma <- rep(0, k)
    } else {
      g_unit <- 1 / vapply(covariates, pa_var_sd, numeric(1))
      v_g <- var(drop(Xc %*% g_unit))
      gamma <- sqrt((signal / (1 - signal)) / v_g) * g_unit
    }
    idxz <- drop(Xc %*% gamma)
  } else {
    gamma <- numeric(0)
    idxz <- rep(0, n_int)
  }
  s2 <- if (k > 0L) var(idxz) else 0

  ## intercept from the counterfactual baseline (focal and moderator at
  ## reference; u is standard normal marginally)
  beta0 <- uniroot(function(b) mean(lf$G(b + idxz)) - baseline,
                   c(-15, 15), extendInt = "upX", tol = 1e-10)$root

  structure(list(
    model = "probit", link = "probit", G = lf$G, g = lf$g,
    focal = focal, moderator = moderator, covariates = covariates, k = k,
    instruments = instruments,
    R = R, R_chol = R_chol,
    baseline = baseline, signal = signal,
    signal_implied = s2 / (1 + s2), s2 = s2,
    gamma = gamma, beta0 = beta0,
    beta_focal = NULL, beta_mod = NULL, beta_int = NULL,
    estimand = NULL, target_est = NULL, p1 = NULL,
    main_focal = NULL, main_moderator = NULL,
    route = "iv", builder = "iv",
    endogeneity = endogeneity, iv_strength = iv_strength,
    first_stage = if (focal$type == "binary") "probit" else "linear",
    pi_z = pi_z, kappa0 = kappa0, sigma_v = sigma_v,
    focal_ref = if (focal$type == "binary") 0 else focal$mean,
    mod_ref = if (has_m) 0 else NULL,
    n_int = as.integer(n_int), seed_int = seed_int,
    spec = list(focal = focal, moderator = moderator, covariates = covariates,
                instruments = instruments, endogeneity = endogeneity,
                iv_strength = iv_strength, correlation = correlation,
                baseline = baseline, signal = signal,
                n_int = n_int, seed_int = seed_int)
  ), class = "powerape_dgp")
}

# One simulated study of size n: exogenous block + first stage + structural
# outcome probabilities, plus everything the CF estimation needs.
draw_x_iv <- function(dgp, n) {
  if (is.null(dgp$beta_focal))
    stop("This DGP has no target effect yet - call set_ape() (or set_aie()) first.",
         call. = FALSE)
  rd <- iv_draw(dgp, n)
  rho <- dgp$endogeneity
  ## u correlated with the (standardized) first-stage error: the point of
  ## the exercise. y must be drawn from the latent index, not from
  ## Bernoulli(Phi(eta)) -- that would integrate u out and lose the
  ## endogeneity.
  v_std <- if (dgp$focal$type == "binary") rd$v else rd$v / dgp$sigma_v
  u <- rho * v_std + sqrt(1 - rho^2) * rnorm(n)
  eta <- dgp$beta0 + dgp$beta_focal * rd$d +
    (if (dgp$k > 0L) drop(rd$X %*% dgp$gamma) else 0)
  if (!is.null(rd$m)) {
    eta <- eta + dgp$beta_mod * rd$m + dgp$beta_int * rd$d * rd$m
  }
  y <- as.integer(eta + u > 0)

  ## second-stage pseudo-true start values: structural coefficients scaled
  ## by 1/sqrt(1 - rho^2); CF coefficient rho/(sigma_v sqrt(1 - rho^2))
  s2s <- sqrt(1 - rho^2)
  lam <- if (dgp$focal$type == "binary") rho / s2s else rho / (dgp$sigma_v * s2s)
  start <- if (is.null(rd$m)) {
    c(c(dgp$beta0, dgp$beta_focal, dgp$gamma) / s2s, lam)
  } else {
    c(c(dgp$beta0, dgp$beta_focal, dgp$beta_mod, dgp$beta_int, dgp$gamma) / s2s,
      lam, 0)
  }

  ## first stage: exogenous covariates, the moderator, the instruments --
  ## and, with a moderator, the instrument-by-moderator interactions
  ## (the general IV requirement for interaction models)
  X1 <- cbind(1, rd$X, rd$m, rd$Zi, if (!is.null(rd$m)) rd$Zi * rd$m)
  list(y = y, d = rd$d, m = rd$m, Zx = if (dgp$k > 0L) rd$X else NULL,
       X1 = X1, start = start)
}
