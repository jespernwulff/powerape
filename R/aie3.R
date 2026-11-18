# Three-way average interaction effects -----------------------------------------
# A parametric design with a focal variable, a
# moderator, and a second moderator w that is binary. The estimand is
# ginteff's three-way average interaction effect: the two-way AIE of
# (focal, moderator) with w set to 1 for every observation minus the same
# with w set to 0. The index, with the focal and the moderator centred at
# their reference values (xd, md), is
#   b0 + bf xd + bm md + bi xd md
#      + w (b2 + bf2 xd + bm2 md + bi3 xd md) + z'gamma,
# and the simulated studies are fitted on the design matrix
#   [1, d, m, d m, w, d w, m w, d m w, Z]   (columns 2..8, Z from 9).
# Every function that treats a design with a moderator as a two-way AIE
# design checks is_aie3() first.

is_aie3 <- function(dgp) !is.null(dgp$moderator2)

# Functions outside the wired set (ape_power, ape_curve, ape_n, ape_mde)
# stop instead of treating a three-way design as a two-way one.
aie3_refuse <- function(dgp, what) {
  if (is_aie3(dgp))
    stop(sprintf(paste("%s is not available for three-way interaction designs",
                       "yet; ape_power(), ape_curve(), ape_n() and ape_mde()",
                       "are."), what), call. = FALSE)
  invisible(TRUE)
}

# "treat x female x public", as named in prints.
aie3_term <- function(dgp) {
  sprintf("%s x %s x %s", dgp$focal$name, dgp$moderator$name, dgp$moderator2$name)
}

# Deterministic draw of all regressors (focal, moderator, second moderator,
# covariates) from the Gaussian copula; reproduced exactly from seed_int.
reg_draw_3way <- function(focal, moderator, moderator2, covariates, R, n_int,
                          seed_int) {
  k <- length(covariates)
  with_seed(seed_int, {
    U <- draw_std_mvn(n_int, chol(R))
    Z <- if (k > 0L) {
      Zm <- U[, 3L + seq_len(k), drop = FALSE]
      for (j in seq_len(k)) Zm[, j] <- transform_u(Zm[, j], covariates[[j]])
      Zm
    } else {
      matrix(numeric(0), n_int, 0L)
    }
    list(xd = transform_u(U[, 1], focal), md = transform_u(U[, 2], moderator),
         wd = transform_u(U[, 3], moderator2), Z = Z)
  })
}

# The builder behind ape_dgp(moderator2 = ). Same anchoring as ape_dgp():
# the baseline is the outcome rate with all three design variables at their
# reference values, averaged over the covariates.
ape_dgp_3way <- function(model, focal, moderator, moderator2, covariates,
                         correlation, baseline, signal, gamma, n_int, seed_int) {
  if (inherits(covariates, "pa_var")) covariates <- list(covariates)
  stopifnot(inherits(focal, "pa_var"),
            all(vapply(covariates, inherits, logical(1), "pa_var")))
  if (is.null(moderator))
    stop(paste("`moderator2` needs a `moderator`: a three-way design has a",
               "focal variable and two moderators."), call. = FALSE)
  stopifnot(inherits(moderator, "pa_var"))
  if (!inherits(moderator2, "pa_var") || !identical(moderator2$type, "binary"))
    stop(paste("`moderator2` must be a binary pa_var(). The three-way effect is",
               "symmetric in its three variables, so place a binary one last;",
               "a design with three continuous variables is not available."),
         call. = FALSE)
  nm <- c(focal$name, moderator$name, moderator2$name)
  if (anyDuplicated(nm))
    stop("`focal`, `moderator` and `moderator2` must be three different variables.",
         call. = FALSE)
  stopifnot(is.numeric(baseline), length(baseline) == 1L,
            baseline > 0, baseline < 1)
  stopifnot(is.numeric(n_int), n_int >= 1e4)

  k <- length(covariates)
  K <- 3L + k
  if (is.matrix(correlation) && !all(dim(correlation) == c(K, K)))
    stop(sprintf(paste("`correlation` must be a %d x %d matrix over (focal,",
                       "moderator, moderator2, covariates)."), K, K), call. = FALSE)
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
    Z <- reg_draw_3way(focal, moderator, moderator2, covariates, R, n_int,
                       seed_int)$Z
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
    focal = focal, moderator = moderator, moderator2 = moderator2,
    covariates = covariates, k = k,
    R = R, R_chol = chol(R),
    baseline = baseline,
    signal = if (is.null(gamma)) signal else NA_real_,
    signal_implied = signal_implied, s2 = s2,
    gamma = gamma_out, beta0 = beta0,
    beta_focal = NULL, beta_mod = NULL, beta_int = NULL,
    beta_mod2 = NULL, beta_int_f2 = NULL, beta_int_m2 = NULL, beta_int3 = NULL,
    estimand = NULL, target_est = NULL, p1 = NULL,
    main3 = NULL, two_way3 = NULL, cells3 = NULL,
    route = "parametric", builder = "parametric",
    focal_ref = if (focal$type == "binary") 0 else focal$mean,
    mod_ref = if (moderator$type == "binary") 0 else moderator$mean,
    n_int = n_int, seed_int = seed_int,
    spec = list(model = model, focal = focal, moderator = moderator,
                moderator2 = moderator2,
                covariates = covariates, correlation = correlation,
                baseline = baseline, signal = signal, gamma = gamma,
                n_int = n_int, seed_int = seed_int)
  ), class = "powerape_dgp")
}

# Integration set: centred focal and moderator draws, the second moderator's
# 0/1 draws, and the nuisance index. Deterministic. With three binary
# variables and no covariates every anchor is a function of eight index
# values, so no draw is needed (the nuisance index is the scalar 0).
aie3_integration <- function(dgp) {
  if (dgp$k == 0L && dgp$focal$type == "binary" && dgp$moderator$type == "binary")
    return(list(xd_c = NULL, md_c = NULL, wd = NULL, idxz = 0))
  rd <- reg_draw_3way(dgp$focal, dgp$moderator, dgp$moderator2, dgp$covariates,
                      dgp$R, dgp$n_int, dgp$seed_int)
  list(xd_c = rd$xd - dgp$focal_ref, md_c = rd$md - dgp$mod_ref, wd = rd$wd,
       idxz = if (dgp$k > 0L) drop(rd$Z %*% dgp$gamma) else rep(0, dgp$n_int))
}

# Two-way AIE of one pair of design variables in the centred
# parameterization: truth_aie() on a copy of the design that carries the
# pair's types and intercept. `x1`, `x2` are the pair's centred draws (not
# used for a binary variable).
aie3_pair <- function(dgp, b0, v1, v2, b1, b2, b12, x1, x2, idxz) {
  dd <- dgp
  dd$beta0 <- b0
  dd$focal <- v1
  dd$moderator <- v2
  truth_aie(dd, b1, b2, b12, list(xd_c = x1, md_c = x2, idxz = idxz))
}

# The two-way AIE of (focal, moderator) with the second moderator set to
# `w` (0 or 1) for everyone, from centred coefficients `cf` (a list with
# b0, bf, bm, bi, b2, bf2, bm2, bi3).
aie3_two_way_at <- function(dgp, cf, w, id) {
  aie3_pair(dgp, cf$b0 + w * cf$b2, dgp$focal, dgp$moderator,
            cf$bf + w * cf$bf2, cf$bm + w * cf$bm2, cf$bi + w * cf$bi3,
            id$xd_c, id$md_c, id$idxz)
}

# Product coefficient of a pair of binary variables (and of the triple
# product when the focal variable and the moderator are both binary): the
# effect is the implied rate of one cell minus known rates, hence strictly
# increasing in the coefficient, with limits f(-Inf) and f(Inf). One
# bracketed root-find instead of the grid scan solve_interaction() needs
# for effects that plateau.
solve_monotone <- function(f, target, label) {
  lo <- f(-40)
  hi <- f(40)
  eps <- 1e-4
  if (target <= lo + eps || target >= hi - eps)
    stop(sprintf(paste(
      "Target %s %.4f is infeasible for this DGP: the cell rate it implies",
      "lies outside (0, 1); the attainable range is (%.4f, %.4f). Adjust the",
      "target, the other anchors, or the baseline."),
      label, target, lo, hi), call. = FALSE)
  uniroot(function(b) f(b) - target, c(-8, 8), extendInt = "upX", tol = 1e-10)$root
}

# The triple-product coefficient from the three-way target, the seven
# lower-order coefficients in `cf` held fixed.
aie3_solve_top <- function(dgp, cf, target, id) {
  cf$bi3 <- 0
  a0 <- aie3_two_way_at(dgp, cf, 0, id)
  f3 <- function(b) { cf$bi3 <- b; aie3_two_way_at(dgp, cf, 1, id) - a0 }
  fb <- dgp$focal$type == "binary"
  mb <- dgp$moderator$type == "binary"
  if (fb && mb) return(solve_monotone(f3, target, "three-way AIE"))
  s_d <- if (fb) 1 else pa_var_sd(dgp$focal)
  s_m <- if (mb) 1 else pa_var_sd(dgp$moderator)
  solve_interaction(f3, target, s_d * s_m, "three-way AIE")
}

# Re-pin a pinned design at another three-way target. The six lower-order
# anchors and their coefficients do not depend on the target, so only the
# triple-product coefficient is solved again (the searches of ape_mde()).
aie3_repin <- function(dgp, target) {
  stopifnot(!is.null(dgp$beta_int3), is.numeric(target), length(target) == 1L,
            is.finite(target))
  id <- aie3_integration(dgp)
  dgp$beta_int3 <- aie3_solve_top(dgp, aie3_coefs(dgp), target, id)
  dgp$target_est <- target
  dgp$cells3 <- if (!is.null(dgp$cells3)) aie3_cells(dgp, id) else NULL
  dgp
}

aie3_coefs <- function(dgp) {
  list(b0 = dgp$beta0, bf = dgp$beta_focal, bm = dgp$beta_mod, bi = dgp$beta_int,
       b2 = dgp$beta_mod2, bf2 = dgp$beta_int_f2, bm2 = dgp$beta_int_m2,
       bi3 = dgp$beta_int3)
}

# All anchors recomputed from the pinned coefficients (tests, print).
aie3_anchors <- function(dgp, id = aie3_integration(dgp)) {
  cf <- aie3_coefs(dgp)
  G <- dgp$G
  g <- dgp$g
  xz <- id$idxz
  main_of <- function(v, b, x) {
    if (v$type == "binary") mean(G(cf$b0 + b + xz)) - mean(G(cf$b0 + xz))
    else b * mean(g(cf$b0 + b * x + xz))
  }
  f <- dgp$focal; m <- dgp$moderator; w <- dgp$moderator2
  list(
    main = c(main_of(f, cf$bf, id$xd_c), main_of(m, cf$bm, id$md_c),
             main_of(w, cf$b2, NULL)),
    two_way = c(
      aie3_pair(dgp, cf$b0, f, m, cf$bf, cf$bm, cf$bi, id$xd_c, id$md_c, xz),
      aie3_pair(dgp, cf$b0, f, w, cf$bf, cf$b2, cf$bf2, id$xd_c, NULL, xz),
      aie3_pair(dgp, cf$b0, m, w, cf$bm, cf$b2, cf$bm2, id$md_c, NULL, xz)),
    target = aie3_two_way_at(dgp, cf, 1, id) - aie3_two_way_at(dgp, cf, 0, id))
}

# The eight counterfactual cell rates of an all-binary design, as an array
# indexed [focal, moderator, moderator2] (first index = 0).
aie3_cells <- function(dgp, id = aie3_integration(dgp)) {
  cf <- aie3_coefs(dgp)
  out <- array(NA_real_, c(2, 2, 2))
  for (d in 0:1) for (m in 0:1) for (w in 0:1) {
    idx <- cf$b0 + cf$bf * d + cf$bm * m + cf$bi * d * m +
      w * (cf$b2 + cf$bf2 * d + cf$bm2 * m + cf$bi3 * d * m)
    out[d + 1L, m + 1L, w + 1L] <- mean(dgp$G(idx + id$idxz))
  }
  out
}

#' Pin the three-way interaction effect of a design with two moderators
#'
#' For a design built with `ape_dgp(moderator = , moderator2 = )`, solves
#' all seven effect coefficients from targets in probability points. The
#' pinned estimand is the **three-way average interaction effect** (as in
#' the `ginteff` package): the two-way average interaction effect of the
#' focal variable and the moderator when the second moderator is set to 1
#' for everyone, minus the same when it is set to 0, each averaged over
#' the covariates. With three binary variables this is the triple
#' difference of the eight counterfactual cell rates. The effect is
#' symmetric in its three variables.
#'
#' All anchors are explicit; there are no silent defaults:
#'
#' * `main`: the average partial effect of each variable with the other
#'   two held at their reference values (0 for binary, the mean for
#'   continuous), in the order focal, moderator, second moderator.
#' * `two_way`: the two-way average interaction effect of each pair with
#'   the third variable held **at its reference value** (not averaged over
#'   it), in the order focal x moderator, focal x second moderator,
#'   moderator x second moderator.
#' * `target`: the three-way effect. The focal-by-moderator effect is
#'   `two_way[1]` when the second moderator is 0 and `two_way[1] + target`
#'   when it is 1.
#'
#' A two-way anchor of 0 says that the pair does not interact in
#' probability points when the third variable is at its reference value.
#' It does not say that the product's coefficient is 0: in a probit or
#' logit, main effects alone produce interaction effects on the
#' probability scale, of every order. Errors report the attainable range
#' when a target is infeasible.
#'
#' Three-way designs run through [ape_power()], [ape_curve()], [ape_n()]
#' and [ape_mde()] (which searches the three-way effect and holds the six
#' other anchors), under the index analysis or, with [set_analysis()], a
#' linear probability analysis whose estimate is the coefficient of the
#' triple product. A three-way effect is expensive: with balanced binary
#' variables and the same effect in probability points, each interaction
#' order needs about four times the observations of the order below.
#'
#' @section Computing time:
#' A design with three binary variables and no covariates is pinned at
#' once. With covariates the anchors are integrals over the covariate
#' distribution, and at the default integration size of [ape_dgp()]
#' (`n_int = 4e5`) pinning takes about 15 seconds with three binary
#' variables, about 40 seconds when the focal variable or the moderator is
#' continuous, and about a minute when both are (effects of continuous
#' variables level off as their coefficient grows, so their coefficients
#' are found by a grid scan). [ape_mde()] solves only the triple product's
#' coefficient again at each candidate effect (about 10 seconds with a
#' continuous variable). A smaller `n_int` speeds up exploration.
#'
#' @param dgp A `powerape_dgp` from [ape_dgp()] with a `moderator` and a
#'   `moderator2`.
#' @param target Assumed true three-way average interaction effect (the
#'   planning value), in probability points.
#' @param main Numeric vector of three: the main-effect anchors (focal,
#'   moderator, second moderator).
#' @param two_way Numeric vector of three: the two-way anchors (focal x
#'   moderator, focal x second moderator, moderator x second moderator).
#'
#' @return The `dgp` with all effect coefficients and `target_est` filled
#'   in; `estimand` becomes `"aie3"`.
#' @examples
#' d <- ape_dgp(focal = pa_var("treat", "binary", p = 0.5),
#'              moderator = pa_var("female", "binary", p = 0.5),
#'              moderator2 = pa_var("public", "binary", p = 0.5),
#'              baseline = 0.30)
#' d <- set_aie3(d, target = 0.10, main = c(0.08, 0.04, 0.02),
#'               two_way = c(0.03, 0, 0))
#' d
#' true_aie3(d)
#' @export
set_aie3 <- function(dgp, target, main, two_way) {
  stopifnot(inherits(dgp, "powerape_dgp"))
  if (!is_aie3(dgp))
    stop(paste("This DGP has no second moderator: build it with",
               "ape_dgp(moderator = , moderator2 = ), or pin a two-way design",
               "with set_aie()."), call. = FALSE)
  if (missing(main) || missing(two_way))
    stop(paste("set_aie3() needs every anchor: `main` (three main effects, each",
               "with the other two variables at their reference values) and",
               "`two_way` (three two-way interaction effects, each with the",
               "third variable at its reference value)."), call. = FALSE)
  stopifnot(is.numeric(target), length(target) == 1L, is.finite(target))
  if (!is.numeric(main) || length(main) != 3L || !all(is.finite(main)))
    stop("`main` must hold three numbers: focal, moderator, second moderator.",
         call. = FALSE)
  if (!is.numeric(two_way) || length(two_way) != 3L || !all(is.finite(two_way)))
    stop(paste("`two_way` must hold three numbers: focal x moderator, focal x",
               "second moderator, moderator x second moderator."), call. = FALSE)
  main <- as.numeric(main)
  two_way <- as.numeric(two_way)

  G <- dgp$G
  id <- aie3_integration(dgp)
  xz <- id$idxz
  b0 <- dgp$beta0
  p0 <- mean(G(b0 + xz))
  f <- dgp$focal; m <- dgp$moderator; w <- dgp$moderator2

  ## 1. main effects, each with the other two variables at reference (all
  ##    product terms drop out)
  solve_main <- function(v, tgt, x) {
    lab <- sprintf("main-effect APE (%s)", v$name)
    if (v$type == "binary") solve_binary_shift(dgp, tgt, xz, p0, lab)
    else solve_cont_main(dgp, tgt, x, xz, pa_var_sd(v), lab)
  }
  bf <- solve_main(f, main[1], id$xd_c)
  bm <- solve_main(m, main[2], id$md_c)
  b2 <- solve_main(w, main[3], NULL)

  ## 2. two-way coefficients, each from its pair's AIE with the third
  ##    variable at reference (the three-way term and the two other
  ##    products drop out)
  solve_pair <- function(v1, v2, b1, b2_, x1, x2, tgt) {
    lab <- sprintf("two-way AIE (%s x %s)", v1$name, v2$name)
    fn <- function(b) aie3_pair(dgp, b0, v1, v2, b1, b2_, b, x1, x2, xz)
    if (v1$type == "binary" && v2$type == "binary")
      return(solve_monotone(fn, tgt, lab))
    s1 <- if (v1$type == "binary") 1 else pa_var_sd(v1)
    s2 <- if (v2$type == "binary") 1 else pa_var_sd(v2)
    solve_interaction(fn, tgt, s1 * s2, lab)
  }
  bi <- solve_pair(f, m, bf, bm, id$xd_c, id$md_c, two_way[1])
  bf2 <- solve_pair(f, w, bf, b2, id$xd_c, NULL, two_way[2])
  bm2 <- solve_pair(m, w, bm, b2, id$md_c, NULL, two_way[3])

  ## 3. the triple-product coefficient from the three-way target
  bi3 <- aie3_solve_top(dgp, list(b0 = b0, bf = bf, bm = bm, bi = bi, b2 = b2,
                                  bf2 = bf2, bm2 = bm2, bi3 = 0), target, id)

  dgp$beta_focal <- bf
  dgp$beta_mod <- bm
  dgp$beta_int <- bi
  dgp$beta_mod2 <- b2
  dgp$beta_int_f2 <- bf2
  dgp$beta_int_m2 <- bm2
  dgp$beta_int3 <- bi3
  dgp$main3 <- main
  dgp$two_way3 <- two_way
  dgp$estimand <- "aie3"
  dgp$target_est <- target
  dgp$p1 <- NULL
  dgp$cells3 <- if (f$type == "binary" && m$type == "binary") aie3_cells(dgp, id)
                else NULL
  dgp
}

#' Recompute the true three-way interaction effect implied by a pinned DGP
#'
#' @param dgp A `powerape_dgp` after [set_aie3()].
#' @return The true three-way average interaction effect (numeric scalar),
#'   evaluated on the deterministic integration set; equals the
#'   [set_aie3()] target up to root-finding tolerance.
#' @export
true_aie3 <- function(dgp) {
  stopifnot(inherits(dgp, "powerape_dgp"))
  if (!is_aie3(dgp) || !identical(dgp$estimand, "aie3"))
    stop("This DGP has no pinned three-way effect - call set_aie3() first.",
         call. = FALSE)
  id <- aie3_integration(dgp)
  cf <- aie3_coefs(dgp)
  aie3_two_way_at(dgp, cf, 1, id) - aie3_two_way_at(dgp, cf, 0, id)
}

# True coefficients in design-matrix (raw, uncentred) order. The base block
# (b0, bf, bm, bi) and the block that switches on with the second moderator
# (b2, bf2, bm2, bi3) absorb the centring by the same linear map.
aie3_beta <- function(dgp) {
  if (is.null(dgp$beta_int3))
    stop("This DGP has no target effect yet - call set_aie3() first.",
         call. = FALSE)
  rd <- dgp$focal_ref
  rm_ <- dgp$mod_ref
  blk <- function(c0, cf, cm, ci) {
    c(c0 - cf * rd - cm * rm_ + ci * rd * rm_, cf - ci * rm_, cm - ci * rd, ci)
  }
  c(blk(dgp$beta0, dgp$beta_focal, dgp$beta_mod, dgp$beta_int),
    blk(dgp$beta_mod2, dgp$beta_int_f2, dgp$beta_int_m2, dgp$beta_int3),
    dgp$gamma)
}

# One simulated design matrix of size n:
#   [1, d, m, d m, w, d w, m w, d m w, Z]
draw_x_aie3 <- function(dgp, n) {
  U <- draw_std_mvn(n, dgp$R_chol)
  d <- transform_u(U[, 1], dgp$focal)
  m <- transform_u(U[, 2], dgp$moderator)
  w <- transform_u(U[, 3], dgp$moderator2)
  Z <- if (dgp$k > 0L) {
    Zm <- U[, 3L + seq_len(dgp$k), drop = FALSE]
    for (j in seq_len(dgp$k)) Zm[, j] <- transform_u(Zm[, j], dgp$covariates[[j]])
    Zm
  } else NULL
  list(X = cbind(1, d, m, d * m, w, d * w, m * w, d * m * w, Z), focal_col = 2L)
}

# Three-way AIE and its gradient from fitted (raw) coefficients: the
# two-way estimator applied to the coefficients collapsed at w = 1 and at
# w = 0, on the columns [1, d, m, d m, Z].
#   T(beta) = A(beta[1:4] + beta[5:8], beta_Z) - A(beta[1:4], beta_Z)
aie3_est <- function(beta, X, link, focal_type, mod_type) {
  p <- length(beta)
  zc <- if (p > 8L) 9:p else integer(0)
  Xc <- X[, c(1:4, zc), drop = FALSE]
  bz <- beta[zc]
  e1 <- aie_est(c(beta[1:4] + beta[5:8], bz), Xc, link, focal_type, mod_type)
  e0 <- aie_est(c(beta[1:4], bz), Xc, link, focal_type, mod_type)
  jac <- numeric(p)
  jac[1:4] <- e1$jac[1:4] - e0$jac[1:4]
  jac[5:8] <- e1$jac[1:4]
  if (length(zc)) jac[zc] <- e1$jac[4L + seq_along(zc)] - e0$jac[4L + seq_along(zc)]
  list(ape = e1$ape - e0$ape, jac = jac)
}

# Separation in the cells that identify a three-way effect: every cell of
# the cross-classification of the binary design variables (the second
# moderator always, the moderator and the focal variable when binary) must
# hold events and non-events. Same conventions as cell_check().
aie3_cell_check <- function(y, X, fb, mb) {
  g <- as.integer(X[, 5L])
  nb <- 2L
  if (mb) { g <- g + nb * as.integer(X[, 3L]); nb <- 2L * nb }
  if (fb) { g <- g + nb * as.integer(X[, 2L]); nb <- 2L * nb }
  gi <- g + 1L
  nn <- tabulate(gi, nbins = nb)
  ev <- tabulate(gi[y == 1], nbins = nb)
  keep <- nn > 0L                  # absent cells make the design rank deficient
  nn <- nn[keep]
  ev <- ev[keep]
  list(sep = any(ev == 0L | ev == nn), min_cell = min(pmin(ev, nn - ev)))
}

# The design print's lines for a three-way design.
print_aie3_effects <- function(x) {
  f <- x$focal$name; m <- x$moderator$name; w <- x$moderator2$name
  cat(sprintf("  true three-way AIE (%s) pinned at %+.4f (beta = %.4f)\n",
              aie3_term(x), x$target_est, x$beta_int3))
  cat(sprintf("  main APEs at reference: %s %+.3f, %s %+.3f, %s %+.3f\n",
              f, x$main3[1], m, x$main3[2], w, x$main3[3]))
  cat(sprintf(paste0("  two-way AIEs with the third variable at reference: ",
                     "%s x %s %+.3f, %s x %s %+.3f, %s x %s %+.3f\n"),
              f, m, x$two_way3[1], f, w, x$two_way3[2], m, w, x$two_way3[3]))
  if (!is.null(x$cells3)) {
    cl <- x$cells3
    lab <- character(0)
    for (wv in 0:1) for (mv in 0:1) for (dv in 0:1)
      lab <- c(lab, sprintf("%d%d%d %.3f", dv, mv, wv, cl[dv + 1L, mv + 1L, wv + 1L]))
    cat(sprintf("  cell rates P(Y=1 | %s, %s, %s): %s\n", f, m, w,
                paste(lab, collapse = " | ")))
  }
  invisible(NULL)
}
