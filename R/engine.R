# Simulation engine and power computation -------------------------------------

# Core loop: nsim replications at design size n; returns CI bounds for the
# APE/AIE, the delta-method SE, and a convergence flag per replication.
# Replications that fail (non-convergence, separation, degenerate y) get
# ok = FALSE and count against every claim downstream (conservative).
sim_ci <- function(dgp, n, nsim, conf, seed = NULL) {
  z <- zcrit(conf)
  bt <- beta_true(dgp)
  is_aie <- identical(dgp$estimand, "aie")
  with_seed(seed, {
    l <- u <- se <- rep(NA_real_, nsim)
    ok <- logical(nsim)
    for (r in seq_len(nsim)) {
      xx <- draw_x(dgp, n)
      pr <- dgp$G(drop(xx$X %*% bt))
      y <- rbinom(n, 1L, pr)
      if (all(y == y[1L])) next
      ft <- fit_index_model(xx$X, y, dgp$link, start = bt)
      if (!ft$ok) next
      est <- if (is_aie) {
        aie_est(ft$fit$coefficients, xx$X, dgp$link,
                dgp$focal$type, dgp$moderator$type)
      } else {
        ape_est(ft$fit$coefficients, xx$X, xx$focal_col, dgp$link,
                dgp$focal$type)
      }
      V <- vcov_from_glmfit(ft$fit)
      s_ <- sqrt(max(0, drop(t(est$jac) %*% V %*% est$jac)))
      if (!is.finite(s_) || s_ <= 0) next
      ok[r] <- TRUE
      se[r] <- s_
      l[r] <- est$ape - z * s_
      u[r] <- est$ape + z * s_
    }
    list(l = l, u = u, se = se, ok = ok)
  })
}

# Guard rails for claim/truth coherence (DESIGN.md section 2.3).
check_coherence <- function(dgp, claim, sesoi, conf) {
  if (is.null(dgp$target_est))
    stop("This DGP has no target effect yet - call set_ape() (or set_aie()) first.",
         call. = FALSE)
  est_lab <- toupper(dgp$estimand %||% "ape")
  t_abs <- abs(dgp$target_est)
  if (claim %in% c("minimum", "equivalence")) {
    if (is.null(sesoi) || !is.numeric(sesoi) || length(sesoi) != 1L || sesoi <= 0)
      stop(sprintf("claim = \"%s\" needs a single positive `sesoi` (in APE units).",
                   claim), call. = FALSE)
  }
  if (claim == "minimum" && t_abs <= sesoi)
    stop(sprintf(paste(
      "Assumed true %s (%.4f in magnitude) does not exceed the SESOI (%.4f).",
      "On or below the null boundary of the minimum-effect test the claim",
      "succeeds with probability ~%.1f%% (the test size) at ANY sample size,",
      "so no n delivers meaningful power. Either assume a planning value above",
      "the SESOI, or run the classic SESOI analysis instead:",
      "claim = \"detect\" with the target set to the SESOI."),
      est_lab, t_abs, sesoi, 100 * (1 - conf) / 2), call. = FALSE)
  if (claim == "equivalence" && t_abs >= sesoi)
    stop(sprintf(paste(
      "Equivalence power requires |true %s| (%.4f) strictly below the SESOI",
      "(%.4f); typically the planning value is target = 0. On or beyond the",
      "boundary the equivalence claim cannot exceed the test size."),
      est_lab, t_abs, sesoi), call. = FALSE)
  invisible(TRUE)
}

# Turn simulated CIs into claim power and the outcome distribution.
# Negative targets are mirrored so claims read in the effect's direction.
summarize_sim <- function(sim, target, claim, sesoi, nsim) {
  l <- sim$l
  u <- sim$u
  if (target < 0) {
    tmp <- l
    l <- -u
    u <- -tmp
  }
  ok <- sim$ok
  det <- ok & (l > 0)
  if (!is.null(sesoi)) {
    mn <- ok & (l > sesoi)
    eq <- ok & (l > -sesoi) & (u < sesoi)
    d_only <- det & !mn & !eq
    inc <- ok & !mn & !eq & !d_only
    outcomes <- c(minimum = mean(mn), detect_only = mean(d_only),
                  inconclusive = mean(inc), equivalence = mean(eq),
                  failed = mean(!ok))
    power <- switch(claim,
                    minimum = mean(mn),
                    detect = mean(det),
                    equivalence = mean(eq))
  } else {
    outcomes <- c(detect = mean(det), inconclusive = mean(ok & !det),
                  failed = mean(!ok))
    power <- mean(det)
  }
  list(power = power,
       mcse = sqrt(power * (1 - power) / nsim),
       outcomes = outcomes,
       n_failed = sum(!ok))
}

# Shared power computation. `enforce = FALSE` skips the coherence guards --
# used by ape_robust(), where scenario drift can legitimately push the
# implied effect across a claim boundary and the point is to SHOW that.
power_once <- function(dgp, n, claim, sesoi, conf, nsim, seed, enforce = TRUE) {
  stopifnot(is.numeric(n), length(n) == 1L, n >= 20)
  stopifnot(is.numeric(conf), length(conf) == 1L, conf > 0.5, conf < 1)
  stopifnot(is.numeric(nsim), length(nsim) == 1L, nsim >= 20)
  if (enforce) check_coherence(dgp, claim, sesoi, conf)
  sim <- sim_ci(dgp, n, nsim, conf, seed)
  res <- summarize_sim(sim, dgp$target_est, claim, sesoi, nsim)
  structure(
    c(res, list(claim = claim, sesoi = sesoi, conf = conf, n = as.integer(n),
                nsim = as.integer(nsim), target = dgp$target_est,
                estimand = dgp$estimand %||% "ape",
                model = dgp$model, dgp = dgp)),
    class = "powerape_power"
  )
}

#' Simulated power for an APE/AIE claim at a given sample size
#'
#' Simulates studies from the pinned DGP, fits the index model by maximum
#' likelihood, computes the APE (or, for moderated DGPs, the AIE) with a
#' delta-method Wald CI, and evaluates the requested confidence-interval
#' claim (Riesthuis, 2024): `"minimum"` (CI lower bound above the SESOI --
#' the default and the package's point), `"detect"` (CI excludes 0,
#' directional), or `"equivalence"` (CI within +/- SESOI). Also reports the
#' full outcome distribution.
#'
#' @param dgp A `powerape_dgp` after [set_ape()] or [set_aie()].
#' @param n Total sample size of the simulated study.
#' @param claim `"minimum"` (default), `"detect"`, or `"equivalence"`.
#' @param sesoi Smallest effect size of interest, in APE units. Required for
#'   `"minimum"` and `"equivalence"`; optional for `"detect"` (if supplied,
#'   the outcome table is still broken out against it).
#' @param conf CI level (default 0.95).
#' @param nsim Number of simulation replications.
#' @param seed Optional seed (the caller's RNG state is preserved).
#'
#' @return A `powerape_power` object: power, Monte Carlo standard error,
#'   outcome distribution, the failed-fit count (failures count against
#'   power, conservatively), and the embedded DGP spec for reproducibility
#'   and [power_statement()].
#' @examples
#' \donttest{
#' d <- ape_dgp(focal = pa_var("treat", "binary", p = 0.5), baseline = 0.30)
#' d <- set_ape(d, target = 0.10)
#' ape_power(d, n = 800, claim = "minimum", sesoi = 0.03, nsim = 500, seed = 1)
#' }
#' @export
ape_power <- function(dgp, n, claim = c("minimum", "detect", "equivalence"),
                      sesoi = NULL, conf = 0.95, nsim = 1000, seed = NULL) {
  claim <- match.arg(claim)
  power_once(dgp, n, claim, sesoi, conf, nsim, seed, enforce = TRUE)
}
