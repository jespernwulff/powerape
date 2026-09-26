# Simulation engine and power computation -------------------------------------

# Error-rate convention (v1.8.0; DESIGN.md section 14). `alpha` is the
# claim's error rate, each claim in its conventional form: detection is the
# two-sided test at alpha (1 - alpha interval, directional counting);
# the minimum-effect claim is a one-sided test at alpha and equivalence a
# TOST at alpha (both read off the 1 - 2 alpha interval). A user-supplied
# `conf` overrides by naming the claim's own interval level directly.
# Since 1.9.0 `alpha` may also be a sample-size rule (a function of n;
# DESIGN.md section 15), evaluated at the design's `n`.
claim_levels <- function(claim, alpha = NULL, conf = NULL, n = NULL) {
  if (is.function(alpha)) {
    if (!is.null(conf))
      stop("Supply either an `alpha` rule or `conf`, not both.", call. = FALSE)
    if (is.null(n))
      stop("A sample-size rule for `alpha` needs `n` to be evaluated.",
           call. = FALSE)
    a <- alpha(n)
    if (!is.numeric(a) || length(a) != 1L || !is.finite(a) || a <= 0 || a >= 0.5)
      stop(sprintf(paste("The `alpha` rule must return a single number in",
                         "(0, 0.5); at n = %s it returned %s."),
                   format(n), paste(format(a), collapse = ", ")),
           call. = FALSE)
    alpha <- a
  }
  if (!is.null(conf)) {
    stopifnot(is.numeric(conf), length(conf) == 1L, conf > 0.5, conf < 1)
    alpha <- if (identical(claim, "detect")) 1 - conf else (1 - conf) / 2
  } else {
    if (is.null(alpha)) alpha <- 0.05
    stopifnot(is.numeric(alpha), length(alpha) == 1L, alpha > 0, alpha < 0.5)
  }
  list(alpha = alpha,
       conf_det = 1 - alpha,
       conf_one = 1 - 2 * alpha,
       conf_claim = if (identical(claim, "detect")) 1 - alpha else 1 - 2 * alpha)
}

# Separation in the cells that identify the estimand (DESIGN.md section 16):
# a level of a binary focal variable -- or, for AIE designs, a cell of the
# binary focal x moderator layout -- whose outcomes are all 0 or all 1. The
# maximum-likelihood estimate of the effect does not exist there: Stata's
# probit/logit drop the perfect predictor and margins reports "not
# estimable", while glm.fit silently "converges" to a large finite
# coefficient whose delta-method APE standard error collapses to the other
# cell's binomial SE (the Hauck-Donner pathology). Returns the flag and the
# smallest cell's count of events or non-events (NA when no binary cell
# identifies the estimand, e.g. a continuous focal without a binary
# moderator).
cell_check <- function(y, b1 = NULL, b2 = NULL) {
  g <- if (!is.null(b1) && !is.null(b2)) {
    2L * as.integer(b1) + as.integer(b2)
  } else if (!is.null(b1)) {
    as.integer(b1)
  } else if (!is.null(b2)) {
    as.integer(b2)
  } else {
    return(list(sep = FALSE, min_cell = NA_real_))
  }
  gi <- g + 1L
  nn <- tabulate(gi, nbins = 4L)
  ev <- tabulate(gi[y == 1], nbins = 4L)
  keep <- nn > 0L                  # absent cells make the design rank deficient
  nn <- nn[keep]
  ev <- ev[keep]
  list(sep = any(ev == 0L | ev == nn), min_cell = min(pmin(ev, nn - ev)))
}

# Backup separation check on a fitted index model: a focal-side coefficient
# (focal, moderator, interaction) whose standardized standard error is
# absurd (> 50 on the latent-index scale) marks a quasi-separated fit the
# cell check cannot see, e.g. complete separation by a continuous focal.
# Legitimate fits sit two orders of magnitude below the threshold.
hauck_donner <- function(fit, X, cols) {
  V0 <- vcov_from_glmfit(fit)
  sds <- apply(X[, cols, drop = FALSE], 2, sd)
  any(sqrt(pmax(diag(V0)[cols], 0)) * sds > 50)
}

# Core loop: nsim replications at design size n; returns the point estimate
# and delta-method SE of the APE/AIE per replication (plus CI bounds at
# `conf` for callers that want them), a usable-estimate flag `ok`, a
# separation flag `sep`, and each replication's smallest identifying-cell
# count of events or non-events (`min_cell`). Replications that fail
# (non-convergence, rank deficiency, degenerate y, runaway coefficients)
# get ok = FALSE and count against every claim downstream (conservative).
# Replications with separation count as failed under separation = "fail"
# (the default since 1.11.0: the field's reference analysis refuses such
# fits) and are kept with their degenerate Wald intervals under "keep"
# (R's glm + marginaleffects, Stata's `asis`); either way they are flagged.
sim_ci <- function(dgp, n, nsim, conf, seed = NULL, se_type = "model",
                   separation = "fail") {
  z <- zcrit(conf)
  is_aie <- identical(dgp$estimand, "aie")
  is_panel <- identical(dgp$route, "panel")
  is_iv <- identical(dgp$route, "iv")
  sep_fail <- !identical(separation, "keep")
  fb <- identical(dgp$focal$type, "binary")
  mb <- !is.null(dgp$moderator) && identical(dgp$moderator$type, "binary")
  hd_cols <- if (!is.null(dgp$moderator)) 2:4 else 2L
  bt <- if (is_panel || is_iv) NULL else beta_true(dgp)
  with_seed(seed, {
    l <- u <- se <- pt <- mc <- rep(NA_real_, nsim)
    ok <- sep <- logical(nsim)
    for (r in seq_len(nsim)) {
      xx <- draw_x(dgp, n)
      if (is_iv) {
        ## y is drawn from the latent index inside draw_x_iv (the
        ## endogeneity lives in the joint (u, v) draw)
        y <- xx$y
        cc <- cell_check(y, if (fb) xx$d, if (mb) xx$m)
        mc[r] <- cc$min_cell
        if (all(y == y[1L])) next      # degenerate outcome: failed, not "separated"
        if (cc$sep) {
          sep[r] <- TRUE
          if (sep_fail) next
        }
        ft <- cf_fit(y, xx$d, m = xx$m, Z = xx$Zx, X1 = xx$X1,
                     first = dgp$first_stage, start = xx$start)
        if (!isTRUE(ft$ok)) next
        est <- if (is_aie) {
          cf_aie_est(ft$theta, ft$W, dgp$focal$type, ft$cf_col, ft$cfm_col)
        } else {
          cf_ape_est(ft$theta, ft$W, dgp$focal$type, ft$cf_col, ft$cfm_col)
        }
        V <- ft$V_main
      } else {
        pr <- if (is_panel) xx$pr else dgp$G(drop(xx$X %*% bt))
        y <- rbinom(length(pr), 1L, pr)
        cc <- cell_check(y, if (fb) xx$X[, 2L], if (mb) xx$X[, 3L])
        mc[r] <- cc$min_cell
        if (all(y == y[1L])) next      # degenerate outcome: failed, not "separated"
        if (cc$sep) {
          sep[r] <- TRUE
          if (sep_fail) next
        }
        ft <- fit_index_model(xx$X, y, dgp$link,
                              start = if (is_panel) xx$start else bt)
        if (!ft$ok) next
        if (!sep[r] && hauck_donner(ft$fit, xx$X, hd_cols)) {
          sep[r] <- TRUE
          if (sep_fail) next
        }
        est <- if (is_aie) {
          aie_est(ft$fit$coefficients, xx$X, dgp$link,
                  dgp$focal$type, dgp$moderator$type)
        } else {
          ape_est(ft$fit$coefficients, xx$X, xx$focal_col, dgp$link,
                  dgp$focal$type)
        }
        V <- if (is_panel) {
          vcov_cluster(ft$fit, xx$X, xx$id, dgp$link)
        } else if (identical(se_type, "robust")) {
          vcov_sandwich(ft$fit, xx$X, dgp$link)
        } else {
          vcov_from_glmfit(ft$fit)
        }
      }
      s_ <- sqrt(max(0, drop(t(est$jac) %*% V %*% est$jac)))
      if (!is.finite(s_) || s_ <= 0) next
      ok[r] <- TRUE
      se[r] <- s_
      pt[r] <- est$ape
      l[r] <- est$ape - z * s_
      u[r] <- est$ape + z * s_
    }
    list(est = pt, l = l, u = u, se = se, ok = ok, sep = sep, min_cell = mc)
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
  } else if (!is.null(sesoi) &&
             (!is.numeric(sesoi) || length(sesoi) != 1L || sesoi <= 0)) {
    ## a sesoi supplied alongside claim = "detect" only breaks out the
    ## outcome table, but a bad value would corrupt that table silently
    stop("`sesoi`, when supplied, must be a single positive number (in APE units).",
         call. = FALSE)
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
summarize_sim <- function(sim, target, claim, sesoi, nsim, alpha = 0.05) {
  ## two interval widths from one estimate and SE (DESIGN.md section 14):
  ## detection reads the two-sided 1 - alpha interval, the one-sided SESOI
  ## claims read the 1 - 2 alpha interval (minimum: one-sided test at alpha;
  ## equivalence: TOST at alpha). Sign is flipped for negative targets so
  ## the hypothesized direction is "positive".
  est <- if (target < 0) -sim$est else sim$est
  se <- sim$se
  l_det <- est - qnorm(1 - alpha / 2) * se
  z1 <- qnorm(1 - alpha)
  l <- est - z1 * se
  u <- est + z1 * se
  ok <- sim$ok
  det <- ok & (l_det > 0)
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
  ## separation diagnostics (absent from synthetic draws)
  sep <- sim$sep %||% logical(length(ok))
  mc <- sim$min_cell
  list(power = power,
       mcse = sqrt(power * (1 - power) / nsim),
       outcomes = outcomes,
       n_failed = sum(!ok),
       separated = mean(sep),
       n_separated = sum(sep),
       min_cell = if (is.null(mc) || all(is.na(mc))) NA_real_ else mean(mc, na.rm = TRUE))
}

# One warning for a design whose simulated studies have separation or
# sparse identifying cells (DESIGN.md section 16). `p` carries the
# `separated` share, the mean smallest-cell count `min_cell`, and the
# `separation` convention; `where` names the design (e.g. "At n = 523").
sparse_warning <- function(p, where = NULL, aie = FALSE) {
  cell <- if (aie) "focal-by-moderator cell" else "focal cell"
  msg <- character()
  if (isTRUE(p$separated > 0.01))
    msg <- c(msg, sprintf(paste0(
      "%.1f%% of simulated studies had separation (a %s with no events or ",
      "no non-events, where the maximum-likelihood effect does not exist)%s"),
      100 * p$separated, cell,
      if (identical(p$separation, "keep"))
        "; they were kept with their degenerate Wald intervals (separation = \"keep\")"
      else "; they count as failed (separation = \"fail\", the default; Stata's probit also refuses such fits)"))
  if (isTRUE(!is.na(p$min_cell) && p$min_cell < 10))
    msg <- c(msg, sprintf(paste0(
      "the smallest %s averages %.1f events (or non-events) per simulated ",
      "study; delta-method Wald inference is unreliable at this sparsity ",
      "(one-sided error rates can run up to about twice nominal, in the ",
      "direction in which the sparse cell has the lower rate)"),
      cell, p$min_cell))
  if (length(msg))
    warning(paste0(if (!is.null(where)) paste0(where, ": ") else "",
                   paste(msg, collapse = "; "),
                   ". See ?ape_power, section 'Separation and sparse cells'."),
            call. = FALSE)
  invisible(length(msg) > 0L)
}

# Shared power computation. `enforce = FALSE` skips the coherence guards --
# used by ape_robust(), where scenario drift can legitimately push the
# implied effect across a claim boundary and the point is to SHOW that.
power_once <- function(dgp, n, claim, sesoi, conf = NULL, nsim, seed,
                       enforce = TRUE, se = "model", alpha = NULL,
                       keep_draws = FALSE, separation = "fail") {
  stopifnot(is.numeric(n), length(n) == 1L, n >= 20)
  if (identical(dgp$route, "panel") && n < 30)
    warning("Fewer than 30 units (clusters): cluster-robust inference is unreliable at this size.")
  lv <- claim_levels(claim, alpha, conf, n)
  rule <- if (is.function(alpha)) alpha else NULL
  conf <- lv$conf_claim
  stopifnot(is.numeric(nsim), length(nsim) == 1L, nsim >= 20)
  if (!identical(se, "model") && dgp$route %in% c("panel", "iv")) {
    warning(sprintf(paste("`se` is fixed by the route: panel designs use",
                          "unit-clustered SEs, IV designs the stacked robust",
                          "sandwich; `se = \"%s\"` is ignored."), se))
    se <- "model"
  }
  if (enforce) check_coherence(dgp, claim, sesoi, conf)
  sim <- sim_ci(dgp, n, nsim, conf, seed, se_type = se, separation = separation)
  res <- summarize_sim(sim, dgp$target_est, claim, sesoi, nsim, lv$alpha)
  structure(
    c(res, list(claim = claim, sesoi = sesoi, conf = conf, alpha = lv$alpha,
                alpha_rule = rule,
                n = as.integer(n),
                nsim = as.integer(nsim), target = dgp$target_est,
                estimand = dgp$estimand %||% "ape",
                se = se, separation = separation, model = dgp$model, dgp = dgp,
                draws = if (isTRUE(keep_draws)) sim[c("est", "se", "ok", "sep")] else NULL)),
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
#' **Error-rate convention.** Every claim is tested at `alpha` (default
#' 0.05) in its conventional form. Detection is the two-sided test at
#' `alpha`, read off the `1 - alpha` (95%) interval (only correctly signed
#' rejections count, which costs no power). The minimum-effect claim is a
#' one-sided test at `alpha` and equivalence a two-one-sided-tests (TOST)
#' procedure at `alpha`; both read off the `1 - 2 alpha` (90%) interval,
#' the convention of Lakens (2017), Lakens et al. (2018), TOSTER, and
#' Riesthuis (2024) for equivalence. The outcome distribution uses both
#' intervals accordingly. To reproduce the more conservative choice of a
#' 95% interval for the minimum-effect test (Riesthuis, 2024) pass
#' `conf = 0.95`, which sets that claim's one-sided error rate to 2.5%.
#' Two routes to an error rate that is chosen rather than inherited are
#' described in `vignette("justified-alpha")`: a sample-size rule passed
#' as `alpha`, and [ape_alpha()] for the error-cost optimum of Maier and
#' Lakens (2022), computed from this function's stored draws.
#'
#' @section Separation and sparse cells:
#' With rare outcomes, small samples, or a small treated group, a simulated
#' study can contain a cell that identifies the effect -- a level of a
#' binary focal variable, or for AIE designs a focal-by-moderator cell --
#' with no events or no non-events. The maximum-likelihood effect does not
#' exist there (quasi-complete separation): Stata's `probit`/`logit` drop
#' the perfect predictor and `margins` reports the effect as not
#' estimable, whereas R's `glm()` stops silently at a large finite
#' coefficient whose delta-method standard error collapses to the other
#' cell's binomial SE, so the replication usually counts as a detection.
#' By default (`separation = "fail"`) such replications count as failed,
#' against every claim, like any other replication without a usable
#' estimate; `separation = "keep"` retains them with their degenerate Wald
#' intervals (the R `glm()` + marginaleffects analysis, and Stata's `asis`
#' option; the behavior before powerape 1.11.0). Either way the share is
#' reported (`separated`), and the function warns when it exceeds 1%, or
#' when the smallest identifying cell averages fewer than 10 events (or
#' non-events) per simulated study: there the Wald test's one-sided error
#' rates can run up to about twice their nominal level in the direction in
#' which the sparse cell has the lower rate, and score, Fisher, or Firth
#' analyses have different power than the Wald analysis powered here.
#'
#' @param dgp A `powerape_dgp` after [set_ape()] or [set_aie()].
#' @param n Total sample size of the simulated study.
#' @param claim `"minimum"` (default), `"detect"`, or `"equivalence"`.
#' @param sesoi Smallest effect size of interest, in APE units. Required for
#'   `"minimum"` and `"equivalence"`; optional for `"detect"` (if supplied,
#'   the outcome table is still broken out against it).
#' @param alpha The claim's error rate (default 0.05): two-sided for
#'   `"detect"`, one-sided for `"minimum"`, TOST for `"equivalence"`.
#'   May also be a **sample-size rule**: a function of `n` returning the
#'   error rate to use at that `n`, for example
#'   `function(n) alphaN::alphaN(n, BF = 3)`, the Bayes-factor
#'   calibration of Wulff and Taylor (2024) that lowers alpha as `n`
#'   grows. The searching functions evaluate the rule at every candidate
#'   `n`, so [ape_n()] designs jointly over the pair (alpha(n), n); results
#'   store the rule and the level it produced.
#' @param conf Optional override: the interval level used for the claim
#'   itself (`1 - alpha` for detection, `1 - 2 alpha` otherwise). Supplying
#'   `conf = 0.95` for a minimum-effect or equivalence claim reproduces the
#'   pre-1.8.0 behavior (one-sided error rate 2.5%).
#' @param nsim Number of simulation replications.
#' @param seed Optional seed (the caller's RNG state is preserved).
#' @param se Standard errors for the exogenous cross-sectional routes:
#'   `"model"` (default, expected-information ML) or `"robust"`
#'   (heteroskedasticity-robust HC0 sandwich, as in the sandwich package).
#'   Panel designs always use unit-clustered SEs and IV designs the
#'   stacked method-of-moments robust sandwich; `se` is ignored there.
#' @param keep_draws Keep the per-replication estimates, standard errors,
#'   and convergence flags in the result (default TRUE; a few kilobytes).
#'   They let [power_at()] and [ape_alpha()] re-evaluate the claim at any
#'   error rate without simulating again.
#' @param separation How simulated studies with separation in an
#'   identifying cell are scored: `"fail"` (default; counted as failed, as
#'   the field's reference analysis refuses them) or `"keep"` (retained with
#'   their degenerate Wald intervals). See the section 'Separation and
#'   sparse cells'.
#'
#' @return A `powerape_power` object: power, Monte Carlo standard error,
#'   outcome distribution, the failed-fit count (failures count against
#'   power, conservatively), the share of simulated studies with separation
#'   (`separated`) and the average smallest identifying-cell count of
#'   events or non-events (`min_cell`, `NA` without binary cells), the
#'   separation convention used, the error rate used (`alpha`, and
#'   `alpha_rule` when it came from a sample-size rule), the stored draws,
#'   and the embedded DGP spec for reproducibility and
#'   [power_statement()].
#' @examples
#' \donttest{
#' d <- ape_dgp(focal = pa_var("treat", "binary", p = 0.5), baseline = 0.30)
#' d <- set_ape(d, target = 0.10)
#' ape_power(d, n = 800, claim = "minimum", sesoi = 0.03, nsim = 500, seed = 1)
#' }
#' @export
ape_power <- function(dgp, n, claim = c("minimum", "detect", "equivalence"),
                      sesoi = NULL, alpha = 0.05, conf = NULL, nsim = 1000,
                      seed = NULL, se = c("model", "robust"),
                      keep_draws = TRUE, separation = c("fail", "keep")) {
  claim <- match.arg(claim)
  se <- match.arg(se)
  separation <- match.arg(separation)
  res <- power_once(dgp, n, claim, sesoi, conf, nsim, seed, enforce = TRUE,
                    se = se, alpha = alpha, keep_draws = keep_draws,
                    separation = separation)
  sparse_warning(res, aie = identical(res$estimand, "aie"))
  res
}
