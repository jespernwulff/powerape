# Power curves and required sample size ---------------------------------------

#' Power curve over a grid of sample sizes
#'
#' Runs the power simulation at each element of `n` and collects the
#' results. Works for APE and AIE designs alike.
#'
#' @inheritParams ape_power
#' @param n Integer vector (length >= 2) of total sample sizes.
#' @param nsim Replications per grid point.
#' @param seed Optional; grid point i uses `seed + i - 1`.
#'
#' @return A `powerape_curve` object with a `results` data frame
#'   (`n`, `alpha`, `power`, `mcse`, `failed`); `alpha` varies along the
#'   grid when it is a sample-size rule.
#' @examples
#' \donttest{
#' d <- ape_dgp(focal = pa_var("treat", "binary", p = 0.5), baseline = 0.30)
#' d <- set_ape(d, target = 0.10)
#' cv <- ape_curve(d, n = c(300, 600, 900), claim = "detect", nsim = 300, seed = 1)
#' plot(cv, target_power = 0.8)
#' }
#' @export
ape_curve <- function(dgp, n, claim = c("minimum", "detect", "equivalence"),
                      sesoi = NULL, alpha = 0.05, conf = NULL, nsim = 1000,
                      seed = NULL, se = c("model", "robust")) {
  claim <- match.arg(claim)
  se <- match.arg(se)
  n <- sort(unique(as.integer(n)))
  stopifnot(length(n) >= 2L, all(n >= 20))
  lv1 <- claim_levels(claim, alpha, conf, n[1L])
  check_coherence(dgp, claim, sesoi, lv1$conf_claim)
  rows <- vector("list", length(n))
  for (i in seq_along(n)) {
    si <- if (is.null(seed)) NULL else seed + i - 1L
    p <- power_once(dgp, n[i], claim, sesoi, conf, nsim, si, enforce = FALSE,
                    se = se, alpha = alpha)
    rows[[i]] <- data.frame(n = n[i], alpha = p$alpha, power = p$power,
                            mcse = p$mcse, failed = p$outcomes[["failed"]])
  }
  rule <- is.function(alpha)
  structure(list(results = do.call(rbind, rows), claim = claim, sesoi = sesoi,
                 conf = if (rule) NA_real_ else lv1$conf_claim,
                 alpha = if (rule) NA_real_ else lv1$alpha,
                 alpha_rule = rule, nsim = nsim, target = dgp$target_est,
                 estimand = dgp$estimand %||% "ape", model = dgp$model),
            class = "powerape_curve")
}

#' @param x A `powerape_curve` object.
#' @param target_power Optional horizontal reference line; when the curve
#'   crosses it, the crossing n is marked and annotated.
#' @param ... Passed to the base plot call.
#' @rdname ape_curve
#' @export
plot.powerape_curve <- function(x, target_power = NULL, ...) {
  d <- x$results
  plot(NA, xlim = range(d$n), ylim = c(0, 1),
       xlab = "Sample size n",
       ylab = sprintf("Power (%s claim)", x$claim), ...)
  abline(h = seq(0, 1, 0.2), col = "grey92", lwd = 0.8)
  ## under an n-dependent alpha rule power need not be monotone in n, so
  ## isotonic smoothing only applies to fixed-alpha curves
  yf <- if (isTRUE(x$alpha_rule)) d$power else isoreg(d$n, d$power)$yf
  if (!is.null(target_power)) {
    abline(h = target_power, lty = 3)
    i <- which(yf >= target_power)[1L]
    if (!is.na(i) && i > 1L) {
      n_cross <- if (yf[i] > yf[i - 1L]) {
        d$n[i - 1L] + (target_power - yf[i - 1L]) *
          (d$n[i] - d$n[i - 1L]) / (yf[i] - yf[i - 1L])
      } else {
        d$n[i]
      }
      abline(v = n_cross, lty = 3)
      mtext(sprintf("n = %.0f", n_cross), side = 3, at = n_cross,
            cex = 0.8, line = 0.1)
    }
  }
  half <- 1.96 * d$mcse
  suppressWarnings(
    arrows(d$n, pmax(0, d$power - half), d$n, pmin(1, d$power + half),
           angle = 90, code = 3, length = 0.03, col = "grey40")
  )
  lines(d$n, yf, lwd = 1.5)
  points(d$n, d$power, pch = 19)
  invisible(x)
}

#' Required sample size for a target power
#'
#' Search plus confirmation. A pilot run estimates the SE scale
#' (SE ~ c/sqrt(n)), a normal-approximation formula proposes n, and
#' simulation at the proposal verifies and adjusts until the estimated power
#' is within Monte Carlo precision of the goal. By default a
#' **confirmation stage** then re-measures the candidate at `nsim_confirm`
#' replications and accepts only if the confirmed power is no more than
#' `max(0.005, 1.5 * MCSE)` below the goal, pushing n upward otherwise (up
#' to three rounds). The confirmation never trims n on overshoot -- slight
#' overshoot is the conservative direction for sample-size planning. The
#' returned power and MCSE come from the confirmation run. Works for APE
#' and AIE designs alike.
#'
#' @inheritParams ape_power
#' @param power Target power (default 0.90).
#' @param nsim Replications per search step.
#' @param n_range Admissible range for the answer; hitting the ceiling warns.
#' @param max_iter Maximum search iterations.
#' @param confirm Run the high-precision confirmation stage (default TRUE).
#'   Set to FALSE for speed in exploratory loops; the returned power then
#'   carries search-stage precision only.
#' @param nsim_confirm Replications per confirmation round (default
#'   `4 * nsim`).
#'
#' @return A `powerape_n` object: the required `n`, the confirmed (or, with
#'   `confirm = FALSE`, search-stage) power and MCSE at that `n`, the
#'   search history (`stage` column distinguishes search and confirmation
#'   rounds), and the embedded DGP.
#' @examples
#' \donttest{
#' d <- ape_dgp(focal = pa_var("treat", "binary", p = 0.5), baseline = 0.30)
#' d <- set_ape(d, target = 0.10)
#' ape_n(d, power = 0.8, claim = "detect", nsim = 600, seed = 1)
#' }
#' @export
ape_n <- function(dgp, power = 0.90, claim = c("minimum", "detect", "equivalence"),
                  sesoi = NULL, alpha = 0.05, conf = NULL, nsim = 1500,
                  seed = NULL, n_range = c(30, 2e6), max_iter = 5,
                  confirm = TRUE, nsim_confirm = 4 * nsim,
                  se = c("model", "robust")) {
  claim <- match.arg(claim)
  se <- match.arg(se)
  rule <- is.function(alpha)
  n_pilot <- 2000L
  lv_p <- claim_levels(claim, alpha, conf, n_pilot)
  stopifnot(is.numeric(power), length(power) == 1L, power > 0.5, power < 0.999)
  check_coherence(dgp, claim, sesoi, lv_p$conf_claim)
  t_abs <- abs(dgp$target_est)
  dist <- switch(claim,
                 detect = t_abs,
                 minimum = t_abs - sesoi,
                 equivalence = sesoi - t_abs)
  if (dist <= 0)
    stop("No solvable sample size: the assumed truth sits on the claim boundary.",
         call. = FALSE)
  ## the critical value belongs to the n it is evaluated at when `alpha` is
  ## a sample-size rule; otherwise it is the same at every n
  z_at <- function(n_eval) zcrit(claim_levels(claim, alpha, conf, n_eval)$conf_claim)
  zg <- qnorm(power)
  newton_n <- function(n_from, p_hat) {
    p_cl <- min(max(p_hat, 0.02), 0.998)
    denom <- max(z_at(n_from) + qnorm(p_cl), 0.1)
    n_new <- ceiling(n_from * ((z_at(n_from) + zg) / denom)^2)
    if (rule) {
      ## fixed-point pass: the proposal's own alpha(n) sets its critical value
      for (k in seq_len(8L)) {
        n_fp <- ceiling(n_from * ((z_at(max(n_new, 20L)) + zg) / denom)^2)
        if (n_fp == n_new) break
        n_new <- n_fp
      }
    }
    n_new
  }

  ps <- sim_ci(dgp, n_pilot, 300L, lv_p$conf_claim,
               seed = if (is.null(seed)) NULL else seed - 1L,
               se_type = if (dgp$route %in% c("panel", "iv")) "model" else se)
  if (mean(ps$ok) < 0.5)
    stop("Most pilot replications failed to fit; the design looks degenerate (e.g. very rare outcome).",
         call. = FALSE)
  c_hat <- mean(ps$se[ps$ok]) * sqrt(n_pilot)

  n_i <- ceiling((c_hat * (z_at(n_pilot) + zg) / dist)^2)
  if (rule) {
    for (k in seq_len(10L)) {
      n_fp <- ceiling((c_hat * (z_at(max(n_i, 20L)) + zg) / dist)^2)
      if (n_fp == n_i) break
      n_i <- n_fp
    }
  }
  hist <- data.frame(stage = character(), iter = integer(), n = integer(),
                     alpha = numeric(), power = numeric(), mcse = numeric())
  best <- NULL
  for (it in seq_len(max_iter)) {
    n_i <- as.integer(max(n_range[1], min(n_range[2], n_i)))
    si <- if (is.null(seed)) NULL else seed + it
    pw <- power_once(dgp, n_i, claim, sesoi, conf, nsim, si, enforce = FALSE,
                     se = se, alpha = alpha)
    hist <- rbind(hist, data.frame(stage = "search", iter = it, n = n_i,
                                   alpha = pw$alpha,
                                   power = pw$power, mcse = pw$mcse))
    if (is.null(best) || abs(pw$power - power) < abs(best$power - power))
      best <- list(n = n_i, power = pw$power, mcse = pw$mcse)
    if (abs(pw$power - power) <= max(0.01, 2 * pw$mcse)) break
    n_i <- newton_n(n_i, pw$power)
    if (n_i %in% hist$n) break
  }

  confirmed <- FALSE
  if (confirm) {
    stopifnot(is.numeric(nsim_confirm), length(nsim_confirm) == 1L,
              nsim_confirm >= nsim)
    n_c <- best$n
    for (cr in seq_len(3L)) {
      sc <- if (is.null(seed)) NULL else seed + 100L + cr
      pwc <- power_once(dgp, n_c, claim, sesoi, conf, nsim_confirm, sc,
                        enforce = FALSE, se = se, alpha = alpha)
      hist <- rbind(hist, data.frame(stage = "confirm", iter = cr, n = n_c,
                                     alpha = pwc$alpha,
                                     power = pwc$power, mcse = pwc$mcse))
      best <- list(n = n_c, power = pwc$power, mcse = pwc$mcse)
      if (pwc$power >= power - max(0.005, 1.5 * pwc$mcse)) {
        confirmed <- TRUE
        break
      }
      if (cr == 3L) break
      n_c <- as.integer(min(n_range[2],
                            max(newton_n(n_c, pwc$power),
                                ceiling(n_c * 1.02), n_c + 1L)))
    }
    if (!confirmed)
      warning(sprintf(paste(
        "The confirmation stage did not reach the target within tolerance;",
        "the returned n = %d has high-precision power %.3f. Increase",
        "`nsim_confirm` or treat the answer as a lower bound."),
        best$n, best$power))
  }
  if (best$n >= n_range[2])
    warning("Required n reached the n_range ceiling; treat the result as a lower bound.")
  lv_b <- claim_levels(claim, alpha, conf, best$n)
  structure(list(n = best$n, power = best$power, mcse = best$mcse,
                 goal = power, claim = claim, sesoi = sesoi,
                 conf = lv_b$conf_claim, alpha = lv_b$alpha,
                 alpha_rule = if (rule) alpha else NULL,
                 nsim = as.integer(nsim),
                 nsim_confirm = if (confirm) as.integer(nsim_confirm) else NA_integer_,
                 confirm_requested = isTRUE(confirm),
                 confirmed = confirmed,
                 history = hist,
                 target = dgp$target_est, estimand = dgp$estimand %||% "ape",
                 model = dgp$model, dgp = dgp),
            class = "powerape_n")
}
