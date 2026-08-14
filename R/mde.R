# Minimum detectable effect at a given sample size ----------------------------
# The inverse of ape_n(): the design (n, world, claim) is fixed and the
# question is the smallest effect it can reliably conclude. Because the
# package separates the world from the hypothesis, the search simply
# re-runs the inversion (set_ape()/set_aie()) at each candidate target --
# so every route (parametric, empirical, pilot-model, panel, IV) and both
# estimands work unchanged, and n keeps its route meaning (units for
# panel designs).

# Re-pin a DGP at a candidate target, whatever the estimand.
repin <- function(dgp, target, mf = NULL, mm = NULL) {
  if (identical(dgp$estimand, "aie") || !is.null(dgp$moderator)) {
    set_aie(dgp, target, main_focal = mf, main_moderator = mm)
  } else {
    set_ape(dgp, target)
  }
}

#' Minimum detectable APE/AIE at a given sample size
#'
#' The inverse of [ape_n()]: the sample size is fixed -- an archive with
#' so many firms, a panel with so many units and waves -- and the
#' question is the smallest effect the design can reliably conclude.
#' For `claim = "detect"` and `claim = "minimum"` the function searches
#' over the assumed true effect, re-running the APE/AIE inversion at
#' every candidate, until simulated power at `n` hits the goal; for
#' `claim = "equivalence"` the target is held at the DGP's pinned value
#' (typically 0) and the search is over the margin instead, returning
#' the tightest equivalence bounds the design can expect to establish.
#'
#' By default the answer is verified the way [ape_n()] verifies its n: a
#' high-precision confirmation stage re-measures power at the candidate
#' and pushes the effect upward (never downward) if it falls short, so
#' the reported MDE errs on the conservative side.
#'
#' For a DGP with a moderator the searched effect is the AIE; the two
#' conditional-at-reference main-effect anchors are taken from the
#' pinned DGP (or passed via `main_focal`/`main_moderator`) and held
#' fixed across candidates.
#'
#' @inheritParams ape_power
#' @param n Sample size of the planned study (units for panel designs).
#' @param power Target power for the claim (default 0.80).
#' @param main_focal,main_moderator Main-effect anchors for AIE designs;
#'   defaults to the values stored by a previous [set_aie()] call.
#' @param max_iter Maximum search refinements.
#' @param confirm,nsim_confirm Confirmation stage as in [ape_n()].
#'
#' @return A `powerape_mde` object: `mde` (the minimum detectable effect,
#'   or for equivalence the smallest establishable margin), the confirmed
#'   `power` and `mcse` at that effect, the search `history`, and the DGP
#'   re-pinned at the answer (so the object feeds [power_statement()]).
#' @examples
#' \donttest{
#' d <- ape_dgp(focal = pa_var("treat", "binary", p = 0.5), baseline = 0.30)
#' ape_mde(d, n = 712, claim = "detect", nsim = 600, seed = 1)
#' }
#' @export
ape_mde <- function(dgp, n, power = 0.80,
                    claim = c("detect", "minimum", "equivalence"),
                    sesoi = NULL, conf = 0.95, nsim = 1000, seed = NULL,
                    main_focal = NULL, main_moderator = NULL,
                    max_iter = 6, confirm = TRUE, nsim_confirm = 4 * nsim,
                    se = c("model", "robust")) {
  claim <- match.arg(claim)
  se <- match.arg(se)
  stopifnot(inherits(dgp, "powerape_dgp"),
            is.numeric(n), length(n) == 1L, n >= 20,
            is.numeric(power), length(power) == 1L, power > 0.5, power < 0.999,
            is.numeric(conf), length(conf) == 1L, conf > 0.5, conf < 1,
            is.numeric(nsim), length(nsim) == 1L, nsim >= 20)
  if (claim == "minimum" &&
      (is.null(sesoi) || !is.numeric(sesoi) || length(sesoi) != 1L || sesoi <= 0))
    stop("claim = \"minimum\" needs a single positive `sesoi` (in APE units).",
         call. = FALSE)
  is_aie <- !is.null(dgp$moderator)
  mf <- main_focal %||% dgp$main_focal
  mm <- main_moderator %||% dgp$main_moderator
  if (is_aie && (is.null(mf) || is.null(mm)))
    stop(paste("An AIE MDE needs the main-effect anchors: pin the DGP with",
               "set_aie() first, or pass `main_focal` and `main_moderator`."),
         call. = FALSE)
  est_lab <- if (is_aie) "AIE" else "APE"
  z <- zcrit(conf)
  zg <- qnorm(power)

  ## --- equivalence: search the margin at the pinned target ------------------
  if (claim == "equivalence") {
    if (is.null(dgp$target_est))
      stop(paste("An equivalence MDE holds the truth fixed: pin the DGP first",
                 "(typically at target = 0) with set_ape() or set_aie()."),
           call. = FALSE)
    t_abs <- abs(dgp$target_est)
    ps <- sim_ci(dgp, n, 300L, conf,
                 seed = if (is.null(seed)) NULL else seed - 1L,
                 se_type = if (dgp$route %in% c("panel", "iv")) "model" else se)
    if (mean(ps$ok) < 0.5)
      stop("Most pilot replications failed to fit; the design looks degenerate.",
           call. = FALSE)
    se_hat <- mean(ps$se[ps$ok])
    d_i <- t_abs + se_hat * (z + zg)
    hist <- data.frame(stage = character(), iter = integer(),
                       value = numeric(), power = numeric(), mcse = numeric())
    best <- NULL
    for (it in seq_len(max_iter)) {
      pw <- power_once(dgp, n, "equivalence", d_i, conf, nsim,
                       if (is.null(seed)) NULL else seed + it,
                       enforce = FALSE, se = se)
      hist <- rbind(hist, data.frame(stage = "search", iter = it, value = d_i,
                                     power = pw$power, mcse = pw$mcse))
      if (is.null(best) || abs(pw$power - power) < abs(best$power - power))
        best <- list(value = d_i, power = pw$power, mcse = pw$mcse)
      if (abs(pw$power - power) <= max(0.01, 2 * pw$mcse)) break
      p_cl <- min(max(pw$power, 0.02), 0.998)
      se_impl <- (d_i - t_abs) / max(z + qnorm(p_cl), 0.1)
      d_i <- t_abs + se_impl * (z + zg)
    }
    res <- finish_mde(dgp, n, best, hist, power, "equivalence", best$value,
                      conf, nsim, seed, confirm, nsim_confirm, se,
                      mf, mm, is_margin = TRUE, t_abs = t_abs)
    return(res)
  }

  ## --- detect / minimum: search the target ----------------------------------
  delta0 <- if (claim == "minimum") sesoi else 0

  ## pilot SE at a reference effect (reuse an existing pin when present)
  d_ref <- if (!is.null(dgp$target_est)) {
    dgp
  } else {
    t_try <- c(0.5 * min(dgp$baseline, 1 - dgp$baseline), 0.05, 0.02)
    d0 <- NULL
    for (tt in t_try) {
      d0 <- tryCatch(repin(dgp, tt, mf, mm), error = function(e) NULL)
      if (!is.null(d0)) break
    }
    if (is.null(d0))
      stop("Could not pin a reference effect for the pilot stage.", call. = FALSE)
    d0
  }
  ps <- sim_ci(d_ref, n, 300L, conf,
               seed = if (is.null(seed)) NULL else seed - 1L,
               se_type = if (dgp$route %in% c("panel", "iv")) "model" else se)
  if (mean(ps$ok) < 0.5)
    stop("Most pilot replications failed to fit; the design looks degenerate.",
         call. = FALSE)
  se_hat <- mean(ps$se[ps$ok])

  t_i <- delta0 + se_hat * (z + zg)
  hist <- data.frame(stage = character(), iter = integer(),
                     value = numeric(), power = numeric(), mcse = numeric())
  best <- NULL
  eval_at <- function(t, nsim_use, seed_use) {
    d_t <- tryCatch(repin(dgp, t, mf, mm), error = function(e) e)
    if (inherits(d_t, "error")) return(d_t)
    power_once(d_t, n, claim, sesoi, conf, nsim_use, seed_use,
               enforce = FALSE, se = se)
  }
  for (it in seq_len(max_iter)) {
    pw <- eval_at(t_i, nsim, if (is.null(seed)) NULL else seed + it)
    if (inherits(pw, "error")) {
      ## candidate beyond the feasible range: report the ceiling honestly
      ceil <- tryCatch({
        f <- function(t) inherits(tryCatch(repin(dgp, t, mf, mm),
                                           error = function(e) e), "error")
        lo <- delta0 + 1e-4; hi <- t_i
        for (k in 1:30) { mid <- (lo + hi) / 2; if (f(mid)) hi <- mid else lo <- mid }
        lo
      }, error = function(e) NA_real_)
      p_ceil <- if (is.finite(ceil)) {
        pc <- eval_at(0.98 * ceil, nsim, if (is.null(seed)) NULL else seed + 50L)
        if (inherits(pc, "error")) NA_real_ else pc$power
      } else NA_real_
      stop(sprintf(paste(
        "No detectable %s at %.0f%% power within the feasible range:",
        "the largest attainable effect is about %.3f, where simulated power",
        "reaches only %.2f. Increase n (or loosen the goal)."),
        est_lab, 100 * power, ceil, p_ceil), call. = FALSE)
    }
    hist <- rbind(hist, data.frame(stage = "search", iter = it, value = t_i,
                                   power = pw$power, mcse = pw$mcse))
    if (is.null(best) || abs(pw$power - power) < abs(best$power - power))
      best <- list(value = t_i, power = pw$power, mcse = pw$mcse)
    if (abs(pw$power - power) <= max(0.01, 2 * pw$mcse)) break
    p_cl <- min(max(pw$power, 0.02), 0.998)
    se_impl <- (t_i - delta0) / max(z + qnorm(p_cl), 0.1)
    t_i <- delta0 + se_impl * (z + zg)
  }

  finish_mde(dgp, n, best, hist, power, claim, sesoi, conf, nsim, seed,
             confirm, nsim_confirm, se, mf, mm,
             is_margin = FALSE, delta0 = delta0)
}

# Confirmation stage + object assembly, shared by both branches. Accepts
# only if confirmed power is within tolerance below the goal; otherwise
# pushes the effect (or margin) UP -- the conservative direction.
finish_mde <- function(dgp, n, best, hist, power, claim, sesoi, conf, nsim,
                       seed, confirm, nsim_confirm, se, mf, mm,
                       is_margin = FALSE, delta0 = 0, t_abs = 0) {
  z <- zcrit(conf)
  zg <- qnorm(power)
  anchor <- if (is_margin) t_abs else delta0
  confirmed <- FALSE
  if (confirm) {
    stopifnot(is.numeric(nsim_confirm), length(nsim_confirm) == 1L,
              nsim_confirm >= nsim)
    v_c <- best$value
    for (cr in seq_len(3L)) {
      sc <- if (is.null(seed)) NULL else seed + 100L + cr
      pwc <- if (is_margin) {
        power_once(dgp, n, "equivalence", v_c, conf, nsim_confirm, sc,
                   enforce = FALSE, se = se)
      } else {
        d_t <- tryCatch(repin(dgp, v_c, mf, mm), error = function(e) e)
        if (inherits(d_t, "error")) break
        power_once(d_t, n, claim, sesoi, conf, nsim_confirm, sc,
                   enforce = FALSE, se = se)
      }
      hist <- rbind(hist, data.frame(stage = "confirm", iter = cr, value = v_c,
                                     power = pwc$power, mcse = pwc$mcse))
      best <- list(value = v_c, power = pwc$power, mcse = pwc$mcse)
      if (pwc$power >= power - max(0.005, 1.5 * pwc$mcse)) {
        confirmed <- TRUE
        break
      }
      if (cr == 3L) break
      p_cl <- min(max(pwc$power, 0.02), 0.998)
      se_impl <- (v_c - anchor) / max(z + qnorm(p_cl), 0.1)
      v_c <- anchor + se_impl * (z + zg)
    }
  }

  dgp_at <- if (is_margin) dgp else repin(dgp, best$value, mf, mm)
  structure(list(
    mde = best$value, power = best$power, mcse = best$mcse,
    goal = power, claim = claim, sesoi = if (is_margin) best$value else sesoi,
    conf = conf, n = as.integer(n), nsim = as.integer(nsim),
    nsim_confirm = as.integer(nsim_confirm),
    confirmed = confirmed, confirm_requested = isTRUE(confirm),
    is_margin = is_margin, se = se,
    estimand = if (!is.null(dgp$moderator)) "aie" else "ape",
    model = dgp$model, history = hist, dgp = dgp_at
  ), class = "powerape_mde")
}

#' @export
print.powerape_mde <- function(x, ...) {
  unit_lab <- if (identical(x$dgp$route, "panel")) {
    sprintf("n = %d units (x %d periods)", x$n, x$dgp$n_periods)
  } else {
    sprintf("n = %d", x$n)
  }
  if (x$is_margin) {
    cat(sprintf("powerape minimum establishable equivalence margin -- %s\n",
                toupper(x$estimand)))
    cat(sprintf("  bounds +/-%.4f at %s for %.0f%% target power (%s %.3f, MCSE %.3f)\n",
                x$mde, unit_lab, 100 * x$goal,
                if (isTRUE(x$confirmed)) "confirmed" else "achieved",
                x$power, x$mcse))
  } else {
    cat(sprintf("powerape minimum detectable %s -- %s claim\n",
                toupper(x$estimand), x$claim))
    cat(sprintf("  MDE = %.4f at %s for %.0f%% target power (%s %.3f, MCSE %.3f)\n",
                x$mde, unit_lab, 100 * x$goal,
                if (isTRUE(x$confirmed)) "confirmed" else "achieved",
                x$power, x$mcse))
    if (x$claim == "minimum")
      cat(sprintf("  smallest effect demonstrably above sesoi %.3f\n", x$sesoi))
  }
  n_search <- sum(x$history$stage == "search")
  n_conf <- sum(x$history$stage == "confirm")
  if (isTRUE(x$confirm_requested)) {
    if (isTRUE(x$confirmed)) {
      cat(sprintf("  search: %d step(s); confirmed in %d round(s) at nsim = %d.\n",
                  n_search, n_conf, x$nsim_confirm))
    } else {
      cat(sprintf("  search: %d step(s); confirmation (%d round(s), nsim = %d) fell short -- treat as a lower bound.\n",
                  n_search, n_conf, x$nsim_confirm))
    }
  } else {
    cat(sprintf("  search: %d step(s); confirm with ape_power() at a larger nsim.\n",
                n_search))
  }
  invisible(x)
}
