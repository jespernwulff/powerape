# Justified alpha -------------------------------------------------------------
# Two routes to an error rate that is chosen rather than inherited (DESIGN.md
# section 15). (1) A sample-size rule: `alpha` may be a function of n in
# every public function (for example alphaN::alphaN), evaluated at each
# candidate n, so ape_n()/ape_curve()/ape_mde()/ape_robust() design jointly
# over (alpha(n), n). (2) Error-cost optimization on stored draws: the engine
# keeps the estimate and SE of every replication, so power at any alpha is a
# re-threshold of ONE simulation -- exact on those draws, monotone in alpha,
# free of Monte Carlo jitter across alpha -- and ape_alpha() minimizes or
# balances the weighted combined error rate of Maier & Lakens (2022) over
# that curve, capped at the conventional .05.

# Re-threshold stored draws at a vector of alphas.
power_at_draws <- function(draws, target, claim, sesoi, nsim, alpha) {
  rows <- lapply(alpha, function(a) {
    r <- summarize_sim(draws, target, claim, sesoi, nsim, a)
    c(alpha = a, conf = claim_levels(claim, a)$conf_claim,
      power = r$power, mcse = r$mcse, r$outcomes)
  })
  out <- as.data.frame(do.call(rbind, rows))
  rownames(out) <- NULL
  out
}

# The world in which the claim sits exactly on its null boundary: true
# effect = +/-SESOI for the SESOI claims (in the hypothesized direction),
# 0 for detection. Used to measure the REALIZED size at a chosen alpha.
boundary_dgp <- function(dgp, claim, sesoi, target) {
  sgn <- if (target < 0) -1 else 1
  b <- switch(claim, minimum = sgn * sesoi, equivalence = sgn * sesoi,
              detect = 0)
  if (identical(dgp$estimand, "aie")) {
    set_aie(dgp, b, main_focal = dgp$main_focal,
            main_moderator = dgp$main_moderator)
  } else {
    set_ape(dgp, b)
  }
}

#' Power at other error rates from one stored simulation
#'
#' [ape_power()] keeps the estimate and standard error of every replication
#' (`keep_draws = TRUE`, the default), so the claim can be re-evaluated at
#' any error rate without simulating again: power at `alpha` is a
#' re-threshold of the same draws, exact on those draws and monotone in
#' `alpha`. This is the curve that [ape_alpha()] optimizes.
#'
#' @param x A `powerape_power` object from [ape_power()] with stored draws.
#' @param alpha Numeric vector of error rates in (0, 0.5), each read in the
#'   claim's conventional form (two-sided for detection, one-sided for the
#'   minimum-effect claim, TOST for equivalence; see [ape_power()]).
#' @return A data frame with one row per `alpha`: `alpha`, the interval
#'   level `conf` it implies for the claim, `power`, its Monte Carlo
#'   standard error `mcse`, and the outcome distribution.
#' @examples
#' \donttest{
#' d <- set_ape(ape_dgp(focal = pa_var("treat", "binary", p = 0.5),
#'                      baseline = 0.30), 0.10)
#' pw <- ape_power(d, n = 800, claim = "minimum", sesoi = 0.03,
#'                 nsim = 500, seed = 1)
#' power_at(pw, alpha = c(0.01, 0.025, 0.05, 0.10))
#' }
#' @export
power_at <- function(x, alpha) {
  if (!inherits(x, "powerape_power"))
    stop("`x` must be an ape_power() result.", call. = FALSE)
  if (is.null(x$draws))
    stop("`x` carries no stored draws; rerun ape_power() with keep_draws = TRUE.",
         call. = FALSE)
  stopifnot(is.numeric(alpha), length(alpha) >= 1L, all(is.finite(alpha)),
            all(alpha > 0), all(alpha < 0.5))
  power_at_draws(x$draws, x$target, x$claim, x$sesoi, x$nsim, alpha)
}

#' Justify the error rate: minimize or balance weighted Type I and Type II errors
#'
#' Chooses the claim's error rate by the decision-theoretic route of Mudge
#' et al. (2012) and Maier and Lakens (2022): for the design in `x` (its
#' world, claim, planning value, and sample size), the weighted combined
#' error rate
#' \deqn{w(\alpha) = \frac{cost \cdot \alpha + prior \cdot \beta(\alpha)}{cost + prior}}
#' is evaluated on a fine grid of `alpha` from the stored draws of one
#' simulation (see [power_at()]) and either minimized (`error =
#' "minimize"`) or balanced so that `cost * alpha = prior * beta` (`error =
#' "balance"`). Here `cost` is the relative cost of a Type I error against
#' a Type II error (Cohen's 4:1 convention produces the familiar .05/.20
#' pair at 80% power) and `prior` the prior odds that the planning value
#' rather than the claim boundary is true. For the SESOI claims the Type I
#' error is the size at the claim's boundary (true effect = SESOI for the
#' minimum-effect claim, = +/-SESOI for equivalence), for detection the
#' directional size at zero; the Type II error is one minus power at the
#' planning value.
#'
#' **The cap.** The search is confined to `alpha <= cap`, .05 by default.
#' Maier and Lakens (2022) argue that raising the error rate above the
#' conventional level requires a justification of its own (direct
#' decision use, a cost--benefit rationale, a low prior on the null, no
#' room to add data) and that journals may accept only lowering; the
#' function therefore reports the capped optimum and states when the
#' unconstrained optimum lies above the cap.
#'
#' **Precision.** Power along the grid is a step function of the stored
#' draws, so the objective is flat near its optimum. The result reports,
#' besides the optimum, the range of error rates whose objective lies
#' within one Monte Carlo standard error of it; report that range, not
#' the point, unless `nsim` is large.
#'
#' **Realized size.** With `size = TRUE` the function simulates the world
#' in which the claim sits exactly on its null boundary and reports the
#' realized error rate at the chosen `alpha`, which the analytic tools
#' behind compromise power analysis cannot do. On the standard routes it
#' matches the nominal level within Monte Carlo error; on clustered panels
#' with few units it may not, and then the realized number is the one to
#' report.
#'
#' @param x A `powerape_power` object from [ape_power()] with stored draws.
#' @param cost Relative cost of a Type I error against a Type II error
#'   (default 1: equally costly).
#' @param prior Prior odds of the planning value against the claim
#'   boundary (default 1: equally likely).
#' @param error `"minimize"` (default) the weighted combined error rate, or
#'   `"balance"` the cost-weighted Type I and Type II error rates.
#' @param cap Largest admissible `alpha` (default 0.05). Values above .05
#'   warn.
#' @param size Simulate the claim-boundary world to measure the realized
#'   error rate at the chosen `alpha` (default FALSE).
#' @param nsim_size Replications for the size simulation (default: the
#'   `nsim` of `x`).
#' @param seed Optional seed for the size simulation.
#'
#' @return A `powerape_alpha` object: `alpha`, the chosen error rate;
#'   `power`, `mcse`, and `beta` at it; `wcer`, the weighted combined error
#'   rate there and `wcer_default` at alpha = .05; `range`, the flat
#'   range; `capped`, whether the unconstrained optimum (`alpha_uncapped`,
#'   `wcer_uncapped`) lies above the cap; `size`, the realized-size
#'   simulation when requested; `grid`, the full alpha/power/beta/wcer
#'   curve; and the embedded DGP. Has `print`, `plot`, and
#'   [power_statement()] methods.
#' @references Maier, M., & Lakens, D. (2022). Justify your alpha: A primer
#'   on two practical approaches. *Advances in Methods and Practices in
#'   Psychological Science, 5*(2). Mudge, J. F., Baker, L. F., Edge, C. B.,
#'   & Houlahan, J. E. (2012). Setting an optimal alpha that minimizes
#'   errors in null hypothesis significance tests. *PLOS ONE, 7*(2).
#' @examples
#' \donttest{
#' d <- set_ape(ape_dgp(focal = pa_var("treat", "binary", p = 0.5),
#'                      covariates = list(pa_var("z", "normal")),
#'                      baseline = 0.30, signal = 0.10), 0.10)
#' pw <- ape_power(d, n = 2000, claim = "minimum", sesoi = 0.05,
#'                 nsim = 1000, seed = 1)
#' ja <- ape_alpha(pw, cost = 4, prior = 1)   # Cohen's 4:1 weighting
#' ja
#' plot(ja)
#' power_statement(ja)
#' }
#' @export
ape_alpha <- function(x, cost = 1, prior = 1, error = c("minimize", "balance"),
                      cap = 0.05, size = FALSE, nsim_size = NULL, seed = NULL) {
  error <- match.arg(error)
  if (!inherits(x, "powerape_power"))
    stop("`x` must be an ape_power() result.", call. = FALSE)
  if (is.null(x$draws))
    stop("`x` carries no stored draws; rerun ape_power() with keep_draws = TRUE.",
         call. = FALSE)
  stopifnot(is.numeric(cost), length(cost) == 1L, is.finite(cost), cost > 0,
            is.numeric(prior), length(prior) == 1L, is.finite(prior), prior > 0,
            is.numeric(cap), length(cap) == 1L, cap > 0, cap < 0.5)
  if (cap > 0.05 + 1e-9)
    warning(paste("`cap` exceeds .05: an error rate above the conventional",
                  "level needs a justification of its own (Maier & Lakens, 2022)."),
            call. = FALSE)

  grid <- c(seq(1e-4, 0.01, by = 1e-4), seq(0.0105, 0.4995, by = 5e-4))
  pa <- power_at_draws(x$draws, x$target, x$claim, x$sesoi, x$nsim, grid)
  pw <- pa$power
  beta <- 1 - pw
  wcer <- (cost * grid + prior * beta) / (cost + prior)
  obj <- if (error == "minimize") wcer else abs(cost * grid - prior * beta)
  i_free <- which.min(obj)
  in_cap <- grid <= cap + 1e-9
  idx_cap <- which(in_cap)
  i_cap <- idx_cap[which.min(obj[idx_cap])]
  capped <- grid[i_free] > cap + 1e-9
  a_star <- grid[i_cap]
  p_star <- pw[i_cap]
  mcse_star <- sqrt(p_star * (1 - p_star) / x$nsim)
  ## error rates (within the cap) whose objective lies within one Monte
  ## Carlo SE of the optimum's -- alpha itself is nominal, so the only
  ## noise is in beta
  tol <- if (error == "minimize") prior * mcse_star / (cost + prior) else prior * mcse_star
  flat <- in_cap & (obj <= obj[i_cap] + tol)
  i_def <- which.min(abs(grid - 0.05))

  size_res <- NULL
  if (isTRUE(size)) {
    ns <- as.integer(nsim_size %||% x$nsim)
    stopifnot(ns >= 20L)
    d_b <- boundary_dgp(x$dgp, x$claim, x$sesoi, x$target)
    sgn <- if (x$target < 0) -1 else 1
    se_type <- if (x$dgp$route %in% c("panel", "iv")) "model" else (x$se %||% "model")
    sim <- sim_ci(d_b, x$n, ns, claim_levels(x$claim, a_star)$conf_claim, seed,
                  se_type = se_type, separation = x$separation %||% "fail")
    ## summarize_sim() reads only the sign of its target argument; pass the
    ## hypothesized direction so the boundary world is scored the same way
    r <- summarize_sim(sim, sgn, x$claim, x$sesoi, ns, a_star)
    size_res <- list(size = r$power, mcse = sqrt(r$power * (1 - r$power) / ns),
                     nsim = ns, truth = d_b$target_est, n_failed = r$n_failed,
                     separated = r$separated)
  }

  structure(list(
    alpha = a_star, beta = beta[i_cap], power = p_star, mcse = mcse_star,
    wcer = wcer[i_cap], power_default = pw[i_def], wcer_default = wcer[i_def],
    range = range(grid[flat]), capped = capped, cap = cap,
    alpha_uncapped = grid[i_free], beta_uncapped = beta[i_free],
    wcer_uncapped = wcer[i_free],
    cost = cost, prior = prior, error = error,
    claim = x$claim, sesoi = x$sesoi, n = x$n, nsim = x$nsim,
    target = x$target, estimand = x$estimand %||% "ape", model = x$model,
    se = x$se, size = size_res,
    grid = data.frame(alpha = grid, power = pw, beta = beta, wcer = wcer),
    dgp = x$dgp
  ), class = "powerape_alpha")
}

#' @export
print.powerape_alpha <- function(x, ...) {
  claim_lab <- switch(x$claim,
    minimum = sprintf("minimum-effect claim (CI lower bound > %.3f)", x$sesoi),
    detect = "detection claim (CI excludes 0, directional)",
    equivalence = sprintf("equivalence claim (CI within +/-%.3f)", x$sesoi))
  cat(sprintf("powerape justified alpha -- %s, n = %d, nsim = %d\n",
              claim_lab, x$n, x$nsim))
  cat(sprintf("  objective: %s (Type I cost x%s, prior odds H1:H0 = %s)\n",
              if (x$error == "minimize") "minimize the weighted combined error rate"
              else "balance the cost-weighted Type I and Type II error rates",
              format(x$cost), format(x$prior)))
  cat(sprintf("  alpha* = %.4f -> power %.3f (MCSE %.3f), beta %.3f, weighted combined error %.4f\n",
              x$alpha, x$power, x$mcse, x$beta, x$wcer))
  cat(sprintf("  at alpha = 0.05: power %.3f, weighted combined error %.4f\n",
              x$power_default, x$wcer_default))
  cat(sprintf("  flat range: alpha in [%.4f, %.4f] is indistinguishable from the optimum at this nsim\n",
              x$range[1], x$range[2]))
  if (isTRUE(x$capped)) {
    cat(sprintf("  cap %.3f BINDS: unconstrained optimum alpha = %.4f (weighted combined error %.4f);\n",
                x$cap, x$alpha_uncapped, x$wcer_uncapped))
    cat("    raising alpha above the conventional level needs its own justification (Maier & Lakens, 2022).\n")
  } else {
    cat(sprintf("  cap %.3f not binding (unconstrained optimum %.4f)\n",
                x$cap, x$alpha_uncapped))
  }
  if (!is.null(x$size))
    cat(sprintf("  realized size at the claim boundary (true %s = %.3f): %.4f (MCSE %.4f, nsim %d)\n",
                toupper(x$estimand), x$size$truth, x$size$size, x$size$mcse,
                x$size$nsim))
  invisible(x)
}

#' @param xlim Horizontal range for the plot (default: 0 to the larger of
#'   0.2 and three times the cap).
#' @param ... Passed to the base plot call.
#' @rdname ape_alpha
#' @export
plot.powerape_alpha <- function(x, xlim = NULL, ...) {
  g <- x$grid
  xl <- xlim %||% c(0, max(0.2, min(0.5, 3 * x$cap)))
  keep <- g$alpha <= xl[2]
  plot(g$alpha[keep], g$wcer[keep], type = "l", lwd = 1.5,
       xlab = "alpha", ylab = "Weighted combined error rate",
       xlim = xl, ylim = c(0, max(g$wcer[keep])), ...)
  abline(h = seq(0, 1, 0.1), col = "grey92", lwd = 0.8)
  abline(v = 0.05, lty = 2, col = "grey60")
  abline(v = x$cap, lty = 3, col = "grey40")
  if (isTRUE(x$capped)) {
    abline(v = x$alpha_uncapped, lty = 3, col = "grey70")
    points(x$alpha_uncapped, x$wcer_uncapped, pch = 1)
  }
  points(x$alpha, x$wcer, pch = 19)
  mtext(sprintf("alpha* = %.4f (%s), weighted combined error %.4f",
                x$alpha, if (isTRUE(x$capped)) "cap binds" else "cap not binding",
                x$wcer),
        side = 3, line = 0.2, cex = 0.85)
  invisible(x)
}

#' Closed-form justified alpha at the detectability frontier
#'
#' The normal-approximation companion of [ape_alpha()] for pricing
#' minimum detectable effects: the error rate that minimizes the
#' weighted combined Type I and Type II error rate (Mudge et al., 2012;
#' Maier & Lakens, 2022) for the effect that is *just detectable* at the
#' target power. Under the normal approximation, the minimize-mode
#' first-order condition for an effect-to-standard-error ratio `r` is
#' `m * cost * phi(z) = prior * phi(z - r)`, with `z` the claim's
#' critical value and `m = 2` for the two-sided detection claim (`m = 1`
#' for the one-sided minimum-effect claim). At the detectability
#' frontier the ratio `r` equals `z + qnorm(power)`, so the condition
#' collapses to
#' `phi(z) = prior * phi(qnorm(power)) / (m * cost)` -- one equation
#' involving neither the sample size nor the standard error. The
#' frontier level is therefore the **same at every n**, which makes it
#' the natural companion of an MDE grid: solve the level once, pass it
#' as `alpha` to [ape_mde()] (or [ape_power()], [ape_n()]), and the
#' whole grid is priced under one justified error rate.
#'
#' Two anchors worth knowing. At Cohen's 4:1 weighting and 80% power the
#' two-sided frontier level is 0.0274. At *equal* costs the one-sided
#' frontier level is exactly `1 - power` (the frontier balances
#' `alpha = beta` there).
#'
#' **Relation to [ape_alpha()].** `ape_alpha()` is exact for a specific
#' design: it optimizes over the stored draws of one simulation, applies
#' the conventional `cap = 0.05`, and can verify the realized size.
#' `alpha_frontier()` is its closed-form, scale-free limit at the
#' just-detectable effect; it applies no cap (levels above .05 need a
#' justification of their own; Maier & Lakens, 2022) and warns when the
#' optimum does not lie below 0.5, where it could not be used as an
#' `alpha` anyway.
#'
#' @param cost Relative cost of a Type I error against a Type II error
#'   (default 1).
#' @param power Target power defining the frontier effect (default 0.80).
#' @param prior Prior odds of the just-detectable effect against the
#'   claim boundary (default 1).
#' @param claim `"detect"` (default; two-sided) or `"minimum"`
#'   (one-sided). For equivalence run [ape_alpha()] on stored draws.
#' @return The frontier error rate (a numeric scalar), suitable as the
#'   `alpha` argument of the power functions when it is below 0.5.
#' @references Maier, M., & Lakens, D. (2022). Justify your alpha: A
#'   primer on two practical approaches. *Advances in Methods and
#'   Practices in Psychological Science, 5*(2). Mudge, J. F., Baker,
#'   L. F., Edge, C. B., & Houlahan, J. E. (2012). Setting an optimal
#'   alpha that minimizes errors in null hypothesis significance tests.
#'   *PLOS ONE, 7*(2).
#' @examples
#' alpha_frontier(cost = 4, power = 0.80)            # 0.0274
#' alpha_frontier(cost = 1, power = 0.80, claim = "minimum")  # exactly 0.20
#' @export
alpha_frontier <- function(cost = 1, power = 0.80, prior = 1,
                           claim = c("detect", "minimum")) {
  claim <- match.arg(claim)
  stopifnot(is.numeric(cost), length(cost) == 1L, is.finite(cost), cost > 0,
            is.numeric(prior), length(prior) == 1L, is.finite(prior), prior > 0,
            is.numeric(power), length(power) == 1L, power > 0, power < 1)
  m <- if (claim == "detect") 2 else 1
  rhs <- prior * dnorm(qnorm(power)) / (m * cost)
  if (rhs >= dnorm(0))
    stop(sprintf(paste(
      "No interior optimum at the frontier: with cost = %s and prior = %s,",
      "the weighted combined error rate keeps falling as alpha rises toward",
      "1 for the just-detectable effect, so no level below 1 is optimal.",
      "Weight Type I errors more heavily (larger cost / smaller prior), or",
      "run ape_alpha() on a stored simulation, which reports the capped",
      "optimum for a specific design."),
      format(cost), format(prior)), call. = FALSE)
  a <- m * pnorm(-sqrt(-2 * log(sqrt(2 * pi) * rhs)))
  if (a >= 0.5)
    warning(sprintf(paste(
      "The frontier level (%.3f) is not below 0.5 and cannot be used as an",
      "`alpha` in powerape's power functions; this cost/prior weighting does",
      "not support a conventional test."), a), call. = FALSE)
  a
}

# The justification paragraph for power_statement().
alpha_statement <- function(x) {
  d <- x$dgp
  est_short <- if (identical(x$estimand, "aie")) "AIE" else "APE"
  how <- if (x$error == "minimize") {
    "minimizing the weighted combined Type I and Type II error rate"
  } else {
    "balancing the cost-weighted Type I and Type II error rates"
  }
  n_text <- if (identical(d$route, "panel")) sprintf("n = %d units", x$n) else sprintf("n = %d", x$n)
  s1 <- sprintf(paste0("We chose the error rate for %s by %s (Maier & Lakens, ",
                       "2022; Mudge et al., 2012), weighting a Type I error %s ",
                       "times a Type II error and taking prior odds of %s for ",
                       "the planning value (%s = %.3f) against the claim ",
                       "boundary, using the powerape package (version %s)."),
                claim_text(x$claim, x$sesoi, claim_levels(x$claim, x$alpha)$conf_claim),
                how, format(x$cost), format(x$prior), est_short, x$target,
                as.character(packageVersion("powerape")))
  s2 <- sprintf("The assumed data-generating process was %s.", describe_dgp(d))
  s3 <- sprintf(paste0("At %s the optimum within the conventional cap of alpha ",
                       "= %s was alpha = %.4f: simulated power %.3f (Monte ",
                       "Carlo SE %.3f; %d replications), a Type II error rate ",
                       "of %.3f, and a weighted combined error rate of %.4f ",
                       "versus %.4f at alpha = 0.05; error rates between %.4f ",
                       "and %.4f are indistinguishable from the optimum at ",
                       "this precision."),
                n_text, format(x$cap), x$alpha, x$power, x$mcse, x$nsim, x$beta,
                x$wcer, x$wcer_default, x$range[1], x$range[2])
  s4 <- if (isTRUE(x$capped)) {
    sprintf(paste0("The unconstrained optimum, alpha = %.4f (weighted combined ",
                   "error rate %.4f), lies above the cap; because raising the ",
                   "error rate above the conventional level needs a ",
                   "justification of its own (Maier & Lakens, 2022), the capped ",
                   "value is reported."), x$alpha_uncapped, x$wcer_uncapped)
  } else {
    NULL
  }
  s5 <- if (!is.null(x$size)) {
    sprintf(paste0("A verification simulation at the claim boundary (true %s = ",
                   "%.3f; %d replications) gave a realized error rate of %.4f ",
                   "(Monte Carlo SE %.4f) at the chosen alpha."),
            est_short, x$size$truth, x$size$nsim, x$size$size, x$size$mcse)
  } else {
    NULL
  }
  structure(paste(c(s1, s2, s3, s4, s5), collapse = " "),
            class = "powerape_statement")
}
