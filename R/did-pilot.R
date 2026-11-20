# A difference-in-differences design from pre-treatment pilot data -------------
# Data observed before the treatment date inform four
# inputs of ape_dgp_did(): the control group's outcome rate (baseline), the
# gap between the groups, the treated share, and, in a panel, the latent
# within-unit correlation of the route's Gaussian copula. They do not
# inform the control group's change from before to after (`trend`), which
# the user supplies.

# P(Y_s = 1, Y_t = 1) for two periods of one unit when both have marginal
# rate p and the latent correlation is rho >= 0 (the route's copula): the
# integral over the shared unit draw of the squared conditional rate.
did_joint_rate <- function(p, rho) {
  if (rho <= 0) return(p^2)
  q <- qnorm(p)
  f <- function(a) dnorm(a) * pnorm((q - sqrt(rho) * a) / sqrt(1 - rho))^2
  integrate(f, -Inf, Inf, rel.tol = 1e-10, abs.tol = 1e-12)$value
}

# Latent within-unit correlation by the method of moments: the rho at which
# the copula's joint rate matches the observed share of within-unit pairs of
# periods with the outcome in both, pooled over the groups with the number
# of pairs as weights. `p`, `m`, `w`: per group, the marginal rate, the
# observed joint rate and the number of pairs. Returns 0 when the data show
# no positive dependence.
did_pilot_rho <- function(p, m, w) {
  f <- function(rho) sum(w * (vapply(p, did_joint_rate, numeric(1), rho = rho) - m))
  if (f(0) >= 0) return(0)
  if (f(0.995) <= 0) return(0.995)
  uniroot(f, c(0, 0.995), tol = 1e-9)$root
}

# The same moment condition for a design with covariates. Covariates are
# constant over time, so they are a second source of dependence between a
# unit's outcomes, and the design's `rho` is the latent correlation that is
# left GIVEN the covariates. The raw pairs of the pilot show the whole
# dependence, so the rho to use is the one at which the design's joint rate,
# averaged over the covariate index q of each group, equals the pilot's.

# An equally weighted grid standing in for the draws of q (quantile
# midpoints): the averages below are smooth in q, and 2,000 points
# reproduce them to about 1e-6.
did_q_grid <- function(q, K = 2000L) {
  if (length(q) <= K) return(q)
  stats::quantile(q, (seq_len(K) - 0.5) / K, names = FALSE, type = 7)
}

# P(Y_s = 1, Y_t = 1) for two periods of one unit whose rate is
# G(cc + q), averaged over the values `q`, at latent correlation `rho`
# given q.
did_joint_rate_cov <- function(lf, cc, q, rho) {
  p <- lf$G(cc + q)
  if (rho <= 0) return(mean(p^2))
  thr <- qnorm(p)
  f <- function(a) {
    dnorm(a) * vapply(a, function(ai)
      mean(pnorm((thr - sqrt(rho) * ai) / sqrt(1 - rho))^2), numeric(1))
  }
  integrate(f, -Inf, Inf, rel.tol = 1e-9, abs.tol = 1e-12)$value
}

# rho given the covariates. `p`, `m`, `w` as in did_pilot_rho(); `qs`: the
# covariate index of the control and of the treated group. The group
# intercepts are solved at the pilot's own rates. Returns the rho and
# whether the covariates alone already imply at least the dependence the
# pilot shows (then rho = 0 and the design has more dependence than the
# pilot).
did_pilot_rho_cov <- function(lf, p, m, w, qs) {
  qs <- lapply(qs, did_q_grid)
  use <- which(w > 0)
  cs <- vapply(1:2, function(g) did_solve(lf, qs[[g]], p[g]), numeric(1))
  f <- function(rho) {
    sum(vapply(use, function(g)
      w[g] * (did_joint_rate_cov(lf, cs[g], qs[[g]], rho) - m[g]), numeric(1)))
  }
  if (f(0) >= 0) return(list(rho = 0, exceeds = TRUE))
  if (f(0.995) <= 0) return(list(rho = 0.995, exceeds = FALSE))
  list(rho = uniroot(f, c(0, 0.995), tol = 1e-9)$root, exceeds = FALSE)
}

#' A difference-in-differences design from pre-treatment pilot data
#'
#' Estimates, from data observed **before** the treatment date, the inputs
#' of [ape_dgp_did()] that such data can inform, and builds the design:
#'
#' * `baseline`: the control group's outcome rate;
#' * `gap`: the treated group's rate minus the control group's;
#' * `p_treated`: the share of units (rows, without `id`) in the treated
#'   group;
#' * `rho`, for a panel: the latent within-unit correlation of the route's
#'   Gaussian copula, by the method of moments (the value at which the
#'   model's probability of the outcome in two periods of one unit equals
#'   the observed share of such pairs); with covariates, the correlation
#'   given the covariates (section 'Covariates and `rho`').
#'
#' All pre-treatment periods are pooled, because the route holds rates flat
#' before the treatment date. Pre-treatment data cannot tell how the control
#' group will change from before to after treatment, so `trend` must be
#' supplied; when the pilot has several periods, the control group's
#' average change from one to the next is reported as a guide. Covariates
#' are not taken from the pilot: pass them as [pa_var()] objects.
#'
#' Pooling assumes that the rates are flat over the pilot's periods. When
#' they drift, the estimated `rho` is too low, because the pooled rate
#' overstates the product of two periods' rates: with rates that change by
#' ten points from one period to the next, a true `rho` of 0.40 came out as
#' 0.34. The error lowers planned power, so it is on the cautious side; if
#' the reported change between pilot periods is large, consider supplying
#' `rho` yourself.
#'
#' @section Covariates and `rho`:
#' In [ape_dgp_did()], `rho` is the dependence between a unit's outcomes
#' that is left given the covariates; covariates are constant over time and
#' add dependence of their own. The pairs of outcomes in a pilot show the
#' whole dependence. When covariates with a nonzero `signal` (or `gamma`)
#' are passed, the function therefore solves for the `rho` at which the
#' design, with those covariates, has the pilot's share of pairs with the
#' outcome in both periods. That value is lower than the raw correlation
#' of the pilot, which is reported next to it (`$pilot$rho_raw`). If the
#' covariates as specified imply more dependence than the pilot shows,
#' `rho` is set to 0 with a warning: the stated `signal` is then too high
#' for these data. (Versions up to 1.12.0 used the raw correlation together
#' with the covariates, which counted part of the dependence twice and
#' overstated power.)
#'
#' The estimates are used as if they were known. They are not: vary them
#' with [ape_robust()] (`baseline`, `gap`, `rho`, `p_treated`), in
#' particular for a small pilot.
#'
#' @param data A data frame of pre-treatment observations.
#' @param outcome,group Names of the binary outcome and of the indicator of
#'   the group that will be treated (0/1 or logical).
#' @param period Optional name of the period variable; used to report the
#'   control group's change between pre-treatment periods. Periods are
#'   ordered by value when they are numbers or dates, by level when the
#'   variable is a factor, and alphabetically when it is other text.
#' @param id Optional name of the unit identifier. With it the design is a
#'   panel and `rho` is estimated (units need two or more periods);
#'   without it, or with `design = "cross-sections"`, the design is
#'   repeated cross-sections.
#' @param trend The control group's change in the outcome rate from before
#'   to after treatment, in probability points. Required.
#' @param ... Further arguments of [ape_dgp_did()], such as `model`,
#'   `n_pre`, `n_post`, `trends`, `analysis`, `covariates`. Naming
#'   `baseline`, `gap`, `p_treated`, `rho` or `design` here overrides the
#'   estimate. A `subgroup` is not supported.
#'
#' @return A `powerape_dgp` as from [ape_dgp_did()], with the estimates and
#'   the pilot's size in `$pilot`. Pin the effect with [set_att()].
#' @examples
#' ## a pilot panel: 300 units observed in two pre-treatment periods
#' set.seed(1)
#' u <- rnorm(300)
#' treated <- rbinom(300, 1, 0.5)
#' pilot <- data.frame(
#'   firm = rep(1:300, 2), year = rep(1:2, each = 300),
#'   treated = rep(treated, 2),
#'   y = rbinom(600, 1, pnorm(-0.7 + 0.3 * rep(treated, 2) + 0.8 * rep(u, 2))))
#' d <- suppressWarnings(
#'   ape_dgp_did_pilot(pilot, outcome = "y", group = "treated", id = "firm",
#'                     period = "year", trend = 0.05, model = "logit"))
#' d <- set_att(d, 0.08)
#' d
#' @export
ape_dgp_did_pilot <- function(data, outcome, group, id = NULL, period = NULL,
                              trend, ...) {
  stopifnot(is.data.frame(data), is.character(outcome), length(outcome) == 1L,
            is.character(group), length(group) == 1L)
  args <- list(...)
  if (!is.null(args$subgroup))
    stop(paste("ape_dgp_did_pilot() does not estimate subgroup designs; build",
               "those with ape_dgp_did()."), call. = FALSE)
  need <- c(outcome, group, id, period)
  miss <- setdiff(need, names(data))
  if (length(miss))
    stop(sprintf("Not in `data`: %s.", paste(miss, collapse = ", ")), call. = FALSE)
  if (anyNA(data[need]))
    stop("`data` has missing values in the variables used; remove those rows first.",
         call. = FALSE)
  y <- as.numeric(data[[outcome]])
  g <- as.numeric(data[[group]])
  if (!all(y %in% c(0, 1)) || !all(g %in% c(0, 1)))
    stop("`outcome` and `group` must be binary (0/1 or logical).", call. = FALSE)
  if (!any(g == 0) || !any(g == 1))
    stop("The pilot needs observations of both groups.", call. = FALSE)
  panel <- !is.null(id)
  p0 <- mean(y[g == 0])
  p1 <- mean(y[g == 1])
  if (p0 <= 0 || p0 >= 1 || p1 <= 0 || p1 >= 1)
    stop("A group has no events or no non-events in the pilot: its rate cannot be planned from these data.",
         call. = FALSE)
  est <- list(baseline = p0, gap = p1 - p0, design = if (panel) "panel" else "cross-sections")
  rho_hat <- NA_real_
  n_pairs <- 0
  n_units <- NA_integer_
  if (panel) {
    u <- data[[id]]
    gi <- tapply(g, u, function(v) v[1L])
    if (any(tapply(g, u, function(v) length(unique(v))) > 1L))
      stop("`group` must be constant within a unit.", call. = FALSE)
    n_units <- length(gi)
    est$p_treated <- mean(gi)
    Ti <- tapply(y, u, length)
    si <- tapply(y, u, sum)
    w <- vapply(0:1, function(k) sum((Ti * (Ti - 1))[gi == k]), numeric(1))
    n_pairs <- sum(w) / 2
    if (n_pairs > 0) {
      m <- vapply(0:1, function(k) {
        if (w[k + 1L] == 0) 0 else sum((si * (si - 1))[gi == k]) / w[k + 1L]
      }, numeric(1))
      rho_hat <- did_pilot_rho(c(p0, p1), m, w)
      est$rho <- rho_hat
    } else if (is.null(args$rho)) {
      stop(paste("Every unit has one pre-treatment period, so the dependence",
                 "between a unit's periods cannot be estimated: supply `rho`."),
           call. = FALSE)
    }
  } else {
    est$p_treated <- mean(g)
  }
  ## the control group's change between pre-treatment periods, as a guide
  ctrl_change <- NA_real_
  n_per <- NA_integer_
  if (!is.null(period)) {
    ## periods in their own order: numbers and dates by value, a factor by
    ## its levels, text that reads as numbers by those numbers, other text
    ## alphabetically
    per <- data[[period]]
    if (is.character(per)) {
      num <- suppressWarnings(as.numeric(per))
      per <- if (!anyNA(num)) num else factor(per)
    }
    r0 <- tapply(y[g == 0], per[g == 0], mean)
    r0 <- r0[!is.na(r0)]
    n_per <- length(unique(data[[period]]))
    if (length(r0) >= 2L) ctrl_change <- mean(diff(as.numeric(r0)))
  }
  if (missing(trend))
    stop(sprintf(paste0(
      "`trend` is required: pre-treatment data do not show how the control ",
      "group will change from before to after treatment.%s"),
      if (is.na(ctrl_change)) ""
      else sprintf(paste0(" (Between the pilot's pre-treatment periods the ",
                          "control group's rate changed by %+.3f on average.)"),
                   ctrl_change)), call. = FALSE)
  overridden <- intersect(names(est), names(args))
  est[overridden] <- NULL
  ## asked for repeated cross-sections although the pilot has units: the
  ## estimated dependence has no use there and is not passed on
  if (identical(args$design, "cross-sections")) est$rho <- NULL
  if (nrow(data) < 500L)
    warning(paste("The pilot has fewer than 500 rows: the estimated rates and",
                  "the within-unit correlation are imprecise. Vary them with",
                  "ape_robust() before relying on a sample size."), call. = FALSE)
  d <- do.call(ape_dgp_did, c(est, list(trend = trend), args))
  ## covariates add within-unit dependence of their own: the design's rho
  ## is the part left given them (see 'Covariates and rho' above). The
  ## calibration does not depend on rho, so the design is updated in place.
  rho_raw <- rho_hat
  if (!is.na(rho_hat) && !("rho" %in% overridden) && did_panel(d) &&
      d$k > 0L && isTRUE(d$s2 > 0)) {
    qd <- did_integration(d)
    rc <- did_pilot_rho_cov(link_funs(d$model), c(p0, p1), m, w,
                            list(qd$q0, qd$q1))
    if (rc$exceeds)
      warning(sprintf(paste(
        "The covariates as specified imply more within-unit dependence than",
        "the pilot shows (its raw latent correlation is %.3f): `rho` given",
        "the covariates is set to 0, and the design still has more dependence",
        "than the pilot. The stated `signal` looks too high for these data."),
        rho_raw), call. = FALSE)
    rho_hat <- rc$rho
    d$rho <- rho_hat
    d$spec$rho <- rho_hat
  }
  d$pilot <- list(n_rows = nrow(data), n_units = n_units, n_periods = n_per,
                  n_pairs = n_pairs, baseline = p0, gap = p1 - p0,
                  p_treated = if (panel) mean(gi) else mean(g), rho = rho_hat,
                  rho_raw = rho_raw,
                  control_change = ctrl_change, overridden = overridden)
  d
}

# One line for the design print, and one clause for the methods paragraph.
did_pilot_line <- function(x) {
  pl <- x$pilot
  est <- setdiff(c("baseline", "gap", "treated share", if (!is.na(pl$rho)) "rho"),
                 sub("p_treated", "treated share", pl$overridden))
  sprintf("  pilot: %s estimated from %d pre-treatment rows%s%s%s",
          paste(est, collapse = ", "), pl$n_rows,
          if (is.na(pl$n_units)) "" else sprintf(" (%d units)", pl$n_units),
          if (is.na(pl$control_change)) ""
          else sprintf("; control group's change between pilot periods %+.3f", pl$control_change),
          if (did_pilot_converted(pl))
            sprintf("; rho given the covariates %.3f (raw within-unit correlation %.3f)",
                    pl$rho, pl$rho_raw)
          else "")
}
# TRUE when the pilot's rho was solved given the design's covariates.
did_pilot_converted <- function(pl) {
  !is.null(pl$rho_raw) && !is.na(pl$rho_raw) && !("rho" %in% pl$overridden) &&
    !isTRUE(all.equal(pl$rho, pl$rho_raw))
}
did_pilot_clause <- function(d) {
  pl <- d$pilot
  if (is.null(pl)) return(NULL)
  est <- setdiff(c("baseline", "gap", "p_treated", if (!is.na(pl$rho)) "rho"), pl$overridden)
  lab <- c(baseline = "the control group's rate before treatment",
           gap = "the gap between the groups", p_treated = "the treated share",
           rho = if (did_pilot_converted(pl)) "the within-unit correlation given the covariates"
                 else "the within-unit correlation")[est]
  if (!length(lab)) return(NULL)
  sprintf("%s estimated from %d pre-treatment pilot observations",
          if (length(lab) > 1L) paste(paste(lab[-length(lab)], collapse = ", "), "and", lab[length(lab)])
          else lab,
          pl$n_rows)
}
