# A difference-in-differences design from pre-treatment pilot data -------------
# DESIGN.md section 20. Data observed before the treatment date inform four
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
#'   the observed share of such pairs).
#'
#' All pre-treatment periods are pooled, because the route holds rates flat
#' before the treatment date. Pre-treatment data cannot tell how the control
#' group will change from before to after treatment, so `trend` must be
#' supplied; when the pilot has several periods, the control group's
#' average change from one to the next is reported as a guide. Covariates
#' are not taken from the pilot: pass them as [pa_var()] objects.
#'
#' The estimates are used as if they were known. They are not: vary them
#' with [ape_robust()] (`baseline`, `gap`, `rho`, `p_treated`), in
#' particular for a small pilot.
#'
#' @param data A data frame of pre-treatment observations.
#' @param outcome,group Names of the binary outcome and of the indicator of
#'   the group that will be treated (0/1 or logical).
#' @param id Optional name of the unit identifier. With it the design is a
#'   panel and `rho` is estimated (units need two or more periods);
#'   without it the design is repeated cross-sections.
#' @param period Optional name of the period variable; used to report the
#'   control group's change between pre-treatment periods.
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
    r0 <- tapply(y[g == 0], data[[period]][g == 0], mean)
    r0 <- r0[order(as.numeric(names(r0)), names(r0))]
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
  if (nrow(data) < 500L)
    warning(paste("The pilot has fewer than 500 rows: the estimated rates and",
                  "the within-unit correlation are imprecise. Vary them with",
                  "ape_robust() before relying on a sample size."), call. = FALSE)
  d <- do.call(ape_dgp_did, c(est, list(trend = trend), args))
  d$pilot <- list(n_rows = nrow(data), n_units = n_units, n_periods = n_per,
                  n_pairs = n_pairs, baseline = p0, gap = p1 - p0,
                  p_treated = if (panel) mean(gi) else mean(g), rho = rho_hat,
                  control_change = ctrl_change, overridden = overridden)
  d
}

# One line for the design print, and one clause for the methods paragraph.
did_pilot_line <- function(x) {
  pl <- x$pilot
  est <- setdiff(c("baseline", "gap", "treated share", if (!is.na(pl$rho)) "rho"),
                 sub("p_treated", "treated share", pl$overridden))
  sprintf("  pilot: %s estimated from %d pre-treatment rows%s%s",
          paste(est, collapse = ", "), pl$n_rows,
          if (is.na(pl$n_units)) "" else sprintf(" (%d units)", pl$n_units),
          if (is.na(pl$control_change)) ""
          else sprintf("; control group's change between pilot periods %+.3f", pl$control_change))
}
did_pilot_clause <- function(d) {
  pl <- d$pilot
  if (is.null(pl)) return(NULL)
  est <- setdiff(c("baseline", "gap", "p_treated", if (!is.na(pl$rho)) "rho"), pl$overridden)
  lab <- c(baseline = "the control group's rate before treatment",
           gap = "the gap between the groups", p_treated = "the treated share",
           rho = "the within-unit correlation")[est]
  if (!length(lab)) return(NULL)
  sprintf("%s estimated from %d pre-treatment pilot observations",
          if (length(lab) > 1L) paste(paste(lab[-length(lab)], collapse = ", "), "and", lab[length(lab)])
          else lab,
          pl$n_rows)
}
