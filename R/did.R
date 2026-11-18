# Difference-in-differences designs --------------------------------------------
# A treated and a control group observed in n_pre
# periods before and n_post periods after one common treatment date (one
# period each by default). The estimand is the average
# treatment effect on the treated (ATT) in the post period, in probability
# points. The untreated marginal mean is
#   E[Y_t(0) | D, Z] = G(a + b D + l t + th D t + q),   q = Z'gamma,
# with th = 0 when the untreated trends are parallel on the index
# (Wooldridge, 2023) and th solved so that the untreated rates are parallel
# in probabilities otherwise. The treated post-period mean adds `de` to the
# index. The analysis is a pooled probit/logit of Y on [1, D, period
# dummies, W, Z], W = D x post, and the ATT is the average partial effect
# of W over the treated post-period rows (Wooldridge, 2023, eq. 2.34;
# Puhani, 2012). Rates are flat within the pre block and within the post
# block, so the world is the two-by-two world replicated over periods.

is_did <- function(dgp) identical(dgp$route, "did")
did_panel <- function(dgp) is_did(dgp) && identical(dgp$design, "panel")

# The analysis of a difference-in-differences design (section 17.1):
# "index" = pooled probit/logit with the ATT as an average partial effect,
# "linear" = linear probability model (OLS), ATT = coefficient of the
# treatment indicator. An analysis is valid in the simulated world exactly
# when it assumes the scale on which the world's untreated trends are
# parallel.
did_linear <- function(dgp) is_did(dgp) && identical(dgp$analysis, "linear")
did_valid <- function(dgp) identical(dgp$trends, dgp$analysis %||% "index")

# Routes whose standard errors are fixed by the design (`se` is ignored) and
# routes whose `n` counts units observed over several periods.
se_fixed <- function(dgp) {
  dgp$route %in% c("panel", "iv") || did_panel(dgp) || is_linear(dgp)
}
unit_design <- function(dgp) identical(dgp$route, "panel") || did_panel(dgp)

# The fitted model as named in print headers: the link, or the linear
# probability model for a design analysed by OLS (linear.R).
model_label <- function(dgp) {
  if (is_linear(dgp)) "linear probability (OLS)" else dgp$model
}

# One clause naming a world/analysis mismatch (curve print, error messages).
did_mismatch_text <- function(dgp) {
  sprintf("the simulated world has %s-parallel trends and the %s analysis assumes %s-parallel trends",
          dgp$trends, model_label(dgp), dgp$analysis %||% "index")
}

# A search needs an analysis that is valid in the simulated world. When the
# analysis assumes another parallel-trends scale than the world has, it
# centres on a different number than the true ATT, so "the n that reaches
# 80%" has no meaning.
did_need_valid <- function(dgp, what) {
  if (is_did(dgp) && !did_valid(dgp))
    stop(sprintf(paste(
      "%s needs an analysis that is valid in the simulated world. Here %s,",
      "so the analysis is biased for the effect on the treated:",
      "ape_power() and ape_curve() report its rejection rate and the mean",
      "of its estimates. Search with `trends` and `analysis` on the same",
      "scale (both \"index\" or both \"linear\")."),
      what, did_mismatch_text(dgp)), call. = FALSE)
  invisible(TRUE)
}

# Group-conditional draws of the nuisance index q = Z'gamma (deterministic).
did_integration <- function(dgp) {
  if (dgp$k == 0L) return(list(q0 = 0, q1 = 0))
  rd <- reg_draw_parametric(dgp$focal, NULL, dgp$covariates, dgp$R,
                            dgp$n_int, dgp$seed_int)
  q <- drop(rd$Z %*% dgp$gamma)
  list(q0 = q[rd$xd == 0], q1 = q[rd$xd == 1])
}

# Solve c in mean(G(offset + c + q)) = rate (closed form without covariates).
did_solve <- function(lf, q, rate, offset = 0) {
  if (length(q) == 1L) return(lf$Ginv(rate) - offset - q)
  uniroot(function(c) mean(lf$G(offset + c + q)) - rate, c(-15, 15),
          extendInt = "upX", tol = 1e-10)$root
}

#' Specify a difference-in-differences design with a binary outcome
#'
#' Describes a difference-in-differences study with one treatment date: a
#' treated and a control group observed before and after the treated group
#' receives the intervention, in one period each by default or in `n_pre`
#' periods before and `n_post` periods after. The target is the **average
#' treatment effect on the treated (ATT)** after treatment, in probability
#' points, pinned with [set_att()]. The outcome is regressed on the group
#' indicator, period indicators, the treatment indicator (treated group
#' after the treatment date), and any covariates.
#' Under `analysis = "index"` (default) that regression is a pooled probit
#' or logit and the ATT is the average partial effect of the treatment
#' indicator over the treated post-period observations (Wooldridge, 2023),
#' with a delta-method standard error. Under `analysis = "linear"` it is a
#' linear probability model fitted by ordinary least squares and the ATT is
#' the coefficient of the treatment indicator, the standard
#' difference-in-differences estimate.
#'
#' The world is described on the probability scale by three rates:
#' `baseline` (control group, before), `gap` (treated minus control,
#' before), and `trend` (the control group's change from before to after).
#' With several periods the rates are the same in every period before and
#' in every period after the treatment date, and so is the effect.
#' What the treated group's rate would have been after the intervention
#' date *without* treatment then depends on the scale on which the two
#' groups' untreated trends are parallel, set by `trends`:
#'
#' * `"index"` (default): parallel on the probit/logit index, the
#'   assumption under which the probit/logit analysis identifies the ATT
#'   (Wooldridge, 2023; Puhani, 2012).
#' * `"linear"`: parallel in probabilities, the assumption under which the
#'   linear analysis identifies the ATT. With covariates the rates are
#'   parallel on average over each group's covariates.
#'
#' Each analysis is unbiased only in the world whose scale it assumes.
#' When `trends` and `analysis` name the same scale, [ape_power()] reports
#' power. When they differ, the analysis centres on another number than
#' the ATT: [ape_power()] and [ape_curve()] then report a **rejection rate**
#' together with the mean of the estimates, and [ape_n()] and [ape_mde()]
#' refuse the design. With two periods and no covariates the data cannot
#' tell the two scales apart (Wooldridge, 2023), and the same four cell
#' rates imply different effects under each: the printed design shows both
#' readings. The linear analysis always centres on the cross difference of
#' the four rates.
#'
#' @section Data structure:
#' `design = "panel"` follows the same units over all periods (a balanced
#' panel); `n` in the power functions counts units, and standard errors
#' are clustered by unit. Outcomes of one unit are dependent through a
#' Gaussian copula with latent correlation `rho` that leaves the marginal
#' rates exact under either link; under probit this is the random-effects
#' probit with `rho = Var(c) / (Var(c) + 1)`. `design = "cross-sections"`
#' draws independent observations, each after the treatment date with
#' probability `p_post` and then in one of that block's periods with equal
#' probability; `n` counts all observations, and `se` in [ape_power()]
#' chooses model-based or robust standard errors for the probit/logit
#' analysis.
#'
#' @section More than two periods:
#' The analysis includes an indicator for every period but the first and
#' one treatment indicator, so it estimates one effect, not one per
#' period. In a panel the dependence between two outcomes of a unit is the
#' same however far apart the periods are (one shared unit component).
#' That is the favourable case for adding periods: averaging over periods
#' then removes period-specific noise only, and the gain is larger than
#' under dependence that fades with the distance between periods, which
#' this version does not offer. As an example, with rates of .20 and .30
#' in the control group, .40 and .60 in the treated group and `rho = 0.5`,
#' two periods on each side instead of one cut the variance of the linear
#' panel estimate to 51%; if the latent correlation halved with every
#' period of distance it would be 71%. Treat the precision bought by extra
#' periods as an upper bound.
#'
#' @section Comparing two subgroups:
#' With `subgroup`, both subgroups are treated in the treated group and
#' the target is the **difference between the two subgroups' effects on
#' the treated**, in probability points (the estimand of a triple
#' difference that compares, say, the effect on women with the effect on
#' men; Caron, 2025). `baseline`, `gap` and `trend` may then hold two
#' numbers, for the reference subgroup (`subgroup = 0`) and for
#' `subgroup = 1`, and [set_att()] takes the two effects. The probit/logit
#' analysis is fully interacted with the subgroup and the estimate is the
#' difference of the two subgroups' average partial effects of treatment;
#' the linear analysis estimates the coefficient of the subgroup-by-
#' treatment product. Trends are parallel within each subgroup, on the
#' scale `trends` names. Three things differ from a single effect. The
#' same eight rates imply different subgroup differences on the two
#' scales, and the printed design shows both. A difference of zero in
#' probability points is not a zero interaction coefficient when the
#' subgroups have different rates. And with equal subgroups a difference
#' as large as the effect itself needs about four times the units the
#' effect needs. [ape_mde()] searches the difference, holding the
#' reference subgroup's effect at its pinned value.
#'
#' @section The linear analysis:
#' A linear probability model is heteroskedastic by construction, so its
#' standard errors are always robust and the `se` argument of the power
#' functions is ignored: HC1 in cross-sections (as Stata's
#' `regress, vce(robust)`) and clustered by unit in panels (with the
#' small-sample factor of Stata's `regress, vce(cluster)`). Ordinary least
#' squares has no separation, so no simulated study is discarded for an
#' all-zero or all-one cell; the sparse-cell warning still applies. In a
#' balanced panel a time-constant covariate changes neither the estimate
#' nor its clustered standard error beyond the degrees of freedom (its
#' contribution cancels in each unit's change), so `signal` buys the linear
#' panel analysis no precision; in cross-sections it does.
#'
#' @section Interval coverage in small samples:
#' Under the index analysis the effect on the treated is a nonlinear
#' function of the estimated rates, and its delta-method interval (the one
#' `margins` reports in Stata for the same model) covers slightly less than
#' its nominal level when the group-by-period cells are small. In the
#' design without covariates the coverage can be computed exactly. For a
#' 95% interval at a true effect of zero, with rates between 25% and 44%,
#' it is 94.0% with 125 observations per cell, 94.5% with 250, 94.7% with
#' 500 and 94.9% with 1,000 under logit (probit is up to 0.2 points closer
#' to 95%): the shortfall halves when the cells double. A detection claim
#' is therefore made for a true zero effect in 6.0% of studies at 125 per
#' cell and in 5.3% at 500, not in 5%. With covariates and in panels the
#' simulated coverage at about 500 per cell is similar (94.4% to 94.7%).
#' The simulated studies use this interval, so the power reported is the
#' power of the analysis a study will run and needs no correction; with
#' small cells, part of it is this excess. The interval of the linear
#' analysis is close to nominal at the same sizes (94.9% to 95.0% for the
#' same rates).
#'
#' @section Limits of this version:
#' One treated group, one treatment date, treatment assigned at the unit
#' level, balanced panels. Staggered adoption, effects that differ by
#' period, pre-trend tests, rates that drift within the before or after
#' block, and treatment assigned to whole clusters are not available.
#' All power functions accept these designs, including [ape_robust()]
#' (which can also vary the scale on which trends are parallel),
#' [power_at()], [ape_alpha()] and [power_statement()]. Functions that need
#' power stop for a design whose analysis assumes another scale than its
#' world, because its result is a rejection rate. [ape_dgp_did_pilot()]
#' builds the design from data observed before the treatment date.
#'
#' @param model `"probit"` (default) or `"logit"`.
#' @param baseline Outcome rate of the control group before treatment, in
#'   (0, 1).
#' @param gap Treated-group rate minus control-group rate before treatment
#'   (probability points; default 0).
#' @param trend Change in the control-group rate from before to after the
#'   treatment date (probability points; default 0).
#' @param n_pre,n_post Number of periods before and after the treatment
#'   date (default one each).
#' @param subgroup Optional binary [pa_var()] splitting both groups into
#'   two subgroups, e.g. `pa_var("female", "binary", p = 0.4)`; see the
#'   section 'Comparing two subgroups'.
#' @param p_treated Share of units (or observations) in the treated group.
#' @param design `"panel"` (default) or `"cross-sections"`.
#' @param rho Panel designs: latent correlation between any two outcomes
#'   of a unit, in `[0, 1)`.
#' @param p_post Cross-section designs: share of observations drawn after
#'   the treatment date (default `n_post / (n_pre + n_post)`).
#' @param trends Scale on which the untreated trends are parallel in the
#'   simulated world: `"index"` (default) or `"linear"`.
#' @param analysis The analysis the simulated studies receive: `"index"`
#'   (default; pooled probit/logit, ATT as an average partial effect) or
#'   `"linear"` (linear probability model, ATT as the coefficient of the
#'   treatment indicator, robust standard errors).
#' @param covariates List of [pa_var()] objects (may be empty);
#'   time-constant.
#' @param correlation `NULL`, a scalar, or a correlation matrix over
#'   `(group, covariates)`, or over `(group, subgroup, covariates)` when
#'   there is a subgroup -- Gaussian-copula dependence, as in [ape_dgp()].
#' @param signal,gamma Nuisance-signal strength or explicit nuisance
#'   coefficients, as in [ape_dgp()].
#' @param n_int,seed_int Size and seed of the deterministic integration
#'   draw used for calibration when there are covariates.
#'
#' @return An object of class `powerape_dgp` (route `"did"`). Pin the
#'   effect with [set_att()], then use [ape_power()], [ape_curve()],
#'   [ape_n()], or [ape_mde()].
#' @references
#' Caron, L. (2025). Triple difference designs with heterogeneous
#' treatment effects. arXiv:2502.19620.
#'
#' Puhani, P. A. (2012). The treatment effect, the cross difference, and
#' the interaction term in nonlinear "difference-in-differences" models.
#' *Economics Letters, 115*(1), 85-87.
#'
#' Wooldridge, J. M. (2023). Simple approaches to nonlinear
#' difference-in-differences with panel data. *The Econometrics Journal,
#' 26*(3), C31-C66.
#' @examples
#' d <- ape_dgp_did(model = "logit", baseline = 0.20, gap = 0.20,
#'                  trend = 0.10, rho = 0.4)
#' d <- set_att(d, target = 0.07)
#' d
#'
#' ## the linear probability analysis in a world with linear-parallel trends
#' dl <- ape_dgp_did(baseline = 0.20, gap = 0.20, trend = 0.10,
#'                   trends = "linear", analysis = "linear")
#' set_att(dl, target = 0.10)
#' @export
ape_dgp_did <- function(model = c("probit", "logit"), baseline, gap = 0,
                        trend = 0, p_treated = 0.5,
                        design = c("panel", "cross-sections"), rho = 0,
                        p_post = NULL, trends = c("index", "linear"),
                        analysis = c("index", "linear"),
                        n_pre = 1, n_post = 1, subgroup = NULL,
                        covariates = list(), correlation = NULL,
                        signal = 0, gamma = NULL,
                        n_int = 4e5, seed_int = 20260812L) {
  model <- match.arg(model)
  design <- match.arg(design)
  trends <- match.arg(trends)
  analysis <- match.arg(analysis)
  ## designs with a subgroup are built in did-subgroup.R (section 17.3)
  if (!is.null(subgroup))
    return(did_subgroup_dgp(model, baseline, gap, trend, p_treated, design,
                            rho, p_post, trends, analysis, n_pre, n_post,
                            subgroup, covariates, correlation, signal, gamma,
                            n_int, seed_int))
  if (inherits(covariates, "pa_var")) covariates <- list(covariates)
  stopifnot(all(vapply(covariates, inherits, logical(1), "pa_var")))
  num1 <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  stopifnot(num1(baseline), baseline > 0, baseline < 1,
            num1(gap), num1(trend),
            num1(p_treated), p_treated > 0, p_treated < 1,
            num1(rho), rho >= 0, rho < 1,
            num1(n_pre), n_pre >= 1, n_pre == round(n_pre),
            num1(n_post), n_post >= 1, n_post == round(n_post),
            is.numeric(n_int), n_int >= 1e4)
  n_pre <- as.integer(n_pre)
  n_post <- as.integer(n_post)
  p_post_in <- p_post
  if (is.null(p_post)) p_post <- n_post / (n_pre + n_post)
  stopifnot(num1(p_post), p_post > 0, p_post < 1)
  eps <- 1e-4
  in01 <- function(x) x > eps && x < 1 - eps
  if (!in01(baseline + gap))
    stop(sprintf(paste("`baseline` + `gap` = %.4f is not a rate strictly",
                       "inside (0, 1): the treated group's pre-period rate",
                       "must be."), baseline + gap), call. = FALSE)
  if (!in01(baseline + trend))
    stop(sprintf(paste("`baseline` + `trend` = %.4f is not a rate strictly",
                       "inside (0, 1): the control group's post-period rate",
                       "must be."), baseline + trend), call. = FALSE)
  if (trends == "linear" && !in01(baseline + gap + trend))
    stop(sprintf(paste(
      "With trends = \"linear\" the treated group's untreated post-period",
      "rate is baseline + gap + trend = %.4f, outside (0, 1). Linear",
      "parallel trends cannot hold for these rates (a bounded outcome",
      "cannot follow the control group's change in probability points);",
      "use trends = \"index\" or other rates."), baseline + gap + trend),
      call. = FALSE)
  if (design == "cross-sections" && rho != 0)
    warning("`rho` is ignored for design = \"cross-sections\" (independent observations).",
            call. = FALSE)

  k <- length(covariates)
  group <- pa_var("treated group", "binary", p = p_treated)
  R <- build_corr(1L + k, correlation)
  lf <- link_funs(model)
  if (is.null(gamma)) {
    stopifnot(num1(signal), signal >= 0, signal < 1)
    if (k == 0L && signal > 0)
      stop("`signal` > 0 requires at least one covariate.", call. = FALSE)
  } else if (!is.numeric(gamma) || length(gamma) != k) {
    stop(sprintf("`gamma` must be a numeric vector of length %d.", k), call. = FALSE)
  }

  if (k == 0L) {
    gamma_out <- numeric(0)
    s2 <- 0
    q0 <- q1 <- 0
  } else {
    rd <- reg_draw_parametric(group, NULL, covariates, R, n_int, seed_int)
    if (is.null(gamma)) {
      if (signal == 0) {
        gamma_out <- rep(0, k)
      } else {
        g_unit <- 1 / vapply(covariates, pa_var_sd, numeric(1))
        v_unit <- var(drop(rd$Z %*% g_unit))
        gamma_out <- sqrt((signal / (1 - signal)) / v_unit) * g_unit
      }
    } else {
      gamma_out <- gamma
    }
    q <- drop(rd$Z %*% gamma_out)
    s2 <- var(q)
    q0 <- q[rd$xd == 0]
    q1 <- q[rd$xd == 1]
  }

  ## calibration on the group-conditional covariate distributions
  a <- did_solve(lf, q0, baseline)
  b <- did_solve(lf, q1, baseline + gap, a)
  l <- did_solve(lf, q0, baseline + trend, a)
  if (trends == "index") {
    th <- 0
    cf <- mean(lf$G(a + b + l + q1))
  } else {
    cf <- baseline + gap + trend
    th <- did_solve(lf, q1, cf, a + b + l)
  }

  structure(list(
    model = model, link = model, G = lf$G, g = lf$g,
    focal = group, moderator = NULL, covariates = covariates, k = k,
    R = R, R_chol = chol(R),
    baseline = baseline,
    signal = if (is.null(gamma)) signal else NA_real_,
    signal_implied = s2 / (1 + s2), s2 = s2,
    gamma = gamma_out, beta0 = a,
    beta_focal = NULL, beta_mod = NULL, beta_int = NULL,
    estimand = NULL, target_est = NULL, p1 = NULL,
    main_focal = NULL, main_moderator = NULL,
    route = "did", builder = "did",
    design = design, trends = trends, analysis = analysis,
    gap = gap, trend = trend,
    p_treated = p_treated, p_post = p_post,
    rho = if (design == "panel") rho else 0,
    did = list(a = a, b = b, l = l, th = th, de = NULL),
    rates = list(c00 = baseline, c10 = baseline + gap,
                 c01 = baseline + trend, cf = cf, c11 = NULL),
    n_pre = n_pre, n_post = n_post,
    n_periods = n_pre + n_post, unbalanced = FALSE,
    focal_ref = 0, mod_ref = NULL,
    n_int = n_int, seed_int = seed_int,
    spec = list(model = model, baseline = baseline, gap = gap, trend = trend,
                p_treated = p_treated, design = design, rho = rho,
                p_post = p_post_in, trends = trends, analysis = analysis,
                n_pre = n_pre, n_post = n_post, covariates = covariates,
                correlation = correlation, signal = signal, gamma = gamma,
                n_int = n_int, seed_int = seed_int)
  ), class = "powerape_dgp")
}

#' Pin the true effect on the treated of a difference-in-differences design
#'
#' Solves for the treatment coefficient such that the design in `dgp` has
#' an average treatment effect on the treated of exactly `target`
#' probability points in the post period: the treated group's post-period
#' rate with treatment minus the rate it would have had without, the
#' latter fixed by `baseline`, `gap`, `trend`, and the `trends` scale of
#' [ape_dgp_did()].
#'
#' For a design with a `subgroup`, `target` holds the two subgroups'
#' effects, reference subgroup first, and the pinned estimand is their
#' difference (second minus first).
#'
#' @param dgp A `powerape_dgp` from [ape_dgp_did()].
#' @param target Assumed true ATT (the planning value), e.g. `0.05` for
#'   five probability points; for a design with a subgroup, the two
#'   subgroups' ATTs, e.g. `c(0.10, 0.04)`.
#'
#' @return The `dgp` with the treatment coefficient, `target_est`, and the
#'   treated group's post-period rate filled in.
#' @examples
#' d <- ape_dgp_did(baseline = 0.30, gap = 0.10, trend = 0.05)
#' d <- set_att(d, target = 0.08)
#' true_att(d)
#' @export
set_att <- function(dgp, target) {
  stopifnot(inherits(dgp, "powerape_dgp"))
  if (did_has_subgroup(dgp)) return(did_subgroup_set_att(dgp, target))
  stopifnot(inherits(dgp, "powerape_dgp"),
            is.numeric(target), length(target) == 1L, is.finite(target))
  if (!is_did(dgp))
    stop(paste("set_att() pins the effect on the treated of a",
               "difference-in-differences design; build one with",
               "ape_dgp_did(). For other designs use set_ape() or set_aie()."),
         call. = FALSE)
  cf <- dgp$rates$cf
  eps <- 1e-4
  if (target <= -cf + eps || target >= (1 - cf) - eps)
    stop(sprintf(paste(
      "Target ATT %.4f is infeasible for this design: without treatment",
      "the treated group's post-period rate is %.3f, so the attainable",
      "range is (%.3f, %.3f). Adjust the target, the rates, or `trends`."),
      target, cf, -cf, 1 - cf), call. = FALSE)
  lf <- link_funs(dgp$model)
  cc <- dgp$did
  q1 <- did_integration(dgp)$q1
  dgp$did$de <- did_solve(lf, q1, cf + target, cc$a + cc$b + cc$l + cc$th)
  dgp$rates$c11 <- cf + target
  dgp$p1 <- cf + target
  dgp$estimand <- "att"
  dgp$target_est <- target
  dgp
}

#' Recompute the true effect on the treated implied by a pinned design
#'
#' @param dgp A `powerape_dgp` from [ape_dgp_did()] after [set_att()].
#' @return The true ATT (numeric scalar); equals the [set_att()] target up
#'   to root-finding tolerance. For a design with a subgroup: the two
#'   subgroups' ATTs and their difference.
#' @export
true_att <- function(dgp) {
  if (did_has_subgroup(dgp)) return(did_subgroup_true_att(dgp))
  if (!is_did(dgp) || !identical(dgp$estimand, "att"))
    stop("This DGP has no pinned ATT - build it with ape_dgp_did() and call set_att().",
         call. = FALSE)
  cc <- dgp$did
  q1 <- did_integration(dgp)$q1
  off <- cc$a + cc$b + cc$l + cc$th
  mean(dgp$G(off + cc$de + q1) - dgp$G(off + q1))
}

# What the design's four observable cell rates would be read as under the
# OTHER parallel-trends scale: the cross difference of probabilities (what a
# linear difference-in-differences centres on), and the index-scale reading
# c11 - G(G^-1(c10) + G^-1(c01) - G^-1(c00)) (what the covariate-free
# probit/logit analysis centres on). Marginal rates throughout.
did_readings <- function(dgp) {
  if (did_has_subgroup(dgp)) return(did_subgroup_readings(dgp))
  r <- dgp$rates
  if (is.null(r$c11)) return(NULL)
  lf <- link_funs(dgp$model)
  list(linear = r$c11 - (r$c10 + r$c01 - r$c00),
       index = r$c11 - lf$G(lf$Ginv(r$c10) + lf$Ginv(r$c01) - lf$Ginv(r$c00)))
}

# True coefficients in design-matrix order [1, D, f_2..f_T, W, Z]: period
# dummies with the first (pre) period as reference, so 0 on the remaining
# pre-periods and `l` on every post-period. The untreated D x post term
# `th` and the treatment term `de` both load on W. With one period each
# this is [1, D, post, W, Z].
did_beta <- function(dgp) {
  if (did_has_subgroup(dgp)) return(did_subgroup_beta(dgp))
  if (is.null(dgp$did$de))
    stop("This DGP has no target effect yet - call set_att() first.", call. = FALSE)
  cc <- dgp$did
  c(cc$a, cc$b, rep(0, dgp$n_pre - 1L), rep(cc$l, dgp$n_post),
    cc$th + cc$de, dgp$gamma)
}

# One simulated study. Panel: n units x T rows (unit-major, periods in
# order), outcome rates conditional on the unit draw of a Gaussian copula
# so that the marginal rates stay exact under either link. Cross-sections:
# n independent rows, post-block membership Bernoulli(p_post) and, with
# more than two periods, a uniform period within the block (drawn only
# then, so two-period studies consume the same random numbers as before).
# Returns the group and post-block indicators next to the design matrix:
# they define the identifying cells, and with several periods `post` is
# not a design column.
draw_x_did <- function(dgp, n) {
  if (did_has_subgroup(dgp)) return(draw_x_did_subgroup(dgp, n))
  bt <- did_beta(dgp)
  U <- draw_std_mvn(n, dgp$R_chol)
  d <- transform_u(U[, 1], dgp$focal)
  Z <- if (dgp$k > 0L) {
    Zm <- U[, 1L + seq_len(dgp$k), drop = FALSE]
    for (j in seq_len(dgp$k)) Zm[, j] <- transform_u(Zm[, j], dgp$covariates[[j]])
    Zm
  } else NULL
  n_pre <- dgp$n_pre
  Tn <- dgp$n_periods
  if (identical(dgp$design, "panel")) {
    id <- rep(seq_len(n), each = Tn)
    D <- d[id]
    period <- rep(seq_len(Tn), times = n)
    post <- as.numeric(period > n_pre)
    if (!is.null(Z)) Z <- Z[id, , drop = FALSE]
  } else {
    id <- NULL
    D <- d
    post <- as.numeric(runif(n) < dgp$p_post)
    period <- if (Tn == 2L) {
      1L + post
    } else {
      ifelse(post == 1, n_pre + sample.int(dgp$n_post, n, replace = TRUE),
             sample.int(n_pre, n, replace = TRUE))
    }
  }
  Fm <- if (Tn == 2L) post else outer(period, 2:Tn, "==") * 1
  X <- cbind(1, D, Fm, D * post, Z)
  p <- dgp$G(drop(X %*% bt))
  pr <- if (!is.null(id) && dgp$rho > 0) {
    a_i <- rnorm(n)
    pnorm((qnorm(p) - sqrt(dgp$rho) * a_i[id]) / sqrt(1 - dgp$rho))
  } else p
  list(X = X, focal_col = Tn + 2L, id = id, pr = pr, start = bt,
       d = D, post = post)
}

# ATT estimator and Jacobian: the average partial effect of the treatment
# indicator (design-matrix column `w_col`) over the rows where it is 1.
#   att = mean_{W=1} [ G(x'b) - G(x'b - b_W) ]
#   jac = mean_{W=1} [ g(a1) x1 - g(a0) x0 ]
# Returns NULL when the study has no treated post-period row.
att_est <- function(beta, X, link, w_col = 4L) {
  rows <- X[, w_col] == 1
  if (!any(rows)) return(NULL)
  X1 <- X[rows, , drop = FALSE]
  X0 <- X1
  X0[, w_col] <- 0
  G <- if (link == "probit") pnorm else plogis
  gd <- if (link == "probit") dnorm else dlogis
  a1 <- drop(X1 %*% beta)
  a0 <- drop(X0 %*% beta)
  list(ape = mean(G(a1) - G(a0)),
       jac = colMeans(gd(a1) * X1 - gd(a0) * X0))
}


# print method for difference-in-differences designs (print.powerape_dgp
# dispatches here).
print_did_dgp <- function(x) {
  if (did_has_subgroup(x)) return(print_did_subgroup(x))
  r <- x$rates
  panel <- identical(x$design, "panel")
  lin <- did_linear(x)
  two <- x$n_periods == 2L
  cat(sprintf("powerape DGP -- %s, difference-in-differences (%s)\n", x$model,
              if (two) "two periods"
              else sprintf("%d periods: %d before, %d after", x$n_periods,
                           x$n_pre, x$n_post)))
  cat(sprintf("  data: %s\n",
              if (panel) {
                sprintf("panel, %s per unit; rho = %.2f (latent within-unit correlation); unit-clustered SEs",
                        if (two) "both periods" else sprintf("all %d periods", x$n_periods),
                        x$rho)
              } else {
                sprintf("repeated cross-sections, %.0f%% of observations in the post %s",
                        100 * x$p_post, if (two) "period" else "periods")
              }))
  cat(sprintf("  analysis: %s\n",
              if (lin) {
                sprintf("linear probability model (OLS); ATT = coefficient of the treatment indicator;%s",
                        sub("^,", "", linear_se_label(x)))
              } else {
                sprintf("pooled %s; ATT = average partial effect of the treatment indicator over treated post-period rows",
                        x$model)
              }))
  cat(sprintf("  treated group: %.0f%% of %s\n", 100 * x$p_treated,
              if (panel) "units" else "observations"))
  if (x$k > 0L) {
    cat("  covariates:",
        paste(vapply(x$covariates, pa_var_label, character(1)), collapse = ", "),
        "\n")
    sig <- if (identical(x$link, "logit") && isTRUE(x$s2 > 0)) {
      sprintf(paste("nuisance signal = %.3f (latent index at unit error variance;",
                    "McKelvey-Zavoina R2 on the logistic scale = %.3f)"),
              x$signal_implied, mz_r2_logit(x$s2))
    } else {
      sprintf("nuisance signal (latent pseudo-R2) = %.3f", x$signal_implied)
    }
    cat(sprintf("  %s\n", sig))
  } else {
    cat("  covariates: none\n")
  }
  cat(sprintf("  untreated trends: parallel %s\n",
              if (identical(x$trends, "index")) {
                if (lin) sprintf("on the %s index; the linear probability analysis assumes parallel trends in probabilities",
                                 x$model)
                else sprintf("on the %s index (the scale the %s analysis assumes)",
                             x$model, x$model)
              } else {
                if (lin) "in probabilities (linear; the scale the linear probability analysis assumes)"
                else sprintf("in probabilities (linear); the %s analysis assumes the index scale",
                             x$model)
              }))
  cat(sprintf("  rates before -> after: control %.3f -> %.3f; treated %.3f -> %s\n",
              r$c00, r$c01, r$c10,
              if (is.null(r$c11)) sprintf("%.3f without treatment", r$cf)
              else sprintf("%.3f (%.3f without treatment)", r$c11, r$cf)))
  if (!is.null(x$pilot)) cat(did_pilot_line(x), "\n", sep = "")
  if (!is.null(x$target_est)) {
    cat(sprintf("  true ATT pinned at %+.4f (treatment coefficient = %.4f)\n",
                x$target_est, x$did$de))
    rd <- did_readings(x)
    if (identical(x$trends, "index")) {
      cat(sprintf("  the same four rates read under linear-parallel trends: %+.4f (their cross difference)\n",
                  rd$linear))
    } else {
      cat(sprintf("  the same four rates read on the %s index: %+.4f (covariate-free reading)\n",
                  x$model, rd$index))
    }
    if (!did_valid(x))
      cat(sprintf("  note: the analysis assumes %s-parallel trends, so it centres on that reading, not on the ATT\n",
                  x$analysis))
  } else {
    cat("  target effect not set yet - call set_att().\n")
  }
  invisible(x)
}
