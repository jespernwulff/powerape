# Difference-in-differences designs with a binary subgroup ---------------------
# Both subgroups are treated in the treated group;
# the estimand is the difference between the two subgroups' effects on the
# treated in probability points (Caron, 2025: the DATT),
#   ATT(subgroup = 1) - ATT(subgroup = 0).
# Each subgroup m has its own untreated index a_m + b_m D + l_m post +
# th_m D post + q and its own treatment shift de_m; trends are parallel
# within each subgroup on the world's scale. The analysis is fully
# interacted with the subgroup: design matrix
#   [1, D, f_2..f_T, W, M, M D, M f_2..M f_T, M W, Z]
# with W in column T + 2 and M W in column 2T + 4. Design objects of these
# designs carry `subgroup`; everything in this file is reached only then.

did_has_subgroup <- function(dgp) is_did(dgp) && !is.null(dgp$subgroup)

# Last focal-side design column (group, period dummies, treatment indicator,
# and their subgroup counterparts): the range the Hauck-Donner backup scans.
did_last_focal_col <- function(dgp) {
  if (did_has_subgroup(dgp)) 2L * dgp$n_periods + 4L else dgp$n_periods + 2L
}

# Nuisance-index draws by group and subgroup: q[[d + 1]][[m + 1]].
did_cells_q <- function(dgp) {
  if (dgp$k == 0L) return(list(list(0, 0), list(0, 0)))
  rd <- reg_draw_parametric(dgp$focal, dgp$subgroup, dgp$covariates, dgp$R,
                            dgp$n_int, dgp$seed_int)
  q <- drop(rd$Z %*% dgp$gamma)
  lapply(0:1, function(d) lapply(0:1, function(m) q[rd$xd == d & rd$md == m]))
}

did_subgroup_dgp <- function(model, baseline, gap, trend, p_treated, design,
                             rho, p_post, trends, analysis, n_pre, n_post,
                             subgroup, covariates, correlation, signal, gamma,
                             n_int, seed_int) {
  if (!inherits(subgroup, "pa_var") || !identical(subgroup$type, "binary"))
    stop("`subgroup` must be a binary pa_var(), e.g. pa_var(\"female\", \"binary\", p = 0.4).",
         call. = FALSE)
  if (inherits(covariates, "pa_var")) covariates <- list(covariates)
  stopifnot(all(vapply(covariates, inherits, logical(1), "pa_var")))
  num1 <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  two <- function(x, nm) {
    if (!is.numeric(x) || !(length(x) %in% 1:2) || any(!is.finite(x)))
      stop(sprintf(paste("`%s` must be one number (the same in both subgroups)",
                         "or two: the reference subgroup (%s = 0), then %s = 1."),
                   nm, subgroup$name, subgroup$name), call. = FALSE)
    rep_len(x, 2L)
  }
  baseline <- two(baseline, "baseline")
  gap <- two(gap, "gap")
  trend <- two(trend, "trend")
  stopifnot(all(baseline > 0), all(baseline < 1),
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
  for (m in 1:2) {
    lab <- sprintf("%s = %d", subgroup$name, m - 1L)
    if (!in01(baseline[m] + gap[m]))
      stop(sprintf(paste("Subgroup %s: `baseline` + `gap` = %.4f is not a rate",
                         "strictly inside (0, 1)."), lab, baseline[m] + gap[m]),
           call. = FALSE)
    if (!in01(baseline[m] + trend[m]))
      stop(sprintf(paste("Subgroup %s: `baseline` + `trend` = %.4f is not a rate",
                         "strictly inside (0, 1)."), lab, baseline[m] + trend[m]),
           call. = FALSE)
    if (trends == "linear" && !in01(baseline[m] + gap[m] + trend[m]))
      stop(sprintf(paste(
        "Subgroup %s: with trends = \"linear\" the treated group's untreated",
        "post-period rate is baseline + gap + trend = %.4f, outside (0, 1).",
        "Linear parallel trends cannot hold for these rates; use",
        "trends = \"index\" or other rates."), lab,
        baseline[m] + gap[m] + trend[m]), call. = FALSE)
  }
  if (design == "cross-sections" && rho != 0)
    warning("`rho` is ignored for design = \"cross-sections\" (independent observations).",
            call. = FALSE)

  k <- length(covariates)
  group <- pa_var("treated group", "binary", p = p_treated)
  R <- build_corr(2L + k, correlation)
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
    qc <- list(list(0, 0), list(0, 0))
  } else {
    rd <- reg_draw_parametric(group, subgroup, covariates, R, n_int, seed_int)
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
    qc <- lapply(0:1, function(d) lapply(0:1, function(m) q[rd$xd == d & rd$md == m]))
    if (any(vapply(unlist(qc, recursive = FALSE), length, integer(1)) < 100L))
      stop("A group-by-subgroup cell is (nearly) empty in the integration draw; check the shares and the correlation.",
           call. = FALSE)
  }

  ## calibration within each subgroup, on its own covariate distributions
  a <- b <- l <- th <- cf <- numeric(2)
  for (m in 1:2) {
    q0 <- qc[[1]][[m]]
    q1 <- qc[[2]][[m]]
    a[m] <- did_solve(lf, q0, baseline[m])
    b[m] <- did_solve(lf, q1, baseline[m] + gap[m], a[m])
    l[m] <- did_solve(lf, q0, baseline[m] + trend[m], a[m])
    if (trends == "index") {
      th[m] <- 0
      cf[m] <- mean(lf$G(a[m] + b[m] + l[m] + q1))
    } else {
      cf[m] <- baseline[m] + gap[m] + trend[m]
      th[m] <- did_solve(lf, q1, cf[m], a[m] + b[m] + l[m])
    }
  }

  structure(list(
    model = model, link = model, G = lf$G, g = lf$g,
    focal = group, moderator = NULL, subgroup = subgroup,
    covariates = covariates, k = k,
    R = R, R_chol = chol(R),
    baseline = baseline,
    signal = if (is.null(gamma)) signal else NA_real_,
    signal_implied = s2 / (1 + s2), s2 = s2,
    gamma = gamma_out, beta0 = a[1],
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
    att_sub = NULL,
    n_pre = n_pre, n_post = n_post,
    n_periods = n_pre + n_post, unbalanced = FALSE,
    focal_ref = 0, mod_ref = NULL,
    n_int = n_int, seed_int = seed_int,
    spec = list(model = model, baseline = baseline, gap = gap, trend = trend,
                p_treated = p_treated, design = design, rho = rho,
                p_post = p_post_in, trends = trends, analysis = analysis,
                n_pre = n_pre, n_post = n_post, subgroup = subgroup,
                covariates = covariates, correlation = correlation,
                signal = signal, gamma = gamma,
                n_int = n_int, seed_int = seed_int)
  ), class = "powerape_dgp")
}

# set_att() for subgroup designs: `target` = the two subgroups' effects on
# the treated (reference subgroup first). The estimand is their difference.
did_subgroup_set_att <- function(dgp, target) {
  nm <- dgp$subgroup$name
  if (!is.numeric(target) || length(target) != 2L || any(!is.finite(target)))
    stop(sprintf(paste(
      "This design has a subgroup: give set_att() the two effects on the",
      "treated, c(effect for %s = 0, effect for %s = 1). The estimand is",
      "their difference."), nm, nm), call. = FALSE)
  lf <- link_funs(dgp$model)
  cc <- dgp$did
  qc <- did_cells_q(dgp)
  eps <- 1e-4
  de <- numeric(2)
  for (m in 1:2) {
    cf <- dgp$rates$cf[m]
    if (target[m] <= -cf + eps || target[m] >= (1 - cf) - eps)
      stop(sprintf(paste(
        "Target ATT %.4f is infeasible for subgroup %s = %d: without treatment",
        "its treated group's post-period rate is %.3f, so the attainable",
        "range is (%.3f, %.3f)."), target[m], nm, m - 1L, cf, -cf, 1 - cf),
        call. = FALSE)
    de[m] <- did_solve(lf, qc[[2]][[m]], cf + target[m],
                       cc$a[m] + cc$b[m] + cc$l[m] + cc$th[m])
  }
  dgp$did$de <- de
  dgp$rates$c11 <- dgp$rates$cf + target
  dgp$att_sub <- target
  dgp$estimand <- "datt"
  dgp$target_est <- target[2] - target[1]
  dgp
}

did_subgroup_true_att <- function(dgp) {
  if (!identical(dgp$estimand, "datt"))
    stop("This design has no pinned effects - call set_att() with the two subgroup effects.",
         call. = FALSE)
  cc <- dgp$did
  qc <- did_cells_q(dgp)
  e <- vapply(1:2, function(m) {
    off <- cc$a[m] + cc$b[m] + cc$l[m] + cc$th[m]
    q1 <- qc[[2]][[m]]
    mean(dgp$G(off + cc$de[m] + q1) - dgp$G(off + q1))
  }, numeric(1))
  c(reference = e[1], subgroup = e[2], difference = e[2] - e[1])
}

# The two readings of the eight block rates, as differences between the
# subgroups: the triple cross difference (linear) and the difference of the
# two index readings.
did_subgroup_readings <- function(dgp) {
  r <- dgp$rates
  if (is.null(r$c11)) return(NULL)
  lf <- link_funs(dgp$model)
  lin <- r$c11 - (r$c10 + r$c01 - r$c00)
  idx <- r$c11 - lf$G(lf$Ginv(r$c10) + lf$Ginv(r$c01) - lf$Ginv(r$c00))
  list(linear = lin[2] - lin[1], index = idx[2] - idx[1],
       linear_by = lin, index_by = idx)
}

# True coefficients in design-matrix order: the reference subgroup's block
# [1, D, f_2..f_T, W], the differences of the other subgroup's block from
# it, then the covariates.
did_subgroup_beta <- function(dgp) {
  if (is.null(dgp$did$de))
    stop("This DGP has no target effect yet - call set_att() first.", call. = FALSE)
  cc <- dgp$did
  blk <- function(m) c(cc$a[m], cc$b[m], rep(0, dgp$n_pre - 1L),
                       rep(cc$l[m], dgp$n_post), cc$th[m] + cc$de[m])
  c(blk(1L), blk(2L) - blk(1L), dgp$gamma)
}

# One simulated study (see draw_x_did): the copula columns are (group,
# subgroup, covariates); the design matrix is the two-way block and its
# products with the subgroup indicator.
draw_x_did_subgroup <- function(dgp, n) {
  bt <- did_subgroup_beta(dgp)
  U <- draw_std_mvn(n, dgp$R_chol)
  d <- transform_u(U[, 1], dgp$focal)
  mm <- transform_u(U[, 2], dgp$subgroup)
  Z <- if (dgp$k > 0L) {
    Zm <- U[, 2L + seq_len(dgp$k), drop = FALSE]
    for (j in seq_len(dgp$k)) Zm[, j] <- transform_u(Zm[, j], dgp$covariates[[j]])
    Zm
  } else NULL
  n_pre <- dgp$n_pre
  Tn <- dgp$n_periods
  if (identical(dgp$design, "panel")) {
    id <- rep(seq_len(n), each = Tn)
    D <- d[id]
    M <- mm[id]
    period <- rep(seq_len(Tn), times = n)
    post <- as.numeric(period > n_pre)
    if (!is.null(Z)) Z <- Z[id, , drop = FALSE]
  } else {
    id <- NULL
    D <- d
    M <- mm
    post <- as.numeric(runif(n) < dgp$p_post)
    period <- if (Tn == 2L) {
      1L + post
    } else {
      ifelse(post == 1, n_pre + sample.int(dgp$n_post, n, replace = TRUE),
             sample.int(n_pre, n, replace = TRUE))
    }
  }
  Fm <- if (Tn == 2L) post else outer(period, 2:Tn, "==") * 1
  base <- cbind(1, D, Fm, D * post)
  X <- cbind(base, M * base, Z)
  p <- dgp$G(drop(X %*% bt))
  pr <- if (!is.null(id) && dgp$rho > 0) {
    a_i <- rnorm(n)
    pnorm((qnorm(p) - sqrt(dgp$rho) * a_i[id]) / sqrt(1 - dgp$rho))
  } else p
  list(X = X, focal_col = 2L * Tn + 4L, w_col = Tn + 2L, id = id, pr = pr,
       start = bt, d = D, post = post, m = M)
}

# Difference between the two subgroups' effects on the treated and its
# Jacobian. Subgroup 1: over its treated post-period rows, with the
# treatment indicator AND its subgroup product switched off; subgroup 0:
# over its treated post-period rows, with the treatment indicator switched
# off. NULL when either subgroup has no treated post-period row.
datt_est <- function(beta, xx, link) {
  X <- xx$X
  w <- xx$w_col
  mw <- xx$focal_col
  r1 <- X[, mw] == 1
  r0 <- X[, w] == 1 & xx$m == 0
  if (!any(r1) || !any(r0)) return(NULL)
  G <- if (link == "probit") pnorm else plogis
  gd <- if (link == "probit") dnorm else dlogis
  part <- function(rows, off) {
    X1 <- X[rows, , drop = FALSE]
    X0 <- X1
    X0[, off] <- 0
    a1 <- drop(X1 %*% beta)
    a0 <- drop(X0 %*% beta)
    list(ape = mean(G(a1) - G(a0)), jac = colMeans(gd(a1) * X1 - gd(a0) * X0))
  }
  e1 <- part(r1, c(w, mw))
  e0 <- part(r0, w)
  list(ape = e1$ape - e0$ape, jac = e1$jac - e0$jac, by = c(e0$ape, e1$ape))
}

# The estimate of a difference-in-differences study under the index
# analysis: the ATT, or the subgroup difference when the draw has one.
did_est <- function(beta, xx, link) {
  if (is.null(xx$m)) att_est(beta, xx$X, link, xx$focal_col)
  else datt_est(beta, xx, link)
}

# Identifying cells: group x before/after, within each subgroup if any.
did_cell_check <- function(y, xx) {
  if (is.null(xx$m)) return(cell_check(y, xx$d, xx$post))
  s0 <- xx$m == 0
  if (!any(s0) || all(s0)) return(list(sep = FALSE, min_cell = NA_real_))
  c0 <- cell_check(y[s0], xx$d[s0], xx$post[s0])
  c1 <- cell_check(y[!s0], xx$d[!s0], xx$post[!s0])
  list(sep = c0$sep || c1$sep, min_cell = min(c0$min_cell, c1$min_cell))
}

# print method for subgroup designs (print_did_dgp dispatches here).
print_did_subgroup <- function(x) {
  r <- x$rates
  nm <- x$subgroup$name
  panel <- identical(x$design, "panel")
  lin <- is_linear(x)
  two <- x$n_periods == 2L
  cat(sprintf("powerape DGP -- %s, difference-in-differences with a subgroup (%s)\n",
              x$model,
              if (two) "two periods"
              else sprintf("%d periods: %d before, %d after", x$n_periods,
                           x$n_pre, x$n_post)))
  cat(sprintf("  data: %s\n",
              if (panel) {
                sprintf("panel, %s per unit; rho = %.2f (%s); unit-clustered SEs",
                        if (two) "both periods" else sprintf("all %d periods", x$n_periods),
                        x$rho, did_rho_label(x))
              } else {
                sprintf("repeated cross-sections, %.0f%% of observations in the post %s",
                        100 * x$p_post, if (two) "period" else "periods")
              }))
  cat(sprintf("  analysis: %s\n",
              if (lin) {
                sprintf("linear probability model (OLS), interacted with %s; estimate = coefficient of %s x treatment;%s",
                        nm, nm, sub("^,", "", linear_se_label(x)))
              } else {
                sprintf("pooled %s, fully interacted with %s; estimate = difference of the two subgroups' average partial effects of treatment",
                        x$model, nm)
              }))
  cat(sprintf("  treated group: %.0f%% of %s; subgroup %s = 1: %.0f%%\n",
              100 * x$p_treated, if (panel) "units" else "observations",
              nm, 100 * x$subgroup$p))
  if (x$k > 0L) {
    cat("  covariates:",
        paste(vapply(x$covariates, pa_var_label, character(1)), collapse = ", "),
        "\n")
    cat(sprintf("  nuisance signal (latent index at unit error variance) = %.3f\n",
                x$signal_implied))
  } else {
    cat("  covariates: none\n")
  }
  cat(sprintf("  untreated trends: parallel within each subgroup, %s\n",
              if (identical(x$trends, "index")) sprintf("on the %s index", x$model)
              else "in probabilities (linear)"))
  for (m in 1:2)
    cat(sprintf("  rates before -> after, %s = %d: control %.3f -> %.3f; treated %.3f -> %s\n",
                nm, m - 1L, r$c00[m], r$c01[m], r$c10[m],
                if (is.null(r$c11)) sprintf("%.3f without treatment", r$cf[m])
                else sprintf("%.3f (%.3f without treatment)", r$c11[m], r$cf[m])))
  if (!is.null(x$target_est)) {
    cat(sprintf("  effects on the treated: %+.4f (%s = 0), %+.4f (%s = 1); true difference (DATT) pinned at %+.4f\n",
                x$att_sub[1], nm, x$att_sub[2], nm, x$target_est))
    rd <- did_subgroup_readings(x)
    if (identical(x$trends, "index")) {
      cat(sprintf("  the same eight rates read under linear-parallel trends: difference %+.4f (their triple cross difference)\n",
                  rd$linear))
    } else {
      cat(sprintf("  the same eight rates read on the %s index: difference %+.4f (covariate-free reading)\n",
                  x$model, rd$index))
    }
    if (!did_valid(x))
      cat(sprintf("  note: the analysis assumes %s-parallel trends, so it centres on that reading, not on the true difference\n",
                  x$analysis))
  } else {
    cat(sprintf("  effects not set yet - call set_att(d, c(effect for %s = 0, effect for %s = 1)).\n",
                nm, nm))
  }
  invisible(x)
}
