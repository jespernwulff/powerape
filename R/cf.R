# Control-function probit: estimator internals ---------------------------------
# DESIGN.md section 12. Replicates Stata 18.5's cfprobit/cfregress (StataNow):
# first-stage regression per endogenous variable (linear / probit / fprobit /
# poisson), control function = the first stage's score residual, second-stage
# probit by ML with the CF (and its interactions) as regressors, and standard
# errors from stacking the exact estimation moments of both stages -- an
# exactly identified method-of-moments system, so no bootstrap is needed.
# margins convention (cfregress_p.ado): inside margins the stored CFs are
# held FIXED, so margins dydx is the ASF-based APE over the empirical CF
# distribution. Our estimators mirror that convention.

# Probit score factor s(q, y) = phi(q)(y - Phi(q)) / {Phi(q)(1 - Phi(q))}
# and its derivative in q. This is simultaneously the probit/fprobit
# generalized residual (Stata's control function) and the main-equation
# "error function" whose products with the regressors are the GMM moments.
probit_score_bits <- function(q, y) {
  g <- dnorm(q)
  G <- pnorm(q)
  B <- pmax(G * (1 - G), 1e-12)
  A <- y - G
  s <- g * A / B
  ds <- (-q * g * A - g^2) / B - g * A * g * (1 - 2 * G) / B^2
  list(s = s, ds = ds)
}

# First-stage fit for one endogenous variable. Returns coefficients, the
# control function c, and dc/dq (q = the first-stage linear index).
cf_first <- function(y2, X1, type = c("linear", "probit", "fprobit", "poisson")) {
  type <- match.arg(type)
  if (type == "linear") {
    ft <- lm.fit(X1, y2)
    q <- drop(X1 %*% ft$coefficients)
    return(list(pi = ft$coefficients, q = q, c = y2 - q,
                dcdq = rep(-1, length(q)), ok = ft$rank == ncol(X1)))
  }
  if (type == "poisson") {
    ft <- suppressWarnings(glm.fit(X1, y2, family = poisson("log"),
                                   control = list(maxit = 100L)))
    q <- drop(X1 %*% ft$coefficients)
    ok <- isTRUE(ft$converged) && ft$rank == ncol(X1) &&
      all(is.finite(ft$coefficients))
    return(list(pi = ft$coefficients, q = q, c = y2 - exp(q),
                dcdq = -exp(q), ok = ok))
  }
  ## probit / fprobit: same fit (Bernoulli QMLE for fractional y) and the
  ## same generalized-residual formula
  ft <- suppressWarnings(glm.fit(X1, y2, family = binomial("probit"),
                                 control = list(maxit = 100L)))
  q <- drop(X1 %*% ft$coefficients)
  sb <- probit_score_bits(q, y2)
  ok <- isTRUE(ft$converged) && ft$rank == ncol(X1) &&
    all(is.finite(ft$coefficients))
  list(pi = ft$coefficients, q = q, c = sb$s, dcdq = sb$ds, ok = ok)
}

# Two-step control-function probit with stacked no-bootstrap vcov.
#
# Main-equation layout (fixed, mirroring the package's conventions):
#   no moderator:  W = [1, y2, Z..., c]                cols: focal 2, cf last
#   moderator:     W = [1, y2, m, y2*m, Z..., c, c*m]  cols: 2,3,4; cf, cf*m last
# X1 is the first-stage design [1, Z..., instruments...] (+ z*m columns if
# supplied). `id` switches the meat from HC (robust) to cluster-robust.
cf_fit <- function(y0, y2, m = NULL, Z = NULL, X1, first = "linear",
                   id = NULL, start = NULL) {
  n <- length(y0)
  fs <- cf_first(y2, X1, first)
  if (!fs$ok) return(list(ok = FALSE))
  has_m <- !is.null(m)
  W <- if (has_m) {
    cbind(1, y2, m, y2 * m, Z, fs$c, fs$c * m)
  } else {
    cbind(1, y2, Z, fs$c)
  }
  cf_col <- ncol(W) - if (has_m) 1L else 0L
  cfm_col <- if (has_m) ncol(W) else NULL

  ft <- fit_index_model(W, y0, "probit", start = start)
  if (!ft$ok) return(list(ok = FALSE))
  theta <- ft$fit$coefficients

  ## stacked moment system: main probit score x W, first-stage score x X1
  idx <- drop(W %*% theta)
  mb <- probit_score_bits(idx, y0)
  k_m <- ncol(W)
  k_1 <- ncol(X1)

  ## A blocks (Jacobian of the summed moments). Stata's convention,
  ## verified to ~1e-7 against cfprobit: observed score derivatives, and
  ## the generated regressors (c, c*m) are treated as FIXED instruments
  ## when differentiating in pi -- only the error function's dependence
  ## on pi (through the index) enters, via d idx / d c = rho + gamma_cm*m.
  A_mm <- crossprod(W * mb$ds, W)
  didx_dc <- theta[[cf_col]] + if (has_m) theta[[cfm_col]] * m else 0
  A_m1 <- crossprod((mb$ds * didx_dc) * W, fs$dcdq * X1)
  A_11 <- crossprod(X1 * fs$dcdq, X1)
  A <- rbind(cbind(A_mm, A_m1),
             cbind(matrix(0, k_1, k_m), A_11))

  ## B: outer product of stacked per-observation scores (robust), or the
  ## raw sum of cluster-aggregated scores -- cfprobit applies NO
  ## small-sample factor to either (gmm-style; verified against Stata)
  S <- cbind(mb$s * W, fs$c * X1)
  B <- if (!is.null(id)) crossprod(rowsum(S, id)) else crossprod(S)
  Ainv <- tryCatch(solve(A), error = function(e) NULL)
  if (is.null(Ainv)) return(list(ok = FALSE))
  V <- Ainv %*% B %*% t(Ainv)
  V_main <- V[seq_len(k_m), seq_len(k_m), drop = FALSE]

  list(ok = TRUE, theta = theta, W = W, V_main = V_main, V = V,
       cf = fs$c, pi = fs$pi, cf_col = cf_col, cfm_col = cfm_col,
       has_m = has_m)
}

# ASF-based APE of the endogenous focal: toggle/differentiate y2 with the
# control function held fixed (Stata's margins convention). Columns that
# are functions of y2 (y2 itself, y2*m) update; c and c*m do not.
cf_ape_est <- function(theta, W, focal_type, cf_col, cfm_col = NULL) {
  k <- length(theta)
  if (focal_type == "binary") {
    W1 <- W; W1[, 2L] <- 1
    W0 <- W; W0[, 2L] <- 0
    if (!is.null(cfm_col)) {   # y2*m column toggles with y2
      W1[, 4L] <- W[, 3L]
      W0[, 4L] <- 0
    }
    a1 <- drop(W1 %*% theta); a0 <- drop(W0 %*% theta)
    ape <- mean(pnorm(a1) - pnorm(a0))
    jac <- colMeans(dnorm(a1) * W1 - dnorm(a0) * W0)
  } else {
    a <- drop(W %*% theta)
    ## d idx / d y2 = alpha (+ gamma_y2m * m)
    slope <- theta[[2L]] + if (!is.null(cfm_col)) theta[[4L]] * W[, 3L] else 0
    ga <- dnorm(a)
    ape <- mean(ga * slope)
    gp <- -a * ga
    jac <- colMeans(gp * slope * W)
    jac[2L] <- jac[2L] + mean(ga)
    if (!is.null(cfm_col)) jac[4L] <- jac[4L] + mean(ga * W[, 3L])
  }
  list(ape = ape, jac = jac)
}

# ASF-based AIE (endogenous focal x exogenous binary moderator): double
# difference over the four (y2, m) cells with c fixed; the c*m column
# follows the m toggle, the y2*m column follows both.
cf_aie_est <- function(theta, W, focal_type, cf_col, cfm_col) {
  cf <- W[, cf_col]
  cell <- function(d, m) {
    Wc <- W
    Wc[, 2L] <- d
    Wc[, 3L] <- m
    Wc[, 4L] <- d * m
    Wc[, cfm_col] <- cf * m
    Wc
  }
  if (focal_type == "binary") {
    W11 <- cell(1, 1); W01 <- cell(0, 1); W10 <- cell(1, 0); W00 <- cell(0, 0)
    a11 <- drop(W11 %*% theta); a01 <- drop(W01 %*% theta)
    a10 <- drop(W10 %*% theta); a00 <- drop(W00 %*% theta)
    aie <- mean(pnorm(a11) - pnorm(a01) - pnorm(a10) + pnorm(a00))
    jac <- colMeans(dnorm(a11) * W11 - dnorm(a01) * W01 -
                      dnorm(a10) * W10 + dnorm(a00) * W00)
    return(list(ape = aie, jac = jac))
  }
  ## continuous endogenous focal: difference of average derivatives at
  ## m = 1 vs m = 0, y2 at observed values, c fixed
  d_obs <- W[, 2L]
  deriv_at <- function(mval) {
    Wc <- W
    Wc[, 3L] <- mval
    Wc[, 4L] <- d_obs * mval
    Wc[, cfm_col] <- cf * mval
    a <- drop(Wc %*% theta)
    slope <- theta[[2L]] + theta[[4L]] * mval
    ga <- dnorm(a)
    gp <- -a * ga
    jac <- colMeans(gp * slope * Wc)
    jac[2L] <- jac[2L] + mean(ga)
    jac[4L] <- jac[4L] + mean(ga) * mval
    list(v = mean(ga * slope), jac = jac)
  }
  d1 <- deriv_at(1); d0 <- deriv_at(0)
  list(ape = d1$v - d0$v, jac = d1$jac - d0$jac)
}
