# Estimation internals: ML fit, APE/AIE estimators, delta-method SE -----------

# Maximum-likelihood fit of the index model. `start` at the true coefficients
# speeds up IRLS; the MLE does not depend on the start (concave likelihood).
# The runaway guard works on standardized slopes so it is invariant to the
# regressors' scales.
fit_index_model <- function(X, y, link, start = NULL) {
  fam <- binomial(link = link)
  fit <- suppressWarnings(
    glm.fit(X, y, family = fam, start = start, control = list(maxit = 100L))
  )
  runaway <- FALSE
  if (ncol(X) > 1L && all(is.finite(fit$coefficients))) {
    sds <- apply(X[, -1L, drop = FALSE], 2, sd)
    runaway <- any(abs(fit$coefficients[-1L]) * pmax(sds, 1e-12) > 15)
  }
  ok <- isTRUE(fit$converged) && fit$rank == ncol(X) &&
    all(is.finite(fit$coefficients)) && !runaway
  list(fit = fit, ok = ok)
}

# Unscaled covariance of the MLE from the IRLS QR decomposition
# (same computation as summary.glm; dispersion is 1 for binomial).
vcov_from_glmfit <- function(fit) {
  p1 <- seq_len(fit$rank)
  chol2inv(fit$qr$qr[p1, p1, drop = FALSE])
}

# Heteroskedasticity-robust (HC0) sandwich for the probit/logit MLE:
# expected-information bread (as in R's sandwich package), score outer
# products as the meat. No small-sample adjustment.
vcov_sandwich <- function(fit, X, link) {
  G <- if (link == "probit") pnorm else plogis
  gd <- if (link == "probit") dnorm else dlogis
  eta <- drop(X %*% fit$coefficients)
  mu <- G(eta)
  w <- gd(eta) * (fit$y - mu) / pmax(mu * (1 - mu), 1e-12)
  B <- vcov_from_glmfit(fit)
  B %*% crossprod(X * w) %*% B
}

# Cluster-robust (CR) sandwich for the pooled probit/logit MLE: score-based
# meat aggregated by cluster, bread from the IRLS QR, with the G/(G-1)
# multiplier that matches Stata's vce(cluster) convention for ML estimators
# (verified against Stata margins in the validation battery).
vcov_cluster <- function(fit, X, id, link) {
  G <- if (link == "probit") pnorm else plogis
  gd <- if (link == "probit") dnorm else dlogis
  eta <- drop(X %*% fit$coefficients)
  mu <- G(eta)
  w <- gd(eta) * (fit$y - mu) / pmax(mu * (1 - mu), 1e-12)
  U <- rowsum(X * w, id)
  B <- vcov_from_glmfit(fit)
  Gn <- nrow(U)
  (Gn / (Gn - 1)) * (B %*% crossprod(U) %*% B)
}

# First and second derivatives of the density for the two links:
#   probit: g'(z) = -z phi(z),          g''(z) = (z^2 - 1) phi(z)
#   logit:  g'(z) = g(1 - 2G),          g''(z) = g(1 - 2G)^2 - 2 g^2
gprime_fun <- function(link) {
  if (link == "probit") function(z) -z * dnorm(z)
  else function(z) dlogis(z) * (1 - 2 * plogis(z))
}

gpp_fun <- function(link) {
  if (link == "probit") {
    function(z) (z^2 - 1) * dnorm(z)
  } else {
    function(z) {
      g <- dlogis(z)
      g * (1 - 2 * plogis(z))^2 - 2 * g^2
    }
  }
}

# APE of a binary focal regressor and its gradient w.r.t. beta:
#   ape = mean_i [ G(x_i'b | d=1) - G(x_i'b | d=0) ]
#   jac = mean_i [ g(a_i1) x_i1 - g(a_i0) x_i0 ]
ape_binary_est <- function(beta, X, focal_col, link) {
  X1 <- X; X1[, focal_col] <- 1
  X0 <- X; X0[, focal_col] <- 0
  G <- if (link == "probit") pnorm else plogis
  gd <- if (link == "probit") dnorm else dlogis
  a1 <- drop(X1 %*% beta)
  a0 <- drop(X0 %*% beta)
  ape <- mean(G(a1) - G(a0))
  jac <- colMeans(gd(a1) * X1 - gd(a0) * X0)
  list(ape = ape, jac = jac)
}

# APE (average derivative) of a continuous focal regressor:
#   ape = b_d * mean_i g(a_i)
#   jac_j = b_d * mean_i [ g'(a_i) x_ij ],  plus mean_i g(a_i) added at j = d
ape_cont_est <- function(beta, X, focal_col, link) {
  gd <- if (link == "probit") dnorm else dlogis
  gp <- gprime_fun(link)
  a <- drop(X %*% beta)
  ga <- gd(a)
  b_d <- beta[[focal_col]]
  ape <- b_d * mean(ga)
  jac <- b_d * colMeans(gp(a) * X)
  jac[focal_col] <- jac[focal_col] + mean(ga)
  list(ape = ape, jac = jac)
}

# Dispatch on the focal variable's type ("binary" or "normal").
ape_est <- function(beta, X, focal_col, link, focal_type) {
  if (focal_type == "binary") ape_binary_est(beta, X, focal_col, link)
  else ape_cont_est(beta, X, focal_col, link)
}

# AIE estimator and jacobian. Design-matrix layout is fixed by draw_x():
# column 2 = focal, 3 = moderator, 4 = focal*moderator interaction.
# Cases follow ginteff: double difference (binary x binary), difference of
# average derivatives (mixed), average cross-partial (continuous x
# continuous). All formulas are parameterization-free: they take the raw
# fitted coefficients.
aie_est <- function(beta, X, link, focal_type, mod_type) {
  G <- if (link == "probit") pnorm else plogis
  gd <- if (link == "probit") dnorm else dlogis
  gp <- gprime_fun(link)
  d_col <- 2L; m_col <- 3L; i_col <- 4L

  if (focal_type == "binary" && mod_type == "binary") {
    X11 <- X; X11[, d_col] <- 1; X11[, m_col] <- 1; X11[, i_col] <- 1
    X01 <- X; X01[, d_col] <- 0; X01[, m_col] <- 1; X01[, i_col] <- 0
    X10 <- X; X10[, d_col] <- 1; X10[, m_col] <- 0; X10[, i_col] <- 0
    X00 <- X; X00[, d_col] <- 0; X00[, m_col] <- 0; X00[, i_col] <- 0
    a11 <- drop(X11 %*% beta); a01 <- drop(X01 %*% beta)
    a10 <- drop(X10 %*% beta); a00 <- drop(X00 %*% beta)
    aie <- mean(G(a11) - G(a01) - G(a10) + G(a00))
    jac <- colMeans(gd(a11) * X11 - gd(a01) * X01 - gd(a10) * X10 + gd(a00) * X00)
    list(ape = aie, jac = jac)
  } else if (focal_type != "binary" && mod_type == "binary") {
    bd <- beta[[d_col]]; bi <- beta[[i_col]]
    X1 <- X; X1[, m_col] <- 1; X1[, i_col] <- X[, d_col]
    X0 <- X; X0[, m_col] <- 0; X0[, i_col] <- 0
    a1 <- drop(X1 %*% beta); a0 <- drop(X0 %*% beta)
    aie <- (bd + bi) * mean(gd(a1)) - bd * mean(gd(a0))
    jac <- (bd + bi) * colMeans(gp(a1) * X1) - bd * colMeans(gp(a0) * X0)
    jac[d_col] <- jac[d_col] + mean(gd(a1)) - mean(gd(a0))
    jac[i_col] <- jac[i_col] + mean(gd(a1))
    list(ape = aie, jac = jac)
  } else if (focal_type == "binary" && mod_type != "binary") {
    bm <- beta[[m_col]]; bi <- beta[[i_col]]
    X1 <- X; X1[, d_col] <- 1; X1[, i_col] <- X[, m_col]
    X0 <- X; X0[, d_col] <- 0; X0[, i_col] <- 0
    a1 <- drop(X1 %*% beta); a0 <- drop(X0 %*% beta)
    aie <- (bm + bi) * mean(gd(a1)) - bm * mean(gd(a0))
    jac <- (bm + bi) * colMeans(gp(a1) * X1) - bm * colMeans(gp(a0) * X0)
    jac[m_col] <- jac[m_col] + mean(gd(a1)) - mean(gd(a0))
    jac[i_col] <- jac[i_col] + mean(gd(a1))
    list(ape = aie, jac = jac)
  } else {
    gpp <- gpp_fun(link)
    bd <- beta[[d_col]]; bm <- beta[[m_col]]; bi <- beta[[i_col]]
    a <- drop(X %*% beta)
    dv <- X[, d_col]; mv <- X[, m_col]
    t1 <- bd + bi * mv
    t2 <- bm + bi * dv
    gpa <- gp(a); ga <- gd(a)
    aie <- mean(gpa * t1 * t2 + ga * bi)
    jac <- colMeans(gpp(a) * t1 * t2 * X) + bi * colMeans(gpa * X)
    jac[d_col] <- jac[d_col] + mean(gpa * t2)
    jac[m_col] <- jac[m_col] + mean(gpa * t1)
    jac[i_col] <- jac[i_col] + mean(gpa * (mv * t2 + dv * t1)) + mean(ga)
    list(ape = aie, jac = jac)
  }
}
