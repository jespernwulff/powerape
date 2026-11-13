# Fixed datasets for the regression tests against Stata
# (test-stata-constants.R). They are built with base R's generators from
# fixed seeds, so they are the same in every session, and the caller's RNG
# state is left alone. The same function produced the files Stata read when
# the stored constants were computed (StataNow 19.5, revision 14 Oct 2026,
# on 2026-11-11).
stata_fixture <- function(name) {
  ws <- getFromNamespace("with_seed", "powerape")
  switch(name,
    ## a clustered cross-section for the linear probability analysis
    lin = ws(20261112L, {
      n <- 600L
      id <- rep(1:150, each = 4L)
      a <- rnorm(150)[id]
      d <- rbinom(n, 1, 0.4)
      z <- rnorm(n) + 0.5 * d
      m <- rbinom(n, 1, 0.5)
      y <- rbinom(n, 1, plogis(-0.8 + 0.5 * d + 0.4 * z + 0.3 * m + 0.4 * d * m + 0.6 * a))
      data.frame(y, d, m, z, id)
    }),
    ## difference-in-differences panel: 400 units, one period before and after
    did = ws(20261113L, {
      G <- 400L
      id <- rep(seq_len(G), each = 2L)
      period <- rep(1:2, G)
      d <- rep(rbinom(G, 1, 0.5), each = 2L)
      z <- rep(rnorm(G), each = 2L)
      a <- rep(rnorm(G), each = 2L)
      w <- d * (period == 2L)
      y <- rbinom(2L * G, 1, plogis(-1 + 0.5 * d + 0.4 * (period == 2L) + 0.6 * w +
                                      0.5 * z + 0.8 * a))
      data.frame(y, id, period, d, w, z)
    }),
    ## repeated cross-sections over four periods, treatment from period 3
    didx = ws(20261114L, {
      n <- 1200L
      d <- rbinom(n, 1, 0.5)
      period <- 1L + as.integer(floor(4 * runif(n)))
      post <- as.integer(period >= 3L)
      w <- d * post
      y <- rbinom(n, 1, plogis(-0.9 + 0.6 * d + 0.1 * (period == 2L) + 0.4 * post +
                                 0.15 * (period == 4L) + 0.5 * w))
      data.frame(y, d, period, w)
    }),
    ## difference-in-differences panel with a binary subgroup m
    sub = ws(20261115L, {
      G <- 600L
      id <- rep(seq_len(G), each = 2L)
      period <- rep(1:2, G)
      d <- rep(rbinom(G, 1, 0.5), each = 2L)
      m <- rep(rbinom(G, 1, 0.4), each = 2L)
      a <- rep(rnorm(G), each = 2L)
      post <- as.integer(period == 2L)
      w <- d * post
      y <- rbinom(2L * G, 1, plogis(-0.8 + 0.4 * d + 0.3 * post + 0.5 * w +
                                      m * (-0.4 + 0.3 * d - 0.1 * post + 0.4 * w) + 0.7 * a))
      data.frame(y, id, period, d, w, m)
    }),
    ## three-way interactions: binary f, m, w and continuous c1, c2
    tw = ws(20261116L, {
      n <- 1500L
      f <- rbinom(n, 1, 0.5); m <- rbinom(n, 1, 0.45); w <- rbinom(n, 1, 0.5)
      c1 <- rnorm(n, 1, 1.2); c2 <- rnorm(n, -0.5, 0.8); x <- rnorm(n)
      idx <- function(a, b) -0.7 + 0.35 * a + 0.2 * b + 0.2 * w + 0.15 * a * b -
        0.1 * a * w + 0.1 * b * w + 0.2 * a * b * w + 0.4 * x
      data.frame(y_bb = rbinom(n, 1, plogis(idx(f, m))),
                 y_cb = rbinom(n, 1, plogis(idx(c1, m))),
                 y_cc = rbinom(n, 1, plogis(idx(c1, c2))),
                 y_bbp = rbinom(n, 1, pnorm(idx(f, m))),
                 f, m, w, c1, c2, x)
    }),
    stop("unknown fixture: ", name))
}
