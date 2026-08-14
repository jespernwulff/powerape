#' Declare a variable of the data-generating process
#'
#' Describes the marginal distribution of one regressor (the focal variable
#' or a covariate) for [ape_dgp()]. Dependence between variables is set via
#' the `correlation` argument of [ape_dgp()] (a Gaussian copula).
#'
#' @param name Variable name (a string; used in printing and error messages).
#' @param type `"binary"` (Bernoulli) or `"normal"`.
#' @param p Success probability, strictly between 0 and 1 (binary only).
#' @param mean,sd Mean and standard deviation (normal only; `sd > 0`).
#' @param icc Panel designs only ([ape_dgp_panel()]): within-unit correlation
#'   of the variable over time, via a shared-component device on the latent
#'   normal scale (`u_it = sqrt(icc) v_i + sqrt(1-icc) e_it`) before the
#'   marginal transform. `icc = 0` redraws each period; `icc = 1` makes the
#'   variable time-constant (e.g., a unit-level treatment). For binary
#'   variables the realized outcome-scale correlation is attenuated relative
#'   to `icc`, as with the Gaussian-copula `correlation`. Ignored by the
#'   cross-sectional constructors.
#'
#' @return An object of class `pa_var`.
#' @examples
#' pa_var("treat", "binary", p = 0.5)
#' pa_var("age", "normal", mean = 45, sd = 12)
#' @export
pa_var <- function(name, type = c("binary", "normal"), p = NULL, mean = 0, sd = 1,
                   icc = NULL) {
  type <- match.arg(type)
  stopifnot(is.character(name), length(name) == 1L, nzchar(name))
  if (type == "binary") {
    if (is.null(p) || !is.numeric(p) || length(p) != 1L || p <= 0 || p >= 1)
      stop("A binary pa_var needs `p` strictly between 0 and 1.", call. = FALSE)
  } else {
    stopifnot(is.numeric(mean), length(mean) == 1L,
              is.numeric(sd), length(sd) == 1L, sd > 0)
  }
  if (!is.null(icc))
    stopifnot(is.numeric(icc), length(icc) == 1L, icc >= 0, icc <= 1)
  structure(list(name = name, type = type, p = p, mean = mean, sd = sd,
                 icc = icc),
            class = "pa_var")
}

# Marginal standard deviation of a pa_var.
pa_var_sd <- function(v) if (v$type == "binary") sqrt(v$p * (1 - v$p)) else v$sd

# Map a standard-normal draw to the variable's marginal (Gaussian copula).
# Binary is increasing in u so a positive copula correlation means positive
# association: P(u > qnorm(1 - p)) = p.
transform_u <- function(u, v) {
  if (v$type == "binary") as.numeric(u > qnorm(1 - v$p)) else v$mean + v$sd * u
}
