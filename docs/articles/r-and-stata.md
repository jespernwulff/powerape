# R and Stata side by side

The R package and the Stata command implement the same procedure: the
same worlds, the same inversion of a target effect into coefficients,
the same claims, searches and sweeps. This page matches their functions,
arguments, subcommands and options, for readers who move between the two
or who want to translate an example.

## The same analysis in both

In R:

``` r

library(powerape)

d <- ape_dgp(focal = pa_var("program", "binary", p = 0.5),
             covariates = list(pa_var("z", "normal")),
             baseline = 0.30, signal = 0.10)
d <- set_ape(d, target = 0.10)

pw <- ape_power(d, n = 1500, claim = "minimum", sesoi = 0.05, seed = 1)
ape_n(d, power = 0.80, claim = "minimum", sesoi = 0.05, seed = 2)
power_statement(pw)
```

In Stata:

``` stata
powerape dgp d, focal(program, binary p(0.5)) covariates((z, normal)) ///
    baseline(0.30) signal(0.10)
powerape setape d, target(0.10)

powerape power d, n(1500) claim(minimum) sesoi(0.05) seed(1)
powerape n d, power(0.80) claim(minimum) sesoi(0.05) seed(2) engine(mata)
powerape statement
```

Three things differ by design. The answers agree only up to simulation
error, because R and Stata have different random-number generators and
the same seed does not give the same simulated studies. In R a design
and a result are objects that you keep, print and pass on, whereas in
Stata they are stored by name, like `estimates store`, and the results
of the last command are in `r()`. And in Stata you choose who fits the
simulated studies: Stata’s own `probit`, `logit` or `cfprobit`
(`engine(stata)`, the default) or a Mata port of the R package’s routine
(`engine(mata)`), which is many times faster and is the one to use for
searches.

## Describing the world

| Purpose | R | Stata |
|----|----|----|
| Probit or logit | `model = "probit"` (default) or `"logit"` | `model(probit)` (default) or `model(logit)` |
| A binary variable | `pa_var("treat", "binary", p = 0.5)` | `(treat, binary p(0.5))` |
| A normal variable | `pa_var("age", "normal", mean = 45, sd = 12)` | `(age, normal mean(45) sd(12))` |
| Focal variable | `focal = pa_var(...)` | `focal(name, ...)` |
| Covariates | `covariates = list(pa_var(...), pa_var(...))` | `covariates((name, ...) (name, ...))` |
| Moderator (interaction designs) | `moderator = pa_var(...)` | `moderator(name, ...)` |
| Latent correlation | `correlation = 0.2`, or a matrix | `correlation(0.2)`, or the name of a matrix |
| Baseline outcome rate | `baseline = 0.30` | `baseline(0.30)` |
| Covariate signal (latent pseudo-R²) | `signal = 0.15` | `signal(0.15)` |
| Explicit nuisance coefficients | `gamma = c(...)` | `gamma(numlist)` |
| Integration draws for the truth | `n_int`, `seed_int` | `nint()`, `seedint()` |
| Covariates resampled from pilot data | `ape_dgp_empirical(data = pilot, focal = "treat", ...)` | `powerape dgp name, empirical focal(treat) covariates(age female) ...`, on the data in memory |
| The world of a fitted pilot model | `ape_dgp_from_fit(fit, focal = "treat")` | `probit y treat age female`, then `powerape dgp name, fromfit focal(treat)` |
| Panel: periods, persistence, unit effect | `ape_dgp_panel(n_periods = 4, rho = 0.3, cre_share = 0.5, ...)`, `pa_var(..., icc = 0.6)` | `powerape dgp name, panel nperiods(4) rho(0.3) creshare(0.5) ...`, `icc(0.6)` in a variable |
| Unbalanced panel | `n_periods = c(2, 4, 6)` with `p_periods`, or `retention = 0.85` | `nperiods(2 4 6)` with `pperiods()`, or `retention(0.85)` |
| Endogenous focal variable | `ape_dgp_iv(instruments = ..., endogeneity = 0.4, iv_strength = 0.2, ...)` | `powerape dgp name, iv instruments(...) endogeneity(0.4) ivstrength(0.2) ...` |
| Pin the APE | `d <- set_ape(d, target = 0.10)` | `powerape setape d, target(0.10)` |
| Pin the AIE | `set_aie(d, target = 0.08, main_focal = 0.10, main_moderator = 0.03)` | `powerape setaie d, target(0.08) mainfocal(0.10) mainmoderator(0.03)` |
| The true effect of a world | `true_ape(d)`, `true_aie(d)` | `powerape trueape d`, `powerape trueaie d` |
| Show a world | `print(d)` | `powerape describe d` |

[`set_ape()`](https://jespernwulff.github.io/powerape/reference/set_ape.md)
returns a new design and leaves `d` as it was unless you assign the
result back; `powerape setape` changes the stored design in place. To
keep two versions in Stata, describe the world twice under two names.

## Power, searches and sweeps

| Purpose | R | Stata |
|----|----|----|
| Power at one sample size | `ape_power(d, n = 1500, ...)` | `powerape power d, n(1500) ...` |
| Claim | `claim = "minimum"` (default), `"detect"`, `"equivalence"` | `claim(minimum)` (default), `claim(detect)`, `claim(equivalence)` |
| Smallest effect size of interest | `sesoi = 0.05` | `sesoi(0.05)` |
| Error rate | `alpha = 0.05` (default) | `alpha(0.05)` (default) |
| Error rate as a function of n | `alpha = function(n) ...` | `alpharule(expression in @n)` |
| Interval level, overriding the claim’s | `conf = 0.95` | `conf(0.95)` |
| Replications, seed | `nsim = 1000`, `seed = 1` | `nsim(1000)`, `seed(1)` |
| Robust standard errors | `se = "robust"` | `se(robust)` |
| Separated replications | `separation = "fail"` (default) or `"keep"` | `separation(fail)` (default) or `separation(keep)` |
| Keep the replication draws (for alpha) | `keep_draws = TRUE` (default) | kept by default; `nodraws` |
| Power curve | `ape_curve(d, n = seq(500, 3000, by = 500), ...)` | `powerape curve d, n(500(500)3000) ...` |
| Required sample size | `ape_n(d, power = 0.80, ...)`; default power 0.90 | `powerape n d, power(0.80) ...`; default power 0.90 |
| Search settings | `n_range`, `max_iter`, `confirm = FALSE`, `nsim_confirm` | `nrange()`, `maxiter()`, `noconfirm`, `nsimconfirm()` |
| Minimum detectable effect | `ape_mde(d, n = 1500, ...)`, `direction`, `main_focal`, `main_moderator` | `powerape mde d, n(1500) ...`, `direction()`, `mainfocal()`, `mainmoderator()` |
| Robustness sweep | `ape_robust(d, n = 2200, vary = list(baseline = c(0.2, 0.4)), ...)` | `powerape robust d, n(2200) vary(baseline(0.2 0.4)) ...` |
| Sweep settings | `pin`, `mode = "mde"`, `grid_points`, `nmax = FALSE`, `nmax_power` | `pin()`, `mode(mde)`, `gridpoints()`, `nonmax`, `nmaxpower()` |
| Who fits the simulated studies | the package’s own routine | `engine(stata)` (default) or `engine(mata)`; `effects(closed)` (default) or `effects(margins)` |

## Reporting and justified error rates

| Purpose | R | Stata |
|----|----|----|
| Citable methods paragraph | `power_statement(pw)` | `powerape statement` (the last result) or `powerape statement name` |
| Power at other error rates | `power_at(pw, alpha = c(0.01, 0.05))` | `powerape powerat name, alpha(0.01 0.05)` |
| Error-cost optimum | `ape_alpha(pw, cost = 4, prior = 1)` | `powerape alpha name, cost(4) prior(1)` |
| Size at the claim boundary | `size = TRUE, nsim_size = 500` | `size nsimsize(500)` |
| Optimum at the detectability frontier | `alpha_frontier(cost = 4, power = 0.80)` | `powerape frontier, cost(4) power(0.80)` |
| Plots | `plot(curve)`, `plot(alpha)` | `twoway` on `r(results)` or `r(grid)` |
| Results | elements of the returned object (`$power`, `$n`, `$mde`, …) | `r()`; `name()` stores a result; `powerape list`, `powerape drop` |

## Differences worth knowing

**Variance estimators.** The R package and `engine(mata)` compute model
standard errors from the expected information, robust ones with the HC0
sandwich and clustered ones with the G/(G − 1) factor. `engine(stata)`
uses what Stata reports: the observed information, and its small-sample
factors for `vce(robust)` and `vce(cluster)`. On fixed data the probit
standard errors differ by at most about 0.00006 and the logit ones agree
to numerical precision, so the two reach the same decision in nearly
every simulated study.

**Alpha rules.** In R a rule is any function of *n*, for example
`function(n) alphaN::alphaN(n, BF = 3)` from the alphaN package. In
Stata it is an expression in `@n`; alphaN’s rule with `BF = 3` is
`alpharule(chi2tail(1, ln(@n) + 2*ln(3)))`, which reproduces the R
function to machine precision. [Justified
alpha](https://jespernwulff.github.io/powerape/articles/justified-alpha.md)
uses both.

**Interaction effects with `effects(margins)`.** In Stata the AIE of
every simulated study can be computed by `ginteff` itself (Radean,
2023), which must then be installed from the Stata Journal:
`net install st0711, from(http://www.stata-journal.com/software/sj23-2)`.
On the control-function route the AIE is computed by `margins` after
`cfprobit`, with the control function held fixed.

**Stata versions.** The command needs Stata 16 or later; the
control-function route under the default engine uses `cfprobit`, which
needs Stata 18 or later (`engine(mata)` covers it on any version).

**Results are stored by name.** Every DGP and every result stays
available until you drop it, and the subcommands that work on a result
(`statement`, `powerat`, `alpha`) take its name; by default they use the
last result.

## How the two agree

The Stata command is validated against the R package in four tiers:
exact power by enumeration, the estimators on fixed data sets (and
against `probit`, `logit`, `margins`, `ginteff` and `cfprobit`), the
calibration constants, and simulated power, sample sizes and detectable
effects in matched designs. See [How powerape is
validated](https://jespernwulff.github.io/powerape/articles/validation.md).

## References

Radean, M. (2023). ginteff: A generalized command for computing
interaction effects. *The Stata Journal*, 23(2), 301–335.
<https://doi.org/10.1177/1536867X231175253>
