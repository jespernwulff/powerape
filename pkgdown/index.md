<!-- Generated from index.Rmd.orig by precompute/precompile.R. Edit the .orig file, not this one. -->



# powerape

<p class="lead">Power analysis when the effect size is an average partial effect or an average interaction effect of a probit or logit model, stated in probability points. An R package and a Stata command.</p>

<p class="home-actions">
<a class="btn btn-primary" href="articles/minimum-effect.html">Get started in R</a>
<a class="btn btn-primary" href="articles/stata.html">Get started in Stata</a>
<a class="btn btn-outline-primary" href="articles/index.html">Workflows</a>
</p>

## What it does

In a probit or logit model the coefficient is not the quantity anyone
interprets. Conclusions are stated in probability points: the program
raises firm survival by ten percentage points. That quantity is the
**average partial effect (APE)**, and for a moderation question the
**average interaction effect (AIE)**. Both depend on all coefficients and
on the distribution of the covariates, so power tools that work in
coefficient units cannot say how many observations a study needs to find
them.

`powerape` takes the effect in those units. You describe the world on
friendly scales (a baseline outcome rate and the share of latent variation
the covariates explain), pin the effect you expect, and the package
simulates the whole study: it draws the data, fits the model, computes the
APE or AIE with a delta-method confidence interval, and checks a claim on
that interval.

## Why it exists

Powering a study to *detect* an effect is not the same as powering it to
show that the effect is *large enough to matter*. Following Riesthuis
(2024), `powerape` prices three claims, each read off a confidence interval:
that the effect exceeds the smallest effect size of interest (the SESOI;
the default), that it differs from zero, or that it is too small to matter.
It reports how the simulated studies end, so the share that would detect the
effect without being able to call it meaningful is in plain view. Around
this sit the tools a design discussion needs: required sample sizes and
minimum detectable effects, robustness sweeps over the contextual guesses
(Hancock & Feng, 2026), justified error rates (Maier & Lakens, 2022; Wulff
& Taylor, 2024), and a citable methods paragraph. The R package and the
Stata command implement the same procedure and are documented side by side
on this site.

## Install

<div class="lang-tabs">
<ul class="nav nav-tabs" role="tablist">
<li class="nav-item" role="presentation"><button class="nav-link active" id="install-r-tab" data-bs-toggle="tab" data-bs-target="#install-r" type="button" role="tab" aria-controls="install-r" aria-selected="true" data-lang="r">R</button></li>
<li class="nav-item" role="presentation"><button class="nav-link" id="install-stata-tab" data-bs-toggle="tab" data-bs-target="#install-stata" type="button" role="tab" aria-controls="install-stata" aria-selected="false" data-lang="stata">Stata</button></li>
</ul>
<div class="tab-content">
<div class="tab-pane active" id="install-r" role="tabpanel" aria-labelledby="install-r-tab">

```r
# install.packages("remotes")
remotes::install_github("jespernwulff/powerape")
```

</div>
<div class="tab-pane" id="install-stata" role="tabpanel" aria-labelledby="install-stata-tab">

```stata
net install powerape, from("https://raw.githubusercontent.com/jespernwulff/powerape-stata/main/")
mata: mata mlib index
```

</div>
</div>
</div>

## One design question

A support program might raise firm survival, which runs near 30 percent
without it, by ten percentage points; it would justify its cost from five
points upward. Covariates explain a modest share of the latent variation.
How many firms must a study observe to show that the program is worth its
cost, and how many to show merely that it does something?

<div class="lang-tabs">
<ul class="nav nav-tabs" role="tablist">
<li class="nav-item" role="presentation"><button class="nav-link active" id="design-r-tab" data-bs-toggle="tab" data-bs-target="#design-r" type="button" role="tab" aria-controls="design-r" aria-selected="true" data-lang="r">R</button></li>
<li class="nav-item" role="presentation"><button class="nav-link" id="design-stata-tab" data-bs-toggle="tab" data-bs-target="#design-stata" type="button" role="tab" aria-controls="design-stata" aria-selected="false" data-lang="stata">Stata</button></li>
</ul>
<div class="tab-content">
<div class="tab-pane active" id="design-r" role="tabpanel" aria-labelledby="design-r-tab">


``` r
library(powerape)

d <- ape_dgp(focal = pa_var("program", "binary", p = 0.5),
             covariates = list(pa_var("z", "normal")),
             baseline = 0.30, signal = 0.10)
d <- set_ape(d, target = 0.10)

n_min <- ape_n(d, power = 0.80, claim = "minimum", sesoi = 0.05, nsim = 1000, seed = 502)
n_det <- ape_n(d, power = 0.80, claim = "detect", nsim = 1000, seed = 503)
n_min
#> powerape required sample size -- minimum claim
#>   n = 2091 for 80% target power (confirmed 0.808, MCSE 0.006)
#>   assumed true APE +0.1000, sesoi 0.050, 90% CI (alpha = 0.05), probit
#>   search: 1 step(s); confirmed in 1 round(s) at nsim = 4000.
n_det
#> powerape required sample size -- detect claim
#>   n = 663 for 80% target power (confirmed 0.803, MCSE 0.006)
#>   assumed true APE +0.1000, sesoi -, 95% CI (alpha = 0.05), probit
#>   search: 1 step(s); confirmed in 1 round(s) at nsim = 4000.
```

</div>
<div class="tab-pane" id="design-stata" role="tabpanel" aria-labelledby="design-stata-tab">

```stata
powerape dgp program, focal(program, binary p(0.5)) covariates((z, normal)) ///
    baseline(0.30) signal(0.10)
powerape setape program, target(0.10)

powerape n program, power(0.80) claim(minimum) sesoi(0.05)               ///
    nsim(1000) seed(502) engine(mata)

powerape n program, power(0.80) claim(detect)                            ///
    nsim(1000) seed(503) engine(mata)
```

```{.stata-output}
powerape required sample size -- minimum claim
  n = 2084 for 80% target power (confirmed 0.800, MCSE 0.006)
  assumed true APE + 0.1000, sesoi 0.050,     90% CI (alpha =    .05), probit
  fit: Mata engine (Fisher scoring)
  search: 1 step(s); confirmed in 1 round(s) at nsim = 4000.

powerape required sample size -- detect claim
  n = 663 for 80% target power (confirmed 0.796, MCSE 0.006)
  assumed true APE + 0.1000, sesoi -,     95% CI (alpha =    .05), probit
  fit: Mata engine (Fisher scoring)
  search: 1 step(s); confirmed in 1 round(s) at nsim = 4000.
```

</div>
</div>
</div>

Showing that the program is worth its cost takes about
2,091 firms; showing that it does something takes
663 (Stata: 2,084 and
663; R and Stata draw different random numbers,
so the two agree up to simulation error). A study planned for detection
would be about a third of the size it needs for the conclusion its authors
want to draw. The [worked example](articles/worked-example.html) takes this
design all the way through, together with four more: an interaction, a
design built from a pilot, a firm-year panel and an instrumented treatment.

## What you get

<div class="home-features">

**Effects in the units you report.** APEs of binary and continuous
variables and AIEs, the estimand of `ginteff`, pinned in probability
points; impossible targets are refused with the attainable range.

**Three claims, and how studies end.** Minimum-effect, detection and
equivalence power, with the share of studies that detect without being
able to call the effect meaningful.

**Sample size, detectable effect, curves.** Required n and minimum
detectable effects, each confirmed by a high-precision run, and power
curves.

**Insurance against the guesses.** Sweeps over the baseline rate, the
covariate signal and design parameters, with the worst case and n_max.

**Justified error rates.** Alpha as a function of n, or the error-cost
optimum, computed from one stored simulation, with a size check.

**Real designs.** Pilot data or a fitted pilot model, correlated random
effects panels, and endogenous treatments by control functions; in Stata,
every simulated study can be fitted by `probit`, `logit` or `cfprobit`.

</div>





<figure class="home-figure">
<img src="man/figures/home-curves.png" alt="Power against the sample size for the detection claim and the minimum-effect claim in the design above. Detection reaches 80 percent power at about 660 firms, the minimum-effect claim at about 2,100." />
<figcaption>The price of a meaningful conclusion. Power of the detection claim and of the minimum-effect claim against the sample size in the design above; the dashed line marks 80 percent.</figcaption>
</figure>

## Where to go from here

<div class="home-cards">
<a class="home-card" href="articles/minimum-effect.html"><span class="home-card-title">Get started in R</span><span class="home-card-text">Describe the world, pin the APE, and interrogate the design in six steps.</span></a>
<a class="home-card" href="articles/stata.html"><span class="home-card-title">Get started in Stata</span><span class="home-card-text">The same six steps with the Stata command, plus engines and pilot data.</span></a>
<a class="home-card" href="articles/worked-example.html"><span class="home-card-title">Worked example</span><span class="home-card-text">Five designs in both languages, from a treatment effect to an instrument.</span></a>
<a class="home-card" href="articles/classic-sesoi-detection.html"><span class="home-card-title">The classic SESOI analysis</span><span class="home-card-text">What powering detection at the SESOI buys, and what it cannot.</span></a>
<a class="home-card" href="articles/equivalence.html"><span class="home-card-title">Equivalence</span><span class="home-card-text">Showing that an effect is too small to matter.</span></a>
<a class="home-card" href="articles/interaction-effects.html"><span class="home-card-title">Interaction effects</span><span class="home-card-text">Power for the AIE, and why interactions are expensive.</span></a>
<a class="home-card" href="articles/panel-designs.html"><span class="home-card-title">Panel designs</span><span class="home-card-text">Correlated random effects probit, with n counted in units.</span></a>
<a class="home-card" href="articles/iv-designs.html"><span class="home-card-title">IV designs</span><span class="home-card-text">Endogenous treatments by control functions, and the price of instruments.</span></a>
<a class="home-card" href="articles/justified-alpha.html"><span class="home-card-title">Justified alpha</span><span class="home-card-text">Alpha by sample size or by error costs, from one simulation.</span></a>
<a class="home-card" href="articles/r-and-stata.html"><span class="home-card-title">R and Stata side by side</span><span class="home-card-text">Every function and subcommand, matched.</span></a>
</div>

## Citation

If you use `powerape`, please cite the software, and the methods it builds
on where you use them.

Wulff, J. N. (2026). *powerape: Power analysis for average partial effects*
(R package version 1.11.4; Stata command version
0.1.3). <https://jespernwulff.github.io/powerape/>

Hancock, G. R., & Feng, Y. (2026). nmax and the quest to restore caution,
integrity, and practicality to the sample size planning process.
*Psychological Methods*, *31*(1), 174--199. <https://doi.org/10.1037/met0000776>

Maier, M., & Lakens, D. (2022). Justify your alpha: A primer on two
practical approaches. *Advances in Methods and Practices in Psychological
Science*, 5(2). <https://doi.org/10.1177/25152459221080396>

Riesthuis, P. (2024). Simulation-based power analyses for the smallest
effect size of interest: A confidence-interval approach for minimum-effect
and equivalence testing. *Advances in Methods and Practices in
Psychological Science*, 7(2). <https://doi.org/10.1177/25152459241240722>

Wulff, J. N., & Taylor, L. (2024). How and why alpha should depend on
sample size: A Bayesian-frequentist compromise for significance testing.
*Strategic Organization*, 22(3), 550--581.
<https://doi.org/10.1177/14761270231214429>
