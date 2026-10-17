* Runs behind the Stata tabs of Design D in "Worked example: five designs":
* a firm-year panel (four waves, a firm-level adoption), analysed by a
* correlated random effects probit with firm-clustered standard errors. The
* single power calls use the default engine (-probit- with -vce(cluster)- on
* every simulated panel), the searches and the sweep the Mata engine. The
* single power calls use 4,000 and 2,000 replications where the article's R
* code uses 20,000 and 10,000: under -probit- each simulated panel takes
* Stata far longer than R, and those fits take most of this run's time (see
* precompute/README.md). Start from this folder (see _prelude.do).
do "_prelude.do" worked-example-d

* >>> world
powerape dgp d, panel focal(adopt, binary p(0.5) icc(1))                   ///
    covariates((size, normal icc(0.6))) nperiods(4) rho(0.30)              ///
    creshare(0.50) baseline(0.30) signal(0.10)
powerape setape d, target(0.10)
* <<<

tic
* >>> power
powerape power d, n(300) claim(minimum) sesoi(0.05) nsim(4000) seed(501)
* <<<
keep_r power
toc power

* >>> cs-world
powerape dgp cs, focal(adopt, binary p(0.5)) covariates((size, normal))    ///
    baseline(0.30) signal(0.10)
powerape setape cs, target(0.10)
* <<<

tic
* >>> cs-power
powerape power cs, n(1200) claim(minimum) sesoi(0.05) nsim(4000) seed(502)
* <<<
keep_r cs_power
toc cs_power

tic
* >>> n
powerape n d, power(0.80) claim(minimum) sesoi(0.05) nsim(1000) seed(503)  ///
    engine(mata)
* <<<
keep_r n
toc n

tic
* >>> rho
powerape robust d, n(300) claim(minimum) sesoi(0.05) vary(rho(0.1 0.5))    ///
    gridpoints(3) nsim(8000) seed(504) nonmax engine(mata)
* <<<
keep_range rho
toc rho

* >>> attrition-worlds
powerape dgp d85, panel focal(adopt, binary p(0.5) icc(1))                 ///
    covariates((size, normal icc(0.6))) nperiods(4) retention(0.85)        ///
    rho(0.30) creshare(0.50) baseline(0.30) signal(0.10) nint(40000)
powerape setape d85, target(0.10)
powerape dgp d70, panel focal(adopt, binary p(0.5) icc(1))                 ///
    covariates((size, normal icc(0.6))) nperiods(4) retention(0.70)        ///
    rho(0.30) creshare(0.50) baseline(0.30) signal(0.10) nint(40000)
powerape setape d70, target(0.10)
* <<<

tic
* >>> power-85
powerape power d85, n(300) claim(minimum) sesoi(0.05) nsim(2000) seed(71)
* <<<
keep_r power_85
toc power_85

tic
* >>> power-70
powerape power d70, n(300) claim(minimum) sesoi(0.05) nsim(2000) seed(71)
* <<<
keep_r power_70
toc power_70

tic
* >>> n-85
powerape n d85, power(0.80) claim(minimum) sesoi(0.05) nsim(800) seed(72)  ///
    engine(mata)
* <<<
keep_r n_85
toc n_85

tic
* >>> n-70
powerape n d70, power(0.80) claim(minimum) sesoi(0.05) nsim(800) seed(72)  ///
    engine(mata)
* <<<
keep_r n_70
toc n_70

file close res
file close tim
log close pc
exit, clear
