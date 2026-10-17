* Runs behind the Stata tabs of "Panel designs: correlated random effects".
* The default engine fits every simulated panel with -probit- and
* -vce(cluster)- on the unit identifier. Start from this folder (see
* _prelude.do).
do "_prelude.do" panel

* >>> world
powerape dgp d, panel focal(treat, binary p(0.5) icc(1))                   ///
    covariates((size, normal icc(0.6))) nperiods(4) rho(0.30) creshare(0.50) ///
    baseline(0.30) signal(0.10)
* <<<

* >>> pin
powerape setape d, target(0.10)
* <<<

tic
* >>> power
powerape power d, n(300) claim(minimum) sesoi(0.05) nsim(400) seed(1)
* <<<
keep_r power
toc power

tic
* >>> mde
powerape mde d, n(400) claim(minimum) sesoi(0.05) nsim(400) seed(9) engine(mata)
* <<<
keep_r mde
toc mde

* >>> cs-world
powerape dgp cs, focal(treat, binary p(0.5)) covariates((size, normal))   ///
    baseline(0.30) signal(0.10)
powerape setape cs, target(0.10)
* <<<

tic
* >>> cs-power
powerape power cs, n(1200) claim(minimum) sesoi(0.05) nsim(400) seed(2)
* <<<
keep_r cs_power
toc cs_power

tic
* >>> rho
powerape robust d, n(300) claim(minimum) sesoi(0.05) vary(rho(0.1 0.5))    ///
    gridpoints(3) nsim(300) seed(3) nonmax engine(mata)
* <<<
keep_r rho
toc rho

* >>> unbalanced-world
powerape dgp da, panel focal(treat, binary p(0.5) icc(1))                  ///
    covariates((size, normal icc(0.6))) nperiods(6) retention(0.85)          ///
    rho(0.30) creshare(0.50) baseline(0.30) signal(0.10)
* <<<

* >>> unbalanced-pin
powerape setape da, target(0.10)
* <<<

tic
* >>> retention
powerape robust da, n(300) claim(minimum) sesoi(0.05) vary(retention(0.7 1)) ///
    gridpoints(3) nsim(300) seed(5) nonmax engine(mata)
* <<<
keep_r retention
toc retention

* >>> aie-world
powerape dgp di, panel focal(treat, binary p(0.5) icc(1))                  ///
    moderator(public, binary p(0.4) icc(0.7))                              ///
    covariates((size, normal icc(0.6))) nperiods(4) rho(0.30) creshare(0.50) ///
    baseline(0.30) signal(0.10)
* <<<

* >>> aie-pin
powerape setaie di, target(0.06) mainfocal(0.10) mainmoderator(0.05)
* <<<

tic
* >>> aie-power
powerape power di, n(400) claim(detect) nsim(250) seed(4)
* <<<
keep_r aie_power
toc aie_power

file close res
file close tim
log close pc
exit, clear
