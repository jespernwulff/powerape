* Runs behind the Stata tabs of Design E in "Worked example: five designs":
* a treatment that firms choose (a training program), instrumented by a
* subsidy offer and analysed by a control-function probit. The single power
* calls use the default engine (-cfprobit-, or -probit- for the exogenous
* comparison, on every simulated study; Stata 18 or later), the searches and
* the sweep the Mata engine. Start from this folder (see _prelude.do).
do "_prelude.do" worked-example-e

* >>> world
powerape dgp e, iv focal(training, binary p(0.4)) covariates((size, normal)) ///
    instruments((subsidy, binary p(0.5))) endogeneity(0.4) ivstrength(0.2)   ///
    baseline(0.30) signal(0.10)
powerape setape e, target(0.10)
* <<<

* >>> exog-world
powerape dgp ex, focal(training, binary p(0.4)) covariates((size, normal))  ///
    baseline(0.30) signal(0.10)
powerape setape ex, target(0.10)
* <<<

tic
* >>> exog-power
powerape power ex, n(1200) claim(detect) nsim(2000) seed(601)
* <<<
keep_r exog_power
toc exog_power

tic
* >>> power
powerape power e, n(1200) claim(detect) nsim(2000) seed(602)
* <<<
keep_r power
toc power

tic
* >>> n
powerape n e, power(0.80) claim(detect) nsim(2000) nsimconfirm(16000)      ///
    seed(603) engine(mata)
* <<<
keep_r n
toc n

tic
* >>> exog-n
powerape n ex, power(0.80) claim(detect) nsim(1000) seed(606) engine(mata)
* <<<
keep_r exog_n
toc exog_n

tic
* >>> robust
powerape robust e, n(1200) claim(detect) vary(ivstrength(0.1 0.3))         ///
    gridpoints(3) nsim(800) seed(604) nonmax engine(mata)
* <<<
keep_range robust
toc robust

* >>> aie-world
powerape dgp ei, iv focal(training, binary p(0.4))                         ///
    moderator(public, binary p(0.5)) covariates((size, normal))            ///
    instruments((subsidy, binary p(0.5))) endogeneity(0.4)                 ///
    ivstrength(0.25) baseline(0.30) signal(0.10)
powerape setaie ei, target(0.06) mainfocal(0.10) mainmoderator(0.05)
* <<<

tic
* >>> aie-power
powerape power ei, n(2000) claim(detect) nsim(1000) seed(605)
* <<<
keep_r aie_power
toc aie_power

file close res
file close tim
log close pc
exit, clear
