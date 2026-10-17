* Runs behind the Stata tabs of "IV designs: endogenous focal variables via
* control functions". The default engine fits every simulated study with
* Stata's own -cfprobit- (Stata 18 or later). Start from this folder (see
* _prelude.do).
do "_prelude.do" iv

* >>> world
powerape dgp d, iv focal(training, binary p(0.4)) covariates((size, normal)) ///
    instruments((subsidy, binary p(0.5))) endogeneity(0.4) ivstrength(0.2)   ///
    baseline(0.30) signal(0.10)
* <<<

* >>> pin
powerape setape d, target(0.10)
* <<<

tic
* >>> power
powerape power d, n(1200) claim(detect) nsim(400) seed(1)
* <<<
keep_r power
toc power

* >>> exog-world
powerape dgp ds, focal(training, binary p(0.4)) covariates((size, normal)) ///
    baseline(0.30) signal(0.10)
powerape setape ds, target(0.10)
* <<<

tic
* >>> exog-power
powerape power ds, n(1200) claim(detect) nsim(400) seed(2)
* <<<
keep_r exog_power
toc exog_power

tic
* >>> robust
powerape robust d, n(1200) claim(detect) vary(ivstrength(0.1 0.3)) gridpoints(3) ///
    nsim(300) seed(3) nonmax engine(mata)
* <<<
keep_r robust
toc robust

* >>> aie-world
powerape dgp di, iv focal(training, binary p(0.4))                       ///
    moderator(public, binary p(0.5)) covariates((size, normal))           ///
    instruments((subsidy, binary p(0.5))) endogeneity(0.4)                ///
    ivstrength(0.25) baseline(0.30) signal(0.10)
powerape setaie di, target(0.06) mainfocal(0.10) mainmoderator(0.05)
* <<<

tic
* >>> aie-power
powerape power di, n(2000) claim(detect) nsim(250) seed(4)
* <<<
keep_r aie_power
toc aie_power

file close res
file close tim
log close pc
exit, clear
