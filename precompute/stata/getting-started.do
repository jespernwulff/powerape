* Runs behind "Get started in Stata": the six steps of the R article
* "Powering the minimum-effect claim" in the same world, then the engines and
* the pilot-data routes. Start from this folder (see _prelude.do).
do "_prelude.do" getting-started

* ---- The six steps -----------------------------------------------------------

* >>> world
powerape dgp d, focal(treat, binary p(0.5))                              ///
    covariates((age, normal mean(45) sd(12)) (female, binary p(0.55)))   ///
    correlation(0.2) baseline(0.30) signal(0.15)
* <<<

* >>> pin
powerape setape d, target(0.10)
* <<<

tic
* >>> power
powerape power d, n(1500) claim(minimum) sesoi(0.05) nsim(500) seed(1)
* <<<
keep_r power

* >>> stored
return list
* <<<
toc power

tic
* >>> n
powerape n d, power(0.80) claim(minimum) sesoi(0.05) nsim(500) seed(2) engine(mata)
* <<<
keep_r n
toc n

tic
* >>> curve
powerape curve d, n(500(500)3000) claim(minimum) sesoi(0.05)             ///
    nsim(300) seed(3) engine(mata)
* <<<
keep_r curve
toc curve

* >>> curve-plot
matrix C = r(results)
svmat C, names(col)
twoway connected power n, yline(0.8, lpattern(dash))                      ///
    ytitle("Power") xtitle("Sample size") ylabel(0(0.2)1)
* <<<
pc_graph getting-started-curve
drop _all

tic
* >>> mde
powerape mde d, n(1500) claim(minimum) sesoi(0.05) nsim(500) seed(6) engine(mata)
* <<<
keep_r mde
toc mde

tic
* >>> robust
powerape robust d, n(2200) claim(minimum) sesoi(0.05)                    ///
    vary(baseline(0.20 0.40)) gridpoints(3) nsim(300) seed(4) nonmax engine(mata)
* <<<
keep_r robust
toc robust

tic
* >>> statement
powerape power d, n(2200) claim(minimum) sesoi(0.05) nsim(500) seed(5)
powerape statement
* <<<
keep_r statement
toc statement

* ---- Engines -------------------------------------------------------------------

tic
* >>> engine-mata
powerape power d, n(1500) claim(minimum) sesoi(0.05) nsim(500) seed(1) engine(mata)
* <<<
keep_r engine_mata
toc engine_mata

tic
* >>> effects-closed
powerape power d, n(1500) claim(minimum) sesoi(0.05) nsim(100) seed(7)
* <<<
keep_r effects_closed
toc effects_closed

tic
* >>> effects-margins
powerape power d, n(1500) claim(minimum) sesoi(0.05) nsim(100) seed(7)   ///
    effects(margins)
* <<<
keep_r effects_margins
toc effects_margins

* >>> session
powerape engine
powerape engine, engine(mata)
powerape engine, engine(stata) effects(closed)
* <<<

* ---- Pilot data ------------------------------------------------------------------

* >>> pilot-data
set seed 20261016
set obs 500
generate treat  = runiform() < .5
generate age    = rnormal(45, 12)
generate female = runiform() < .55
generate y      = runiform() < normal(-0.9 + 0.3*treat + 0.01*age + 0.2*female)
* <<<

tic
* >>> pilot-empirical
powerape dgp pilot_rows, empirical focal(treat) covariates(age female)    ///
    baseline(0.30) signal(0.15)
powerape setape pilot_rows, target(0.10)
powerape power pilot_rows, n(1500) claim(minimum) sesoi(0.05) nsim(500) seed(11)
* <<<
keep_r pilot_empirical
toc pilot_empirical

tic
* >>> pilot-fromfit
probit y treat age female, nolog
powerape dgp pilot_fit, fromfit focal(treat)
powerape setape pilot_fit, target(0.10)
powerape power pilot_fit, n(1500) claim(minimum) sesoi(0.05) nsim(500) seed(12)
* <<<
keep_r pilot_fromfit
toc pilot_fromfit

* >>> list
powerape list
* <<<

file close res
file close tim
log close pc
exit, clear
