* Runs behind the Stata tabs of "Powering for equivalence: showing an effect is
* too small to matter". Start from this folder (see _prelude.do).
do "_prelude.do" equivalence

* >>> world
powerape dgp d, focal(treat, binary p(0.5))                               ///
    covariates((age, normal mean(45) sd(12))) baseline(0.30) signal(0.10)
powerape setape d, target(0)
* <<<

tic
* >>> power
powerape power d, n(2000) claim(equivalence) sesoi(0.05) nsim(600) seed(1)
* <<<
keep_r power
toc power

tic
* >>> n
powerape n d, power(0.80) claim(equivalence) sesoi(0.05)                 ///
    nsim(500) seed(2) engine(mata)
* <<<
keep_r n
toc n

tic
* >>> curve
powerape curve d, n(500(750)3500) claim(equivalence) sesoi(0.05)         ///
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
pc_graph equivalence-curve
drop _all

tic
* >>> robust
powerape robust d, n(2600) claim(equivalence) sesoi(0.05)                 ///
    vary(baseline(0.20 0.50)) gridpoints(3) nsim(300) seed(4) nonmax engine(mata)
* <<<
keep_r robust
toc robust

tic
* >>> statement
powerape power d, n(2600) claim(equivalence) sesoi(0.05) nsim(500) seed(5)
powerape statement
* <<<
keep_r statement
toc statement

file close res
file close tim
log close pc
exit, clear
