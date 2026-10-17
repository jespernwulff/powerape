* Runs behind the Stata tabs of "Justifying the error rate: alpha by sample
* size, or by error costs". Start from this folder (see _prelude.do).
do "_prelude.do" justified-alpha

* >>> world
powerape dgp d, focal(program, binary p(0.5)) covariates((z, normal))     ///
    baseline(0.30) signal(0.10)
powerape setape d, target(0.10)
* <<<

tic
* >>> default
powerape power d, n(2000) claim(minimum) sesoi(0.05) nsim(400) seed(1)
powerape power d, n(2000) claim(minimum) sesoi(0.05) alpha(0.01) nsim(400) seed(1)
* <<<
keep_r default
toc default

* >>> rule
* alphaN's rule alphaN(n, BF = 3) (Wulff & Taylor, 2024) as an expression in @n
display chi2tail(1, ln(500) + 2*ln(3)), chi2tail(1, ln(2000) + 2*ln(3)), ///
    chi2tail(1, ln(8000) + 2*ln(3))
* <<<

tic
* >>> rule-n
powerape n d, power(0.80) claim(minimum) sesoi(0.05)                     ///
    alpharule(chi2tail(1, ln(@n) + 2*ln(3)))                              ///
    nsim(300) seed(2) noconfirm engine(mata)
* <<<
keep_r rule_n
toc rule_n

tic
* >>> rule-power
powerape power d, n(5000) claim(minimum) sesoi(0.05)                      ///
    alpharule(chi2tail(1, ln(@n) + 2*ln(3))) nsim(400) seed(3)
powerape statement
* <<<
keep_r rule_power
toc rule_power

tic
* >>> costs
powerape power d, n(2000) claim(minimum) sesoi(0.05) nsim(1000) seed(4) name(pw)
powerape powerat pw, alpha(0.005 0.01 0.025 0.05 0.10)
powerape alpha pw, cost(4) prior(1)
* <<<
keep_r costs
toc costs

* >>> costs-plot
local opt = r(alpha)
local lo  = r(range_lo)
local hi  = r(range_hi)
matrix G = r(grid)
svmat G, names(col)
twoway line wcer alpha, xline(`opt') xline(`lo' `hi', lpattern(dash))     ///
    ytitle("Weighted combined error rate") xtitle("Alpha")
* <<<
pc_graph justified-alpha-costs
drop _all

* >>> equal-costs
powerape alpha pw, cost(1) prior(1)
* <<<

* >>> frontier
powerape frontier, cost(4) power(0.80)
powerape frontier, cost(1) power(0.80) claim(minimum)
* <<<

tic
* >>> size
powerape alpha pw, cost(4) prior(1) size nsimsize(500) seed(5)
powerape statement _lastalpha
* <<<
keep_r size
toc size

file close res
file close tim
log close pc
exit, clear
