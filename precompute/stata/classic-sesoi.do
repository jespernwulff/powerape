* Runs behind the Stata tabs of "The classic SESOI power analysis, and what it
* cannot tell you". Start from this folder (see _prelude.do).
do "_prelude.do" classic-sesoi

* >>> world
powerape dgp d, focal(treat, binary p(0.5)) baseline(0.30)
powerape setape d, target(0.05)
* <<<

tic
* >>> recipe
powerape power d, n(5000) claim(detect) sesoi(0.05) nsim(600) seed(1)
* <<<
keep_r recipe
toc recipe

* >>> wall
capture noisily powerape power d, n(5000) claim(minimum) sesoi(0.05)
* <<<

tic
* >>> n-detect
* detection power for truth = SESOI = 0.05
powerape n d, power(0.80) claim(detect) nsim(500) seed(2) engine(mata)
* <<<
keep_r n_detect
toc n_detect

* >>> pin8
* minimum-effect power for truth = 0.08 against SESOI = 0.05
powerape setape d, target(0.08)
* <<<

tic
* >>> n-minimum
powerape n d, power(0.80) claim(minimum) sesoi(0.05) nsim(500) seed(3) engine(mata)
* <<<
keep_r n_minimum
toc n_minimum

* >>> pin0
powerape setape d, target(0)
* <<<

tic
* >>> size
powerape power d, n(2000) claim(detect) nsim(600) seed(4)
* <<<
keep_r size
toc size

file close res
file close tim
log close pc
exit, clear
