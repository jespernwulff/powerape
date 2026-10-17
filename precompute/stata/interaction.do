* Runs behind the Stata tabs of "Power for average interaction effects (AIEs)".
* Start from this folder (see _prelude.do).
do "_prelude.do" interaction

* >>> world
powerape dgp d, focal(treat, binary p(0.5)) moderator(female, binary p(0.55)) ///
    covariates((age, normal mean(45) sd(12))) baseline(0.30) signal(0.10)
* <<<

* >>> pin
powerape setaie d, target(0.08) mainfocal(0.10) mainmoderator(0.03)
* <<<

tic
* >>> power
powerape power d, n(3000) claim(detect) nsim(500) seed(1)
* <<<
keep_r power
toc power

* >>> main-world
powerape dgp d_main, focal(treat, binary p(0.5))                          ///
    covariates((age, normal mean(45) sd(12))) baseline(0.30) signal(0.10)
powerape setape d_main, target(0.08)
* <<<

tic
* >>> main-power
powerape power d_main, n(3000) claim(detect) nsim(500) seed(2)
* <<<
keep_r main_power
toc main_power

tic
* >>> n
powerape n d, power(0.80) claim(detect) nsim(500) seed(3) engine(mata)
* <<<
keep_r n
toc n

tic
* >>> closed25
powerape power d, n(3000) claim(detect) nsim(25) seed(1)
* <<<
keep_r closed25
toc closed25

tic
* >>> ginteff
powerape power d, n(3000) claim(detect) nsim(25) seed(1) effects(margins)
* <<<
keep_r ginteff
toc ginteff

file close res
file close tim
log close pc
exit, clear
