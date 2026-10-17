* Convert the help file of the Stata command to plain text for the website
* article "Stata help file". Start from this folder; POWERAPE_STATA_DIR must
* point to the folder that holds powerape.sthlp.
clear all
set more off
local adodir : environment POWERAPE_STATA_DIR
if "`adodir'" == "" {
    display as error "set POWERAPE_STATA_DIR to the folder with powerape.sthlp"
    exit 198
}
capture mkdir "output"
translator set smcl2txt linesize 92
quietly translate "`adodir'/powerape.sthlp" "output/powerape-help.txt", ///
    translator(smcl2txt) replace
exit, clear
