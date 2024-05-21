* Homework 2
* Author: Jiahui Shui
* Problem 1 (a), (b), (c) and (d) modified code from John (TA session)

* Import data
clear all

global fgcolor "15 60 15"
// ssc install freduse, replace

* Load Data
local tsvar "FEDFUNDS UNRATE GDPDEF USRECM"

foreach v of local tsvar {
	import delimited using "data/`v'.csv", clear case(preserve)
	rename DATE date
	tempfile `v'_dta
	save ``v'_dta', replace
	}

use `FEDFUNDS_dta', clear
keep date
foreach v of local tsvar {
	joinby date using ``v'_dta', unm(b)
	drop _merge
	}
	
* Clean Data
gen daten = date(date, "YMD")
format daten %td
drop if yofd(daten) < 1960  | yofd(daten) > 2023 // data is per quarter in 1947 but per month after 
gen INFL = 100*(GDPDEF - GDPDEF[_n-12])/GDPDEF[_n-12] //year to year inflation
    la var INFL "Inflation Rate"
    la var FEDFUNDS "Federal Funds Rate"
    la var UNRATE "Unemployment Rate"
	la var GDPDEF   "GDP Deflator"
    la var daten Date // re=label date 
tsset daten
local tsvar "FEDFUNDS UNRATE INFL"

* 1(a) Plot the data
tsline FEDFUND UNRATE, yaxis(1) || tsline GDPDEF, yaxis(2) graphregion(fcolor(white)) xlabel(,labsize(small))
graph export "1a_plot.pdf", replace

* 1(b) Aggregate Quarterly Data
gen dateq = qofd(daten)
collapse (mean) `tsvar' (min) daten, by(dateq)
tsset dateq, quarterly
keep if (yofd(daten) >= 1960) & (yofd(daten) <= 2007)
var INFL UNRATE FEDFUNDS, lags(1/4)

* 1(c) See PDF
* 1(d) VAR Estimation
matrix A = (.,0,0 \ .,.,0 \ .,.,.)
matrix B = (1,0,0 \ 0,1,0 \ 0,0,1)
svar INFL UNRATE FEDFUNDS, lags(1/4) aeq(A) beq(B)
	irf create mysirf, set(mysirfs) step(20) replace
	
* Plot and save IRFs
irf graph sirf, irf(mysirf) impulse(INFL UNRATE FEDFUNDS) response(INFL UNRATE FEDFUNDS)  ///
byopts(graphregion(fcolor(white)) yrescale /// 
xrescale note("") legend(pos(3) )) legend(stack col(1) /// 
order(1 "95% CI" 2 "SIRF") symx(*.5) size(vsmall))  xtitle("Quarters after shock")	
graph export "1d_svar_irf.pdf", replace

* 1(f)
predict resid_monetary, residuals equation(FEDFUNDS)
tsline resid_monetary, graphregion(fcolor(white)) 
graph export "1f_money_shocks.pdf", replace

** Problem 2
* 2(a)
rename dateq date
merge 1:1 date using "Data/RR_monetary_shock_quarterly.dta", assert (master match) nogen

replace resid		= 0 if yofd(daten) < 1969
replace resid_romer = 0 if yofd(daten) < 1969
replace resid_full	= 0 if yofd(daten) < 1969

* 2(b)
tsset date
var INFL UNRATE FEDFUNDS, lags(1/8)  exog(L(0/12).resid_full)

irf create rirf, step(20) replace
irf graph dm, impulse(resid_full) irf(rirf) graphregion(fcolor(white))
graph export "2b_var_rirf.pdf", replace

* 2(c)
matrix AA = (.,0,0,0 \ .,.,0,0 \ .,.,.,0 \ .,.,.,.)
matrix BB = (1,0,0,0 \ 0,1,0,0 \ 0,0,1,0 \ 0,0,0,1)

svar resid_full INFL UNRATE FEDFUNDS, lags(1/4) aeq(AA) beq(BB)
irf create svarrr, set(mysvarrr) step(20) replace

irf graph sirf, irf(svarrr) impulse(resid_full INFL UNRATE FEDFUNDS) response(resid_full INFL UNRATE FEDFUNDS) ///
byopts(graphregion(fcolor(white)) yrescale /// 
xrescale note("") legend(pos(3) )) legend(stack col(1) /// 
order(1 "95% CI" 2 "SIRF") symx(*.5) size(vsmall))  xtitle("Quarters after shock")	
graph export "2c_svar_RR_irf.pdf", replace

* 2(g)
predict resid_rr_monetary, residuals equation(FEDFUNDS)

label var resid_rr_monetary "Romer Shocks" ///

tsline resid_rr_monetary || tsline resid_monetary, graphregion(fcolor(white)) 
graph export "2g_money_shocks.pdf", replace
