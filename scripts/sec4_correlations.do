/**********************************************************
* Purpose: producing Table 2 in chunks
*	- bottom portion (quarterly wage growth measures)
*	- first 4 rows, Col 3 (quarterly average of monthly wage growth measures v. UR gap)
*	- remainder of top portion (monthly wage growth measures)
***********************************************************/

clear all

global data ./data
global output ./tables

**************************
* QUARTERLY CORRELATIONS *
**************************

* Import TWIn estimates
import delim ${data}/twin_series.csv, clear
gen qdate = qofd(date(date, "DMY"))
collapse (mean) twin common, by(qdate)

tempfile twinq 
save `twinq', replace

* Import monthly measures
import delim ${data}/monthly_series.csv, clear
gen mdate = mofd(date(date, "DMY"))
tsset mdate
gen vacancies_rate = jolts_openings/(unemployment + employment)
gen dln_pce = (pce - L12.pce)/L12.pce
gen dln_pce_exhousing = (pce_exhousing - L12.pce_exhousing)/L12.pce_exhousing
gen dln_ahe = (ahe - L12.ahe)/L12.ahe
gen qdate = qofd(date(date, "DMY"))
collapse (mean) urate vacancies_rate dln_* atlwgt, by(qdate)
tsset qdate

tempfile monthly
save `monthly', replace

* Import quarterly measures
import delim ${data}/quarterly_series.csv, clear
gen qdate = qofd(date(date, "DMY"))
tsset qdate
format qdate %tq

merge 1:1 qdate using `twinq', nogen
merge 1:1 qdate using `monthly', nogen
gen dln_eci = (eci - L4.eci)/L4.eci
gen gap_urate = urate - ustar

* Label variables 
label var twin "TWIn"
label var common "TWIn (Common)"
label var vacancies_rate "V/LF"
label var urate  "Unemployment Rate"
label var dln_pce "Core PCE"
label var dln_pce_exhousing "PCE Services Ex Housing"
label var dln_eci "ECI"
label var dln_ahe "AHE Prod. 12m"
label var atlwgt "Atlanta Fed Wage Tracker"
label var gap_urate "Unemployment Rate Gap"

* First differences
qui ds qdate, not

foreach var in `r(varlist)' {
	gen `var'_diff = D.`var'
}

foreach var of varlist *_diff {
	local lev = subinstr("`var'", "_diff", "", .)
	local lab : var label `lev'
	label var `var' "`lab'"
}

* Pairwise correlation table
order qdate vacancies_rate_diff urate_diff gap_urate_diff dln_pce_diff dln_pce_exhousing_diff twin_diff common_diff dln_eci_diff atlwgt_diff dln_ahe_diff 
keep qdate vacancies_rate_diff urate_diff gap_urate_diff dln_pce_diff dln_pce_exhousing_diff twin_diff common_diff dln_eci_diff atlwgt_diff dln_ahe_diff 

eststo clear
estpost corr *_diff, m nohalf 
esttab using ${output}/Tab2_quarterly.tex, b(%9.3f) unstack not nonum noobs star(* 0.10 ** 0.05 *** 0.01) keep(twin_diff common_diff dln_eci_diff atlwgt_diff dln_ahe_diff) addnotes(Note: Sample period: 1997Q1-2023Q2, except for tightness ratio starting in 2000Q4. All correlations are for variables in first differences.) compress label replace tex

************************
* MONTHLY CORRELATIONS *
************************

* Import TWIn estimates
import delim ${data}/twin_series.csv, clear
ren date _date
gen date = date(_date, "DMY")
drop _date
tempfile twin
save `twin', replace

* Import monthly series
import delimited ${data}/monthly_series.csv, clear
ren date _date
gen date = date(_date, "DMY")
merge 1:1 date using `twin', nogen
gen mdate = mofd(date)
tsset mdate
gen vacancies_rate = jolts_openings/(unemployment + employment)
gen dln_pce = (pce - L12.pce)/L12.pce
gen dln_pce_exhousing = (pce_exhousing - L12.pce_exhousing)/L12.pce_exhousing
gen dln_ahe = (ahe - L12.ahe)/L12.ahe

* Label variables 
label var twin "TWIn"
label var common "TWIn (Common)"
label var vacancies_rate "V/LF"
label var urate  "Unemployment Rate"
label var dln_pce "Core PCE"
label var dln_pce_exhousing "PCE Services Ex Housing"
label var dln_ahe "AHE Prod. 12m"
label var atlwgt "Atlanta Fed Wage Tracker"

qui ds mdate, not

foreach var in `r(varlist)' {
	gen `var'_diff = D.`var'
}

foreach var of varlist *_diff {
	local lev = subinstr("`var'", "_diff", "", .)
	local lab : var label `lev'
	label var `var' "`lab'"
}

* Pairwise correlation table
order mdate vacancies_rate_diff urate_diff dln_pce_diff dln_pce_exhousing_diff twin_diff common_diff dln_ahe_diff atlwgt_diff  
keep mdate vacancies_rate_diff urate_diff dln_pce_diff dln_pce_exhousing_diff twin_diff common_diff dln_ahe_diff atlwgt_diff

eststo clear
estpost corr *_diff, m nohalf 
esttab using ${output}/Tab2_monthly.tex, b(%9.3f) unstack not nonum noobs star(* 0.10 ** 0.05 *** 0.01) keep(twin_diff common_diff dln_ahe_diff atlwgt_diff) addnotes(Note: Sample period: 1997m1-2023m9, except for tightness ratio starting in 2000m12. All correlations are for variables in first differences.) compress label replace tex

