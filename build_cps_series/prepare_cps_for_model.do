
********************
****PROGRAM INFO****
********************

* This program is a slightly amended version of the code written by 
* John Robertson & Ellyn Terry for the Atlanta Fed Wage Growth Tracker.
* Please cite appropriately.

* It produces the CPS time series for the model.

* NOTE: THIS VERSION OF THE PROGRAM IS FOR THE Harmonized Variable and 
* Longitudinally Matched [Atlanta Federal Reserve] (1976-Present) 
* DOWNLOADED FROM CADRE: https://cps.kansascityfed.org/cps*/

**********************************
****SET PATHS AND READ IN DATA****
**********************************

set more off, permanently
display "$S_TIME"

**** set directories ****
do globals.do

***** set weighting option *********
global weights [pweight=weightern82] // 'global weights' not to weight


*************************
****IMPORT CADRE DATA****
*************************
#delimit ;
use 
personid 
recession76 
age76 
wageperhrclean82 
weightern82 
wagegrowthtracker83 
date 
using "${raw}/${cps_dta}", clear;
#delimit cr

keep if date>=mdy(1,1,1982) & date<=mdy(9,1,2023) & age76>=16

/* merge on the WGT group data created by create_MWT_groups_usingcadre.do */
merge 1:1 personid date using "${tmp}WGT_groups"

// generate various date variables
cap gen _year = year(date)
cap gen month = month(date)
cap gen date_monthly = ym(_year,month(date))
cap format date_monthly %tm

xtset personid date_monthly, monthly

// rename wagegrowth variable to simplify
// NB. WGT definition is wagegrowthtracker83 = 100*(log(wageperhrclean82) - log(l12.wageperhrclean82))
rename wagegrowthtracker83 WGT

// Create number of zero wage changes
gen WGT_zer = 0 if WGT!=.
replace WGT_zer = 100 if abs(WGT)<0.5

*******************************************
****Create (un)weighted WGT time series****
*******************************************
// 1st, drop missing WGT observations from dataset except for 85-86 and 95-96 when ALL WGT obs are missing for some months due to Census masking of identifiers (need to keep those missing months for collapsed dataset)
keep if WGT !=. | inlist(_year,1995,1996,1985,1986)

// For unsmoothed version 
gen WGT_raw = WGT

// average wage quartiles
gen WGT_q1 = WGT if wagegroup == "1st"
gen WGT_q2 = WGT if wagegroup == "2nd"
gen WGT_q3 = WGT if wagegroup == "3rd"
gen WGT_q4 = WGT if wagegroup == "4th"

// metro and non-metro 
gen WGT_ym = WGT if msagroup == "MSA"
gen WGT_nm = WGT if msagroup=="NonMSA"

// 3 age groups 
gen WGT_ya = WGT if agegroup=="16-24"
gen WGT_pa = WGT if agegroup=="25-54"
gen WGT_oa = WGT if agegroup=="55+"

// usually ft/usually pt 
gen WGT_ft = WGT if ftptgroup=="Full-time"
gen WGT_pt = WGT if ftptgroup=="Part-time"

// male/female 
gen WGT_ms = WGT if gengroup=="Male"
gen WGT_ws = WGT if gengroup=="Female"

// degree  education 
gen WGT_he = WGT if edgroup3=="Bachelor+" | edgroup3=="Associates"

// 3 ed groups
gen WGT_de = WGT if edgroup3=="Bachelor+"
gen WGT_ae = WGT if edgroup3=="Associates"
gen WGT_le = WGT if edgroup3=="Nodegree"

// occupation groups
gen WGT_lo = WGT if skillgroup=="LowSkill"
gen WGT_mo = WGT if skillgroup=="MidSkill"
gen WGT_ho = WGT if skillgroup=="HigSkill"

// service and goods industries 
gen WGT_si = WGT if secgroup=="Services"
gen WGT_gi = WGT if secgroup=="Goods"

// White and other race 
gen WGT_wr = WGT if racegroup=="White"
gen WGT_or = WGT if racegroup=="Nonwhite"

// job stayer/switcher 
gen WGT_jst = WGT if jstayergroup == "Job Stayer"
gen WGT_jsw = WGT if jstayergroup == "Job Switcher"

gen WGT_switch_lh= WGT if jstayergroup == "Job Switcher" & indgroup=="Leisure & Hospitality"
gen WGT_switch_low= WGT if jstayergroup == "Job Switcher" & occgroup=="Low"
gen WGT_switch_middle= WGT if jstayergroup == "Job Switcher" & occgroup=="Middle"
gen WGT_switch_high= WGT if jstayergroup == "Job Switcher" & occgroup=="High"

// census division
gen WGT_pac = WGT if cdivgroup == "pac"
gen WGT_esc = WGT if cdivgroup == "esc"
gen WGT_wsc = WGT if cdivgroup == "wsc"
gen WGT_mnt = WGT if cdivgroup == "mnt"
gen WGT_nen = WGT if cdivgroup == "nen"
gen WGT_sat = WGT if cdivgroup == "sat"
gen WGT_wnc = WGT if cdivgroup == "wnc"
gen WGT_enc = WGT if cdivgroup == "enc"
gen WGT_mat = WGT if cdivgroup == "mat"

// industries
gen WGT_cmi = WGT if indgroup=="Construction & Mining"
gen WGT_ehi = WGT if indgroup=="Education & Health"
gen WGT_fpi = WGT if indgroup=="Finance and Business Services"
gen WGT_lhi = WGT if indgroup=="Leisure & Hospitality"
gen WGT_mni = WGT if indgroup=="Manufacturing"
gen WGT_pai = WGT if indgroup=="Public Administration"
gen WGT_tti = WGT if indgroup=="Trade & Transportation"

// hourly
gen WGT_yhr = WGT if hrlygroup=="Hourly"
gen WGT_nhr = WGT if hrlygroup=="Non-Hourly"

// occupation by age
gen WGT_ageocc_ya_lo = WGT if ageoccgroup == "16to24_LowSkill"
gen WGT_ageocc_ya_mo = WGT if ageoccgroup == "16to24_MidSkill"
gen WGT_ageocc_ya_ho = WGT if ageoccgroup == "16to24_HigSkill"

gen WGT_ageocc_pa_lo = WGT if ageoccgroup == "25to54_LowSkill"
gen WGT_ageocc_pa_mo = WGT if ageoccgroup == "25to54_MidSkill"
gen WGT_ageocc_pa_ho = WGT if ageoccgroup == "25to54_HigSkill"

gen WGT_ageocc_oa_lo = WGT if ageoccgroup == "55plus_LowSkill"
gen WGT_ageocc_oa_mo = WGT if ageoccgroup == "55plus_MidSkill"
gen WGT_ageocc_oa_ho = WGT if ageoccgroup == "55plus_HigSkill"

// occupation by age by gender
gen WGT_ageoccgend_ya_lo_ms = WGT if ageoccgendgroup == "16to24_LowSkill_Men"
gen WGT_ageoccgend_ya_mo_ms = WGT if ageoccgendgroup == "16to24_MidSkill_Men"
gen WGT_ageoccgend_ya_ho_ms = WGT if ageoccgendgroup == "16to24_HigSkill_Men"

gen WGT_ageoccgend_pa_lo_ms = WGT if ageoccgendgroup == "25to54_LowSkill_Men"
gen WGT_ageoccgend_pa_mo_ms = WGT if ageoccgendgroup == "25to54_MidSkill_Men"
gen WGT_ageoccgend_pa_ho_ms = WGT if ageoccgendgroup == "25to54_HigSkill_Men"

gen WGT_ageoccgend_oa_lo_ms = WGT if ageoccgendgroup == "55plus_LowSkill_Men"
gen WGT_ageoccgend_oa_mo_ms = WGT if ageoccgendgroup == "55plus_MidSkill_Men"
gen WGT_ageoccgend_oa_ho_ms = WGT if ageoccgendgroup == "55plus_HigSkill_Men"

gen WGT_ageoccgend_ya_lo_ws = WGT if ageoccgendgroup == "16to24_LowSkill_Wom"
gen WGT_ageoccgend_ya_mo_ws = WGT if ageoccgendgroup == "16to24_MidSkill_Wom"
gen WGT_ageoccgend_ya_ho_ws = WGT if ageoccgendgroup == "16to24_HigSkill_Wom"

gen WGT_ageoccgend_pa_lo_ws = WGT if ageoccgendgroup == "25to54_LowSkill_Wom"
gen WGT_ageoccgend_pa_mo_ws = WGT if ageoccgendgroup == "25to54_MidSkill_Wom"
gen WGT_ageoccgend_pa_ho_ws = WGT if ageoccgendgroup == "25to54_HigSkill_Wom"

gen WGT_ageoccgend_oa_lo_ws = WGT if ageoccgendgroup == "55plus_LowSkill_Wom"
gen WGT_ageoccgend_oa_mo_ws = WGT if ageoccgendgroup == "55plus_MidSkill_Wom"
gen WGT_ageoccgend_oa_ho_ws = WGT if ageoccgendgroup == "55plus_HigSkill_Wom"

keep personid _year month date_monthly WGT* *group* recession76 weightern82
save "${tmp}WGT_unweighted.dta", replace  /* these are the unweighted individual level WGT observations for the various cuts */

/* collapse WGT dataset into a time series */
use "${tmp}WGT_unweighted.dta", clear

#delimit ;

collapse
(median)
WGT WGT_raw
WGT_ya WGT_pa WGT_oa
WGT_ft WGT_pt
WGT_ms WGT_ws
WGT_he WGT_de WGT_ae WGT_le
WGT_si WGT_gi
WGT_wr WGT_or WGT_ho WGT_lo WGT_mo
WGT_jst WGT_jsw
WGT_cmi WGT_ehi WGT_fpi WGT_lhi WGT_mni WGT_pai WGT_tti
WGT_pac WGT_esc WGT_wsc WGT_mnt WGT_nen WGT_sat WGT_wnc WGT_enc WGT_mat
WGT_ym WGT_nm
WGT_q1 WGT_q2 WGT_q3 WGT_q4
WGT_yhr WGT_nhr
WGT_ageocc_ya_lo WGT_ageocc_ya_mo WGT_ageocc_ya_ho 
WGT_ageocc_pa_lo WGT_ageocc_pa_mo WGT_ageocc_pa_ho 
WGT_ageocc_oa_lo WGT_ageocc_oa_mo WGT_ageocc_oa_ho 
WGT_ageoccgend_ya_lo_ms WGT_ageoccgend_ya_mo_ms WGT_ageoccgend_ya_ho_ms
WGT_ageoccgend_pa_lo_ms WGT_ageoccgend_pa_mo_ms WGT_ageoccgend_pa_ho_ms
WGT_ageoccgend_oa_lo_ms WGT_ageoccgend_oa_mo_ms WGT_ageoccgend_oa_ho_ms
WGT_ageoccgend_ya_lo_ws WGT_ageoccgend_ya_mo_ws WGT_ageoccgend_ya_ho_ws
WGT_ageoccgend_pa_lo_ws WGT_ageoccgend_pa_mo_ws WGT_ageoccgend_pa_ho_ws
WGT_ageoccgend_oa_lo_ws WGT_ageoccgend_oa_mo_ws WGT_ageoccgend_oa_ho_ws     
(mean)
WGT_avg=WGT ZERO=WGT_zer REC=recession76
(p25) WGT_p25=WGT
(p75) WGT_p75=WGT 
(count)
WGT_n=WGT
WGT_ya_n=WGT_ya WGT_pa_n=WGT_pa WGT_oa_n=WGT_oa
WGT_ft_n=WGT_ft WGT_pt_n=WGT_pt
WGT_ms_n=WGT_ms WGT_ws_n=WGT_ws
WGT_he_n=WGT_he WGT_de_n=WGT_de WGT_ae_n=WGT_ae WGT_le_n=WGT_le
WGT_si_n=WGT_si WGT_gi_n=WGT_gi 
WGT_wr_n=WGT_wr WGT_or_n=WGT_or WGT_ho_n=WGT_ho WGT_lo_n=WGT_lo WGT_mo_n=WGT_mo
WGT_jst_n=WGT_jst WGT_jsw_n=WGT_jsw
WGT_cmi_n=WGT_cmi WGT_ehi_n=WGT_ehi WGT_fpi_n=WGT_fpi WGT_lhi_n=WGT_lhi WGT_mni_n=WGT_mni WGT_pai_n=WGT_pai WGT_tti_n=WGT_tti
WGT_pac_n=WGT_pac WGT_esc_n=WGT_esc WGT_wsc_n=WGT_wsc WGT_mnt_n=WGT_mnt WGT_nen_n=WGT_nen WGT_sat_n=WGT_sat WGT_wnc_n=WGT_wnc WGT_enc_n=WGT_enc WGT_mat_n=WGT_mat
WGT_ym_n=WGT_ym WGT_nm_n=WGT_nm
WGT_q1_n=WGT_q1 WGT_q2_n=WGT_q2 WGT_q3_n=WGT_q3 WGT_q4_n=WGT_q4
WGT_yhr_n=WGT_yhr WGT_nhr_n=WGT_nhr
WGT_ageocc_ya_lo_n=WGT_ageocc_ya_lo WGT_ageocc_ya_mo_n=WGT_ageocc_ya_mo WGT_ageocc_ya_ho_n=WGT_ageocc_ya_ho 
WGT_ageocc_pa_lo_n=WGT_ageocc_pa_lo WGT_ageocc_pa_mo_n=WGT_ageocc_pa_mo WGT_ageocc_pa_ho_n=WGT_ageocc_pa_ho 
WGT_ageocc_oa_lo_n=WGT_ageocc_oa_lo WGT_ageocc_oa_mo_n=WGT_ageocc_oa_mo WGT_ageocc_oa_ho_n=WGT_ageocc_oa_ho
WGT_ageoccgend_ya_lo_ms_n=WGT_ageoccgend_ya_lo_ms WGT_ageoccgend_ya_mo_ms_n=WGT_ageoccgend_ya_mo_ms WGT_ageoccgend_ya_ho_ms_n=WGT_ageoccgend_ya_ho_ms
WGT_ageoccgend_pa_lo_ms_n=WGT_ageoccgend_pa_lo_ms WGT_ageoccgend_pa_mo_ms_n=WGT_ageoccgend_pa_mo_ms WGT_ageoccgend_pa_ho_ms_n=WGT_ageoccgend_pa_ho_ms
WGT_ageoccgend_oa_lo_ms_n=WGT_ageoccgend_oa_lo_ms WGT_ageoccgend_oa_mo_ms_n=WGT_ageoccgend_oa_mo_ms WGT_ageoccgend_oa_ho_ms_n=WGT_ageoccgend_oa_ho_ms
WGT_ageoccgend_ya_lo_ws_n=WGT_ageoccgend_ya_lo_ws WGT_ageoccgend_ya_mo_ws_n=WGT_ageoccgend_ya_mo_ws WGT_ageoccgend_ya_ho_ws_n=WGT_ageoccgend_ya_ho_ws
WGT_ageoccgend_pa_lo_ws_n=WGT_ageoccgend_pa_lo_ws WGT_ageoccgend_pa_mo_ws_n=WGT_ageoccgend_pa_mo_ws WGT_ageoccgend_pa_ho_ws_n=WGT_ageoccgend_pa_ho_ws
WGT_ageoccgend_oa_lo_ws_n=WGT_ageoccgend_oa_lo_ws WGT_ageoccgend_oa_mo_ws_n=WGT_ageoccgend_oa_mo_ws WGT_ageoccgend_oa_ho_ws_n=WGT_ageoccgend_oa_ho_ws      
${weights}, 
by(date_monthly _year month)
;

#delimit cr

/*Create date variable*/
gen date = mdy(month,1,_year)
sort date
format date %tdDD-Mon-CCYY
save "${tmp}WGT_collapsed.dta", replace  


*********************************************
****Prepare time series to estimate TWIn ****
*********************************************
use "${tmp}WGT_collapsed.dta", clear // raw series for model 
keep if date>=mdy(1,1,1997) // when WGT starts

// job switchers vs. stayers
local dirpath "${processed}switchers_vs_stayers/"
cap mkdir `dirpath'
export delimited date WGT_jst WGT_jsw using `dirpath'wageinflation.csv, replace
export delimited date WGT_jst_n WGT_jsw_n using `dirpath'weights.csv, replace

// industries
local dirpath "${processed}industries/"
cap mkdir `dirpath'
export delimited date WGT_cmi WGT_ehi WGT_fpi WGT_lhi WGT_mni WGT_pai WGT_tti using `dirpath'wageinflation.csv, replace
export delimited date WGT_cmi_n WGT_ehi_n WGT_fpi_n WGT_lhi_n WGT_mni_n WGT_pai_n WGT_tti_n using `dirpath'weights.csv, replace

// goods vs. services
local dirpath "${processed}goods_vs_services/"
cap mkdir `dirpath'
export delimited date WGT_gi WGT_si using `dirpath'wageinflation.csv, replace
export delimited date WGT_gi_n WGT_si_n using `dirpath'weights.csv, replace

// male vs female
local dirpath "${processed}gender/"
cap mkdir `dirpath'
export delimited date WGT_ms WGT_ws using `dirpath'wageinflation.csv, replace
export delimited date WGT_ms_n WGT_ws_n using `dirpath'weights.csv, replace

// white vs nonwhite
local dirpath "${processed}race/"
cap mkdir `dirpath'
export delimited date WGT_wr WGT_or using `dirpath'wageinflation.csv, replace
export delimited date WGT_wr_n WGT_or_n using `dirpath'weights.csv, replace

// education
local dirpath "${processed}education/"
cap mkdir `dirpath'
export delimited date WGT_le WGT_ae WGT_de using `dirpath'wageinflation.csv, replace
export delimited date WGT_le_n WGT_ae_n WGT_de_n using `dirpath'weights.csv, replace

// occupation group
local dirpath "${processed}occupation/"
cap mkdir `dirpath'
export delimited date WGT_lo WGT_mo WGT_ho using `dirpath'wageinflation.csv, replace
export delimited date WGT_lo_n WGT_mo_n WGT_ho_n using `dirpath'weights.csv, replace

// age groups
local dirpath "${processed}age/"
cap mkdir `dirpath'
export delimited date WGT_ya WGT_pa WGT_oa using `dirpath'wageinflation.csv, replace
export delimited date WGT_ya_n WGT_pa_n WGT_oa_n using `dirpath'weights.csv, replace

// msa vs nonmsa
local dirpath "${processed}msa_vs_nonmsa/"
cap mkdir `dirpath'
export delimited date WGT_ym WGT_nm using `dirpath'wageinflation.csv, replace
export delimited date WGT_ym_n WGT_nm_n using `dirpath'weights.csv, replace

// census division groups
local dirpath "${processed}region/"
cap mkdir `dirpath'
export delimited date WGT_pac WGT_esc WGT_wsc WGT_mnt WGT_nen WGT_sat WGT_wnc WGT_enc WGT_mat using `dirpath'wageinflation.csv, replace
export delimited date WGT_pac_n WGT_esc_n WGT_wsc_n WGT_mnt_n WGT_nen_n WGT_sat_n WGT_wnc_n WGT_enc_n WGT_mat_n using `dirpath'weights.csv, replace

// wage quartile
local dirpath "${processed}wage_quartile/"
cap mkdir `dirpath'
export delimited date WGT_q1 WGT_q2 WGT_q3 WGT_q4 using `dirpath'wageinflation.csv, replace
export delimited date WGT_q1_n WGT_q2_n WGT_q3_n WGT_q4_n using `dirpath'weights.csv, replace

// age by occupation
local dirpath "${processed}age_occupation/"
cap mkdir `dirpath'
#delimit ;
export delimited date
	WGT_ageocc_ya_lo WGT_ageocc_ya_mo WGT_ageocc_ya_ho 
	WGT_ageocc_pa_lo WGT_ageocc_pa_mo WGT_ageocc_pa_ho 
	WGT_ageocc_oa_lo WGT_ageocc_oa_mo WGT_ageocc_oa_ho 
	using `dirpath'wageinflation.csv, replace;
export delimited date 
	WGT_ageocc_ya_lo_n WGT_ageocc_ya_mo_n WGT_ageocc_ya_ho_n 
	WGT_ageocc_pa_lo_n WGT_ageocc_pa_mo_n WGT_ageocc_pa_ho_n 
	WGT_ageocc_oa_lo_n WGT_ageocc_oa_mo_n WGT_ageocc_oa_ho_n 
	using `dirpath'weights.csv, replace;
#delimit cr

// age by occupation by gender
local dirpath "${processed}age_occupation_gender/"
cap mkdir `dirpath'
#delimit ;
export delimited date 
	WGT_ageoccgend_ya_lo_ms WGT_ageoccgend_ya_mo_ms WGT_ageoccgend_ya_ho_ms
	WGT_ageoccgend_pa_lo_ms WGT_ageoccgend_pa_mo_ms WGT_ageoccgend_pa_ho_ms
	WGT_ageoccgend_oa_lo_ms WGT_ageoccgend_oa_mo_ms WGT_ageoccgend_oa_ho_ms
	WGT_ageoccgend_ya_lo_ws WGT_ageoccgend_ya_mo_ws WGT_ageoccgend_ya_ho_ws
	WGT_ageoccgend_pa_lo_ws WGT_ageoccgend_pa_mo_ws WGT_ageoccgend_pa_ho_ws
	WGT_ageoccgend_oa_lo_ws WGT_ageoccgend_oa_mo_ws WGT_ageoccgend_oa_ho_ws     
	using `dirpath'wageinflation.csv, replace;
export delimited date 
	WGT_ageoccgend_ya_lo_ms_n WGT_ageoccgend_ya_mo_ms_n WGT_ageoccgend_ya_ho_ms_n
	WGT_ageoccgend_pa_lo_ms_n WGT_ageoccgend_pa_mo_ms_n WGT_ageoccgend_pa_ho_ms_n
	WGT_ageoccgend_oa_lo_ms_n WGT_ageoccgend_oa_mo_ms_n WGT_ageoccgend_oa_ho_ms_n
	WGT_ageoccgend_ya_lo_ws_n WGT_ageoccgend_ya_mo_ws_n WGT_ageoccgend_ya_ho_ws_n
	WGT_ageoccgend_pa_lo_ws_n WGT_ageoccgend_pa_mo_ws_n WGT_ageoccgend_pa_ho_ws_n
	WGT_ageoccgend_oa_lo_ws_n WGT_ageoccgend_oa_mo_ws_n WGT_ageoccgend_oa_ho_ws_n
	using `dirpath'weights.csv, replace;
#delimit cr

* Use mean instead of median for collapse *************
use "${tmp}WGT_unweighted.dta", clear

gen date = mdy(month,1,_year)
keep if date>=mdy(1,1,1997) // when WGT starts

#delimit ;
collapse
(mean) WGT_cmi WGT_ehi WGT_fpi WGT_lhi WGT_mni WGT_pai WGT_tti
(count) WGT_cmi_n=WGT_cmi WGT_ehi_n=WGT_ehi WGT_fpi_n=WGT_fpi WGT_lhi_n=WGT_lhi WGT_mni_n=WGT_mni WGT_pai_n=WGT_pai WGT_tti_n=WGT_tti
${weights}, 
by(date_monthly _year month);
#delimit cr

format date %tdDD-Mon-CCYY

// industries
local dirpath "${processed}industries_average/"
cap mkdir `dirpath'
export delimited date WGT_cmi WGT_ehi WGT_fpi WGT_lhi WGT_mni WGT_pai WGT_tti using `dirpath'wageinflation.csv, replace
export delimited date WGT_cmi_n WGT_ehi_n WGT_fpi_n WGT_lhi_n WGT_mni_n WGT_pai_n WGT_tti_n using `dirpath'weights.csv, replace


* Unweighted cuts *************
use "${tmp}WGT_unweighted.dta", clear

gen date = mdy(month,1,_year)
keep if date>=mdy(1,1,1997) // when WGT starts

#delimit ;
collapse
(p50) WGT_cmi WGT_ehi WGT_fpi WGT_lhi WGT_mni WGT_pai WGT_tti
(count) WGT_cmi_n=WGT_cmi WGT_ehi_n=WGT_ehi WGT_fpi_n=WGT_fpi WGT_lhi_n=WGT_lhi WGT_mni_n=WGT_mni WGT_pai_n=WGT_pai WGT_tti_n=WGT_tti, 
by(date_monthly _year month);
#delimit cr

format date %tdDD-Mon-CCYY

// industries
local dirpath "${processed}industries_unweighted/"
cap mkdir `dirpath'
export delimited date WGT_cmi WGT_ehi WGT_fpi WGT_lhi WGT_mni WGT_pai WGT_tti using `dirpath'wageinflation.csv, replace
export delimited date WGT_cmi_n WGT_ehi_n WGT_fpi_n WGT_lhi_n WGT_mni_n WGT_pai_n WGT_tti_n using `dirpath'weights.csv, replace

