* Download CPS cadre data directly from Kansas City Fed

set more off, permanently
display "$S_TIME"

do globals.do

* Import
copy "${download_link}" "${raw}`gzname'", replace
// !curl -L "${download_link}" > "${raw}/${cps_gzip}"
!gunzip -cf "${raw}/${cps_gzip}" > "${raw}/${cps_dta}" 


* Check version of data 
desc using "${raw}/${cps_dta}", short



