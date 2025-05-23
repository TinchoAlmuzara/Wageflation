
* set folders 
global raw "data/raw/"
global tmp "data/tmp/"
global processed "data/twin_inputs/"

!mkdir -p "${raw}"
!mkdir -p "${tmp}"
!mkdir -p "${processed}"

* CPS data name 
global cps_dta "CPS_harmonized_variable_longitudinally_matched_age16plus.dta"
global cps_gzip "${cps_dta}.gz"
global download_link "https://cps.kansascityfed.org/data/${cps_gzip}"


