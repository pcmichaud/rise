*************************************************
* Constructing PSID Panel Statistics                     *
*************************************************

use "$src/clean/psid_extract.dta", clear

sum
keep if rage>=25&rage<=84
sum
	
*birth date
tab rbyr_c

*Health 3 and 5
sum rshlt rshlt3
tab rshlt3
*Education
sum rskill
tab rskill
*Participation
sum rwork
tab rwork

*Income
gen riearn2 =riearn

replace riearn2 = . if riearn>=250000
replace otherinc = . if otherinc>=250000

*Wealth
replace ass = 0 if ass<0

table year, c(mean ass)

sum rage rskill rwork riearn2 otherinc ass 


	
gen logass = log(1+ass)
sum logass
	

