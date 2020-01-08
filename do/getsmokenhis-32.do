
* getting smoking histories 1997 to 2008

forvalues j = 1997/2008 {
	use $src/nhis/samadult-`j'.dta, clear
	gen nhis_yr = `j'
	capture gen fpx = px
	capture drop publicid
	di "`j'"
	tab fmx
	tab fpx
	egen publicid = concat(nhis_yr hhx fmx fpx)
	keep publicid smkstat2 nhis_yr hhx fmx fpx
	save $src/temp/nhis-smoke-`j'.dta, replace
}
use $src/temp/nhis-smoke-1997.dta, clear
forvalues j = 1998/2008 {
	append using $src/temp/nhis-smoke-`j'.dta
}
recode smkstat2 (1/2 5 = 1) (3=2) (4=3) (else=.), gen(smoke3)
recode smoke3 (1/2 = 1) (3=0) (else=.), gen(smokev)
recode smoke3 (1 = 1) (2 3=0) (else=.), gen(smoken)
tab nhis_yr smokev, row nofreq

keep publicid nhis_yr smoke3 smokev smoken hhx fmx fpx

save $src/temp/nhis-smoke.dta, replace
forvalues j = 1997/2008 {
	erase $src/temp/nhis-smoke-`j'.dta
}
