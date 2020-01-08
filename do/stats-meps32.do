*************************************************
* Constructing MEPS and NHIS Panel Statistics                      *
*************************************************

use "$src/clean/meps-nhis_wide_extract.dta", clear
sum
sum year age byr totinc obese smokev mx mepsexp inscov 
sum year age byr totinc hlth obese smokev mx mepsexp inscov 
tab inscov
tab hlth


