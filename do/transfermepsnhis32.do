* file to transfer to Stata all MEPS-NHIS link files

*global mepsdir "/Users/uqam/Documents/data/MEPS/" 

global mepsdir  "$src/MEPS_NHIS"
global mepsout "$src/temp"

* 1996 (panel 1)
infix str dupersid 1-8 str his_id 9-20 using "$mepsdir/MEPS_1996/his_mepx.dat", clear
gen year = 1996
gen nhis_yr = 1995
gen nhyr = 95
gen hhx = his_id
gen str2 fmx = "01"
gen str2 px = "01"
egen publicid = concat(nhyr his_id)
drop nhyr
gen panel = 1
save "$mepsout/his_meps_1996.dta", replace

* 1997 (panel 2)
infix str dupersid 1-8 str hhx 9-18 str px 19-20 linkflag 21-21 panel97 22-22 using "$mepsdir/MEPS_1997/nhmep97x.dat", clear
gen year = 1997
keep if panel97==2
keep if linkflag==1
gen nhis_yr = 1996
gen nhyr = 96 if nhis_yr==1996
gen str2 fmx = "01"
egen publicid = concat(nhyr hhx px)
gen panel = 2
save "$mepsout/his_meps_1997.dta", replace

* 1998 (panel 3)
infix str dupersid 1-8 str hhx 9-18 str px 19-20 linkflag 21-21 panel98 22-22 using "$mepsdir/MEPS_1998/NHMEP98X.DAT", clear
gen year = 1998
keep if linkflag==1
keep if panel98==3
gen nhis_yr = 1997
gen str2 fmx = "01"
egen publicid = concat(nhis_yr hhx fmx px)
gen panel = 3
save "$mepsout/his_meps_1998.dta", replace

* 1999 (panel 4)
infix str dupersid 1-8 str hhx 9-14 str px 15-16 linkflag 17-17 panel99 18-18 using "$mepsdir/MEPS_1999/NHMEP99X.DAT", clear
gen year = 1999
keep if linkflag==1
keep if panel99==4
gen nhis_yr = 1998
gen str2 fmx = "01"
egen publicid = concat(nhis_yr hhx fmx px)
gen panel = 4
save "$mepsout/his_meps_1999.dta", replace

* 2000 (panel 5)
infix str dupersid 1-8 str hhx 9-14 str px 15-16 linkflag 17-17 panel00 18-18 using "$mepsdir/MEPS_2000/NHMEP00X.DAT", clear
gen year = 2000
keep if linkflag==1
keep if panel00==5
gen nhis_yr = 1999
gen str2 fmx = "01"
egen publicid = concat(nhis_yr hhx fmx px)
gen panel = 5
save "$mepsout/his_meps_2000.dta", replace

* 2001 (panel 6)
infix str dupersid 1-8 str hhx 9-14 str px 15-16 linkflag 17-17 panel01 18-18 srvy_yr 19-22 using "$mepsdir/MEPS_2001/NHMEP01X.DAT", clear
gen year = 2001
tab linkflag
keep if linkflag==1
keep if panel01==6
gen nhis_yr = srvy_yr
gen str2 fmx = "01"
egen publicid = concat(nhis_yr hhx fmx px)
gen panel = 6
save "$mepsout/his_meps_2001.dta", replace

* 2002 (panel 7)
infix str dupersid 1-8 str hhx 9-14 str px 15-16 linkflag 17-17 panel02 18-18 srvy_yr 19-22 using "$mepsdir/MEPS_2002/NHMEP02X.DAT", clear
gen year = 2002
keep if linkflag==1
keep if panel02==7
gen nhis_yr = srvy_yr
gen str2 fmx = "01"
egen publicid = concat(nhis_yr hhx fmx px)
gen panel = 7
save "$mepsout/his_meps_2002.dta", replace

* 2003 (panel 8)
infix str dupersid 1-8 str hhx 9-14 str px 15-16 linkflag 17-17 panel03 18-18 srvy_yr 19-22 using "$mepsdir/MEPS_2003/NHMEP03X.DAT", clear
gen year = 2003
keep if linkflag==1
keep if panel03==8
gen nhis_yr = srvy_yr
gen str2 fmx = "01"
egen publicid = concat(nhis_yr hhx fmx px)
gen panel = 8
save "$mepsout/his_meps_2003.dta", replace

* 2004 (panel 9)
infix str dupersid 1-8 str hhx 9-14 str px 15-16 linkflag 17-17 panel04 18-18 srvy_yr 19-22 using "$mepsdir/MEPS_2004/NHMEP04X.DAT", clear
gen year = 2004
keep if linkflag==1
keep if panel04==9
gen nhis_yr = srvy_yr
gen str2 fmx = "01"
egen publicid = concat(nhis_yr hhx fmx px)
gen panel = 9
save "$mepsout/his_meps_2004.dta", replace

* 2005 (panel 10)
infix str dupersid 1-8 str hhx 9-14 str px 15-16 linkflag 17-17 panel05 18-19 srvy_yr 20-23 str fmx 24-25 str fpx 26-27 using "$mepsdir/MEPS_2005/NHMEP05X.DAT", clear
gen year = 2005
keep if linkflag==1
keep if panel05==10
gen nhis_yr = srvy_yr
egen publicid = concat(nhis_yr hhx fmx fpx)
gen panel = 10
save "$mepsout/his_meps_2005.dta", replace

* 2006 (panel 11)
infix str dupersid 1-8 str hhx 9-14 str fmx 15-16 str fpx 17-18 linkflag 19-19 panel06 20-21 srvy_yr 22-25 using "$mepsdir/MEPS_2006/NHMEP06X.DAT", clear
gen year = 2006
tab linkflag
tab srvy_yr
keep if linkflag==1
keep if panel06==11
gen nhis_yr = srvy_yr
egen publicid = concat(nhis_yr hhx fmx fpx)
gen panel = 11
save "$mepsout/his_meps_2006.dta", replace

* 2007 (panel 12)
infix str dupersid 1-8 str hhx 9-14 str fmx 15-16 str fpx 17-18 linkflag 19-19 panel07 20-21 srvy_yr 22-25 using "$mepsdir/MEPS_2007/NHMEP07X.DAT", clear
gen year = 2007
keep if linkflag==1
keep if panel07==12
gen nhis_yr = srvy_yr
egen publicid = concat(nhis_yr hhx fmx fpx)
gen panel = 12
save "$mepsout/his_meps_2007.dta", replace

* 2008 (panel 13)
infix str dupersid 1-8 str hhx 9-14 str fmx 15-16 str fpx 17-18 linkflag 19-19 panel08 20-21 srvy_yr 22-25 using "$mepsdir/MEPS_2008/NHMEP08X.DAT", clear 
gen year = 2008
keep if linkflag==1
keep if panel08==13
gen nhis_yr = srvy_yr
egen publicid = concat(nhis_yr hhx fmx fpx)
gen panel =13
save "$mepsout/his_meps_2008.dta", replace


append using "$mepsout/his_meps_2007.dta"
append using "$mepsout/his_meps_2006.dta"
append using "$mepsout/his_meps_2005.dta"
append using "$mepsout/his_meps_2004.dta"
append using "$mepsout/his_meps_2003.dta"
append using "$mepsout/his_meps_2002.dta"
append using "$mepsout/his_meps_2001.dta"
append using "$mepsout/his_meps_2000.dta"
append using "$mepsout/his_meps_1999.dta"
append using "$mepsout/his_meps_1998.dta"
append using "$mepsout/his_meps_1997.dta"
append using "$mepsout/his_meps_1996.dta"


sort dupersid year
tab year

* merge NHIS smoking data
isid publicid



merge 1:1 publicid using $mepsout/nhis-smoke.dta
drop if _merge ==2
tab nhis_yr smokev
drop _merge
save "$mepsout/his_meps_all.dta", replace


