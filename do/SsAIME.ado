
/* Computation of Averaged Indexed Monthly Earning 
based on Social Security Handbook
						for United States -
Based on discussions with N. Maestas
and Social Security Handbook.
* Version 2. Dec 2007, Add comments by Yuhui
---------------------------------------------------- */


/// 3. Calculation of AIME

 
cap program drop SsAIME
program define SsAIME, sortpreserve 
	version 9
	syntax varlist [if] [in], gen(name) earn(string) fyr(integer) yr(integer) 
	marksample touse
	
	/*
	"name": the variable to generate (AIME)
	"earn": beginning name of earnings history, should be in the format of "earn"YYYY
	"fyr" : A scalar of beginning year of earnings history
	"cyr" : A scalar of the year of calculating AIME
	*/
	
	/// 0.  Confirm existence of variables from "fyr" to "yr", and also two acciliary datasets exist
	foreach t of numlist `fyr'/`yr' {
		  cap confirm var `earn'`t'
		  if _rc {
		  	dis "Missing earning history for year `t'"
		  	exit(333)
		  }
	} 
	foreach f in cnatwage earnlimit{ 
		cap confirm file `f'.dta
			if _rc{
					dis "Please put file `f'.dta under the working directory"
					exit (555)
			}
	}
	
	tokenize `varlist'	
	#delimit;
	tempvar rq rbyr rclyr;
	
	///1. Assign Data to Variables -------------------
	
	qui gen `rq'       = `1'  if `touse';	
	qui gen `rbyr'     = `2'  if `touse';	
	qui gen `rclyr'    = `3'  if `touse';
	
	/*  by YUHUI
			 `rq': total quarters of covered earnings
			 `rbyr': year of birth
			 `rclyr': claim year
	*/
			 
	/// 2.A Important Indicators ----------------------
	tempvar ry60 endyr;
	qui gen `ry60' = `rbyr' + 60 if `touse';
	qui gen `endyr' = min(`rclyr',`yr') if `touse';
	
	/// 2.B Earning History ---------------------------
	capture drop y60;
	gen y60 = `ry60' if `touse';
	sort y60;
	capture drop _merge;
	qui merge y60 using "cnatwage.dta", nokeep;	
	drop _merge y60;
	
	/* by YUHUI
		cnatwage.dta includes wage threshold from yr1960-2060
		"w" is the threshold in a certain year
		"wYYYY" is the wide format of threshold in a certain year, constant across rows
	*/
	
	/// 2.C. Allowed Earnings --------------------------
	tempvar ecap;
	gen yr = `yr' if `touse';
	sort yr;
	capture drop _merge;
	qui merge yr using "earnlimit.dta", nokeep;	
	drop _merge yr cap;
	
	/// 2.D. Indexation --------------------------------
	tempvar yrn;
	qui gen `yrn' = `ry60'-`fyr'+1 if `touse';

foreach t of numlist `fyr'/`yr' {;
    tempvar rw`t';
    qui gen `rw`t'' = `earn'`t' if `touse';
    qui replace `rw`t'' = min(`rw`t'',c`t')*w/w`t' if `touse'&`ry60'>=`t'& `rclyr'>= `t';	
    qui replace `rw`t'' = min(`rw`t'',c`t') if `touse'&`ry60'<`t' & `rclyr' >= `t';	
    qui replace `rw`t'' = 0 if `rclyr'< `t'  ;
  };
  
	/* by YUHUI
	  the "
	  if not older than age 60, taxable earnings will be inflated to the year at which one reaches 60
	  if older than 60, then no inflation is made
	  if year later than claiming year, then earning is zero
	  the year of claiming might include some earnings
	*/

  drop w1* c1* w2* c2* `yrn' w;
  /// 2.E. Computation Years -------------------------- 
  tempvar rcpyr endyr2;
  qui gen `endyr2' = `endyr'-`fyr'+1 if `touse';
  qui gen `rcpyr' = min(35,`endyr2')  if `touse';
  
  global w "";
  forvalues i = `fyr'/`yr'{; 
			global w $w `rw`i'' ; 
	}; 
	

  /// 2.F. Computation AIME   -------------------------- 
  cap drop aime_r;
  mata: caime("$w","`rcpyr'","`endyr2'");
	capture drop `gen';

///	list aime in 1;
	qui gen `gen' = aime_r/12 if `touse';
	drop aime_r ; 
  macro drop w; 
  
///	list `gen' in 1;
	mata: mata clear;
	#delimit cr
end


version 9
 mata:
 void caime(earnings, rcpyr, endyr2)
 {  
    mW = st_data(.,tokens(earnings))
    vC = st_data(.,rcpyr)
    vE = st_data(.,endyr2)
    N = rows(mW)
    mA = J(N,1,0)
    for (i=1;i<=N;i++){									
			endyr = vE[i]
			cyr = vC[i]
			mwi = sort(mW[i,1::endyr]',-1) 
			
			mA[i] =	sum(mwi[1::cyr])/35
 			
    } 
    idx = st_addvar("double","aime_r")
    st_store(range(1,N,1),"aime_r",mA);  			
 }   

 end
 
