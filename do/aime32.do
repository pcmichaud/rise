* AIME replacement rates (as in French)
insheet par using params/input/earnings.csv, clear

* send pars to matrices
mkmat par, matrix(par)
matrix par = par'
matrix colnames par = age age2 one rho sige sigv
matrix list par

* make space for panel
set obs 5000

* horizon T
global T = 2010

* age at start
gen age = 20
gen age2 = age*age

* growth rate in economy
global g = 0.0

* initial eps
gen eps = rnormal(0,par[1,5]/(1-par[1,4]^2)) 

forvalues t=1960/$T {
	qui gen earn`t' = .
		qui replace eps = par[1,4]*eps + sqrt(par[1,5])*rnormal() 
		qui replace earn`t' = exp(par[1,1]*age + par[1,2]*age2 + par[1,3] + eps)
		qui replace earn`t' = min(earn`t',98500)
	qui replace age = age+1
	qui replace age2 = age*age
}

gen dearn = earn2004 - earn2003

// Compute highest 35 years at each year
mata: 
	st_view(W=.,.,tokens("earn*"))
	N = rows(W);
	T = cols(W);
	TrueA = J(N,T,.);
	ApprA = J(N,T,.);
	alpha = J(N,T,0);
	W = W:/12;
	for (i=1;i<=N;i++) {
		for (t=1;t<=T;t++) {
			if (t<35) {
				cap = t
				TrueA[i,t] = sum(W[i,1::t])/35;
			}
			else {
			
				if (W[i,t]>min(W[i,1::t-1])) {
					//(W[i,t],min(W[i,1::t-1]))
				}	
				cap = 35
				x = sort(W[i,1::t]',-1);
				//x
				TrueA[i,t] = sum(x[1::cap])/35;
						
			}		
		
			
			if (t==1) {
				ApprA[i,t] = W[i,t]/35;
			}
			else {
				ApprA[i,t] = (1+$g):*ApprA[i,t-1] + (W[i,t] - (1+$g):*ApprA[i,t-1])/35;
			}
			if (t>1) {
				if (TrueA[i,t]<=TrueA[i,t-1]) {	
					alpha[i,t] = 1 
				}
			}

		}
			
	}
	mTrueA = mean(TrueA)
	mApprA = mean(ApprA)
	alpha = mean(alpha[1::N,36::51])
	//alpha = mean(alpha[1::N,])
	(alpha')
	st_matrix("rep",(alpha'))

end

preserve
	svmat rep
	keep rep1
	order rep1
	keep if _n<=16
	outsheet using params/input/rep.csv, replace comma nonames nolabel
restore

gen id = _n
reshape long earn@, i(id) j(year)
drop age

gen age = year - 1940

table age, content(mean earn)

exit



 
 
