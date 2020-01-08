
integer nages, nh, ne, nz, nf, nd, name, nw  			! grid sizes
integer nsim, ndata, nind								! size of data sets
double precision ratio_obs								! ratio of obs simulated and in data
integer npar											! number of parameters
integer nmoments, nfreemoments, ngroup					! number of parameters and moments
character*80 scenario									! name of scenario (input)
integer scn(4)											! scenario identifiers (for tags in program)
integer nfreepar										! number of free parameters in estimation
integer, allocatable:: ages(:), gridz(:)				! grid for ages
integer agemin, zmin, zmax, gapz, zbase					! starting age
double precision xmin, leisure_endow, rrate				! institutions			
double precision wmin, wmax, amemin, amemax				! bounds for grids
double precision c_wmin, c_wmax							! bounds for grid (transformed)
character*2 alg 										! choice of algorithm
integer nmaxstates, nspace,nstatevars, ndecisionvars	! variables for solving
integer rmin, rmax, mcare								! ages for working and medicare eligibility
integer ngrid											! number of points on grid for continuous vars (ame and w)
integer nobs											! number of observations in simulated dataset
character*80 path 										! path for the project (set in estimation.info)
double precision hours_worked							! hours worked full time
double precision, allocatable :: cashmin(:,:)				! floor for cash-on-hand
double precision  umin									! minimum utility (to avoid undeflow)
integer ier, rank, numprocs, numworkers, status(MPI_STATUS_SIZE)	! parameters for MPI
logical ismaster										! swtich to define who is master process
real*8 starttime, endtime, totaltime					! to record time spent solving	
double precision, allocatable :: gridw(:), c_gridw(:)	! grid for wealth
double precision, allocatable :: gridame(:)				! grid for AME
double precision gapw, gapame							! gaps for both grids
double precision mcdcopay		 						! Medicaid co-pay
double precision medmax									! maximum on how much can spend on health
integer, allocatable :: nstates(:),states(:,:,:)		! arrays to select possible states
double precision, allocatable :: copay(:,:,:),insgen(:,:)	! co-pays for insurance
double precision arf, drc, earntaxgain 					! social security parameters
double precision, allocatable :: sstax(:,:),ssmax(:,:),mctax(:,:),mcmax(:,:) 	! social security parameters
integer npiabend, npiarate, nearntax								! social security parameters
double precision, allocatable :: piabend(:),piarate(:), rep(:), ssgen(:) 	! social security parameters
integer era, nra 													! social security parameters
double precision, allocatable :: earntaxrate(:), taxpar(:,:,:)		! taxes
double precision, allocatable :: hdraws(:,:), edraws(:,:), idraws(:)			! draws for estimation
character*80 buffer													! buffer for reading names
integer origin, now 												! years for scenarios
double precision wages, insurance, other
character*12 fmt													! format for real numbers
! auxiliary processes
double precision parearn(5), parother(4)
integer, allocatable :: progress(:,:)
double precision, allocatable :: probe(:,:), cumprobe(:,:), pte(:)
double precision incsig
double precision, allocatable :: earn(:,:), otherinc(:,:)
double precision par_adjust(3)
double precision, allocatable :: lambda(:,:,:,:), obese(:,:), smoke(:,:), adjust(:)
double precision, allocatable :: parh(:,:), pard(:), theta(:,:),probs(:,:)
double precision parh_risk(2,3) 
double precision targ_medexp, targ_e25
parameter (targ_medexp =1031.57, targ_e25 = 70.0)
integer nparh, npard
logical isave
integer nsave 
integer, allocatable :: saveages(:)
double precision curv
integer ngridc, ngridm												! number of grid points c,m
logical iweight, iload

! type for moments
type mom
	character*10 label					! label of moment
	integer nages						! number of ages for that moment
	integer, allocatable :: ages(:)		! list of ages at which moment computed
	logical nh 							! whether by health or not
	integer nmom 						! total number of moments for this moment (computed)
	logical free						! whether moment used in estimation
	double precision, allocatable :: sim(:,:)	! simulated value
	integer, allocatable :: nsim(:,:)			! number of simulated folks
end type mom

type (mom), allocatable :: moments(:)
double precision, allocatable :: data(:,:),datamean(:), datans(:), optw(:,:), datacov(:,:), select(:,:)

type scn_details
	logical iprogress
	double precision progress
	logical iother
	double precision other	
	logical iwages
	double precision wages
	logical itaxes
	integer taxes_yr
	logical icashmin
	integer cashmin_yr
	logical issben
	integer ssben_yr
	logical icoverage
	double precision coverage
	logical icopay_mc
	double precision copay_mc
	logical icopay_ei
	double precision copay_ei	
	logical iprices
	double precision prices
	logical idwage
	double precision dwage
	logical idcopay
	double precision dcopay
end type scn_details

type (scn_details) scn_sheet					! details of scenario

! type for parameters
type params
	character*20 label					! label of parameter
	logical free						! is free parameter
	integer ipos						! index in vector of free parameters
	double precision lb 				! lower bound for simulated annealing
	double precision ub					! upper bound for simulated annealing
	double precision step 				! step for Nelder-Mead
	double precision value				! estimated value (trial or final)
	double precision serr				! standard error
	double precision pvalue				! p-value
	double precision minci				! minimum value in 95% confidence interval
	double precision maxci				! maximum value in 95% confidence interval
end type params

type (params), allocatable :: g_initpar(:)

type preferences
	double precision beta				! discount factor
	double precision sigma				! risk aversion
	double precision phi(2)				! leisure effect of health
	double precision alpha(4)			! baseline utility levels
	double precision psi 				! share of consumption in utility
	double precision kapa(3)			! bequest function parameters
	double precision age(3)				! age effects on production function
	double precision theta(3)			! productivity on medical expenditures
	double precision theta2(3)			! curvature productivity on medical expenditures
	double precision gamma(3,3)			! effect of current health on transitions
	double precision base(3)			! base constants in production function
	double precision leisure_endow		! leisure endowment
	double precision prog_med			! annual rate progress in productivity
	double precision prog_good			! annual rate progress in productivity
	double precision prog_vgood			! annual rate progress in productivity
	double precision prog_surv			! annual rate progress survival
end type preferences

! optimal decision rules
type rules
	double precision value
	double precision cons
	double precision medexp
	double precision claim
	double precision work	
	double precision cash
	double precision options(4)
end type rules

! simulated data
type sim 
	integer id
	integer age
	integer z
	double precision wealth
	integer work
	integer claim
	double precision income
	integer insurance
	double precision ame
	double precision consumption
	integer dead
	integer health
	integer base_health
	double precision oop
	double precision medexp
	double precision value
	double precision tr
end type sim


! for passing easily scalars inside optimization routine
type spendargs
	integer a
	integer h
	integer e
	integer z
	integer f
	integer d
	double precision cash
	double precision tr
	double precision nextame
	double precision cmin
	double precision cmax
	double precision mmin
	double precision mmax
	integer ff
	integer fff
	integer dd
	integer rr
end type spendargs


