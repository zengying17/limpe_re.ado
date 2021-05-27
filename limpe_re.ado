*! version 2.0 May 27th, 2021.
program define limpe_re, eclass sortpreserve 
version 13.1
syntax varlist(numeric) [if] [in] , ///
   groupid(varname numeric) ///		 
   [avevar(varlist) /// 		 
   noCONStant ///
   from(passthru) ///
   ]
marksample touse, novarlist   
markout `touse' `varlist' `avevar' `groupid'

preserve
qui keep if `touse' 

tempvar grsize grfirst
qui bysort `groupid':egen `grsize'=count(`groupid')
qui bysort `groupid':gen `grfirst'=_n==1
qui count if `grsize'<=1

local singlegr=r(N)     //Drop Singleton groups 
if `singlegr'>0{        
   qui drop if `grsize'<=1
   di "`singlegr' observations ignored because group size equals to 1"
}
//Adjust the Variables.
gettoken depvar indeps : varlist    
_fv_check_depvar `depvar'
		
foreach x of local avevar{				
   tempvar tot`x'
   bysort `groupid':egen double `tot`x''=total(`x')
   gen double _ave_`x'=(`tot`x''-`x')/(`grsize'-1) 		//genearte leave-out-mean of variables
   local avevar2 `avevar2' _ave_`x'				     
}
local indeps "`indeps' `avevar2'"

qui _rmcoll `indeps' if `touse',noconstant expand
local omitk=r(k_omitted)
if `omitk'>=1 di "Note:`omitk' variable dropped due to collinearty."
local indnames `r(varlist)'

if "`constant'" == "" {
tempvar cons
qui gen `cons'=1
local indeps "`indeps' `cons'" 		//Define Independent Variables. 
local indnames "`indnames' _cons"
}
local groupdum "ibn.`groupid'"
//add mean and group deviation
foreach x of local indeps{
	bysort `groupid':egen double _mean_`x'=mean(`x')
	gen double _demean_`x'=`x'-_mean_`x'
	local meanx "`meanx' _mean_`x'"
	local demeanx "`demeanx' _demean_`x'"
}

	
tempvar meany demeany 
bysort `groupid':egen double `meany'=mean(`depvar')

gen double `demeany'=`depvar'-`meany'

//Set Macros for initial value and grid search

qui sum `depvar'
	local sdy=r(sd)	
	tempname init initp 
	if `"`from'"' != "" {
		_mkvec `init', `from' error("from()") 	
	}
	else matrix `init'=(0,0.5*`sdy',0.5*`sdy') 
	mata:`initp'=st_matrix("`init'") 

mata: remain("`meany'","`demeany'","`meanx'","`demeanx'","`grsize'","`grfirst'","`groupdum'","`touse'",`initp') 

tempname b V
matrix `b'=r(b)
matrix `V'=r(V)
matname `b' lambda std_e std_a `indnames',explicit c(.)
matname `V' lambda std_e std_a `indnames',explicit

ereturn post `b' `V' ,esample(`touse')
ereturn local  cmd  "limpe_re"
ereturn scalar R=r(R)
ereturn scalar N=r(N)
ereturn scalar lnf0=r(lnf0)
ereturn scalar lnf=r(lnf)
matrix init=r(init) 
ereturn matrix init=init

/*
matrix V1=r(V1)
matrix V2=r(V2)
matrix V3=r(V3)
matrix g=r(g)
matrix gp=r(gp)
ereturn matrix V1=V1
ereturn matrix V2=V2
ereturn matrix V3=V3
ereturn matrix g=g
ereturn matrix gp=gp
*/
ereturn display
qui capture drop _ave_*
restore
end

//Part 2: Define the Mata reopt for Objective Function. 

mata:
void reopt(todo,p,Ybar,Ystar,Xbar,Xstar,mn,mr,en,er,lnf,g,H)
{
lambda=p[1,1]
vare=p[1,2]^2
vara=p[1,3]^2
N=rows(mn)
R=rows(mr)
P_o=diag(en:/(vare*en+mn*vara))
p_s=((lambda-1)*en+mn):/(mn-en)
beta=invsym(1/vare*Xstar'*Xstar+Xbar'*P_o*Xbar)*(1/vare*Xstar'*(p_s:*Ystar)+(1-lambda)*Xbar'*P_o*Ybar)
ustar=p_s:*Ystar-Xstar*beta
ubar=Ybar*(1-lambda)-Xbar*beta
lnf=-0.5*N*log(2*pi())+R*ln(abs(1-lambda))+(mr-er)'*ln(er+lambda*er:/(mr-er))  ///
    -0.5*N*ln(vare) -0.5*er'*ln(er+vara/vare*mr) ///
	-0.5*(ustar'*ustar/vare+ubar'*P_o*ubar)

   
if (todo>=1){
g1=-lambda/(1-lambda)*mr'*(1:/(mr+(lambda-1)*er))-ustar'*(Ystar:/(mn-en))/vare+ubar'*P_o*Ybar
g2=(-0.5*N/vare+0.5*er'*((vara/vare^2*mr):/(er+vara/vare*mr))+0.5*ustar'*ustar/vare^2+0.5*ubar'*P_o*P_o*ubar)*2*p[1,2]
g3=(-0.5/vare*mr'*(er:/(er+vara/vare*mr))+0.5*ubar'*diag(mn:/((vare*en+mn*vara):^2))*ubar)*2*p[1,3]

g=(g1,g2,g3)
}

}

//Part 3: Define the Mata re_main for Optimization. 

void remain( ///
string scalar meany, 		///  
string scalar demeany, 		/// 			
string scalar meanxs, 		///
string scalar demeanxs, 		/// 
string scalar grsize, 	///
string scalar grfirst, ///
string scalar groupdum, ///
string scalar tousename, 	///
real vector init			
)
{
Ybar=st_data(.,meany,tousename)
Ystar=st_data(.,demeany,tousename)
Xbar=st_data(.,meanxs,tousename)
Xstar=st_data(.,demeanxs,tousename)
C=st_data(.,groupdum,tousename)
mn=st_data(.,grsize,tousename)
firstind=st_data(.,grfirst,tousename)
mr= select(mn, firstind)

N=rows(mn)
R=rows(mr)
en=J(N,1,1)
er=J(R,1,1)


//The Optimization Process
S=optimize_init()
optimize_init_evaluator(S, &reopt())
optimize_init_evaluatortype(S,"gf1")
optimize_init_params(S,init)
optimize_init_which(S,"max")
optimize_init_argument(S,1,Ybar)
optimize_init_argument(S,2,Ystar)
optimize_init_argument(S,3,Xbar)
optimize_init_argument(S,4,Xstar)
optimize_init_argument(S,5,mn)
optimize_init_argument(S,6,mr)
optimize_init_argument(S,7,en)
optimize_init_argument(S,8,er)
p=optimize(S)
V1 =optimize_result_V_oim(S)

lnf=optimize_result_value(S)
lnf0=optimize_result_value0(S)

lambda=p[1,1]
vare=p[1,2]^2
vara=p[1,3]^2

p=(lambda,abs(p[1,2]),abs(p[1,3])) 

P_o=diag(en:/(vare*en+mn*vara))

p_s=((lambda-1)*en+mn):/(mn-en)


beta=invsym(1/vare*Xstar'*Xstar+Xbar'*P_o*Xbar)*(1/vare*Xstar'*(p_s:*Ystar)+(1-lambda)*Xbar'*P_o*Ybar)

Vbeta=invsym(1/vare*Xstar'*Xstar+Xbar'*P_o*Xbar)

ustar=p_s:*Ystar-Xstar*beta
ubar=Ybar*(1-lambda)-Xbar*beta


g1=-lambda/(1-lambda)*mr'*(1:/(mr+(lambda-1)*er))-ustar'*(Ystar:/(mn-en))/vare+ubar'*P_o*Ybar
g2=(-0.5*N/vare+0.5*er'*((vara/vare^2*mr):/(er+vara/vare*mr))+0.5*ustar'*ustar/vare^2+0.5*ubar'*P_o*P_o*ubar)*2*p[1,2]
g3=(-0.5/vare*mr'*(er:/(er+vara/vare*mr))+0.5*ubar'*diag(mn:/((vare*en+mn*vara):^2))*ubar)*2*p[1,3]
g=(g1,g2,g3)

gp1=(-lambda)/(1-lambda)*(mr:/(mr+(lambda-1)*er))-C'*((ustar:*(Ystar:/(mn-en))/vare)-(ubar:*(P_o*Ybar)))
gp2=(-0.5/vare*mr+0.5*vara/vare^2*mr:/(er+vara/vare*mr)+0.5*C'*((ustar:*ustar)/vare^2+(ubar:*(P_o*P_o*ubar))))*2*p[1,2]
gp3=(-0.5*mr:/(vare*er+vara*mr)+0.5*C'*(ubar:*(diag(mn:/((vare*en+mn*vara):^2))*ubar)))*2*p[1,3]
gp=(gp1,gp2,gp3)  //3*R matrix 


V2=invsym(gp'*gp)
V3=V1*(gp'*gp)*V1
st_rclear()
st_matrix("r(V1)",V1)
st_matrix("r(V2)",V2)
st_matrix("r(V3)",V3)
st_matrix("r(g)",g)
st_matrix("r(gp)",gp)

st_matrix("r(b)",(p,beta'))
st_matrix("r(V)",blockdiag(V3,Vbeta))
st_numscalar("r(N)", N)      
st_numscalar("r(R)", R)
st_matrix("r(init)",init)
st_numscalar("r(lnf)", lnf)      
st_numscalar("r(lnf0)", lnf0)            
}
end	
