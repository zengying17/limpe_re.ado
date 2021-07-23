*! version 2.1 Jul 20th, 2021.
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
qui sum `grsize'
local maxm=r(max)

mata: remain("`meany'","`demeany'","`meanx'","`demeanx'","`grsize'","`grfirst'","`touse'",`initp',`maxm') 

tempname b V
matrix `b'=r(b)
matrix `V'=r(V)
matname `b' lambda vare vara `indnames',explicit c(.)
matname `V' lambda vare vara `indnames',explicit

ereturn post `b' `V' ,esample(`touse')
ereturn local  cmd  "limpe_re"
ereturn scalar R=r(R)
ereturn scalar N=r(N)
ereturn scalar lnf0=r(lnf0)
ereturn scalar mu3e=r(mu3e)
ereturn scalar mu4e=r(mu4e)
ereturn scalar mu3a=r(mu3a)
ereturn scalar mu4a=r(mu4a)

matrix init=r(init) 
ereturn matrix init=init


matrix V1=r(V1)
matrix V2=r(V2)
matrix V3=r(V3)

ereturn matrix V1=V1
ereturn matrix V2=V2
ereturn matrix V3=V3

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
lnf=-0.5*N*log(2*pi())+R*ln(abs(1-lambda))+(mr-er)'*ln (er+lambda*er:/(mr-er))  ///
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
string scalar tousename, 	///
real vector init, ///
real scalar maxm			
)
{
Ybar=st_data(.,meany,tousename)
Ystar=st_data(.,demeany,tousename)
Xbar=st_data(.,meanxs,tousename)
Xstar=st_data(.,demeanxs,tousename)
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
ustar=p_s:*Ystar-Xstar*beta
ubar=Ybar*(1-lambda)-Xbar*beta

//First calculate the third and fourth moments
u2star=select(ustar,mn:>=3)
u2bar=select(ubar,mn:>=3)
R2=colsum(mr:>=3) //number of grooups with size greater than or equals to 3
m2n=select(mn,mn:>=3)
e2n=select(en,mn:>=3)

fe3=(u2star:^3):/(m2n-3*e2n+2:/m2n)
mu3e=colsum(fe3)/R2

fa3=(u2bar:^3):/m2n-fe3:/(m2n:^2)
mu3a=colsum(fa3)/R2

fe4=(m2n:^3):/(m2n:^3-4*m2n:^2+6*m2n-3*e2n):*(u2star:^4:/m2n-3*(m2n-e2n):*(2*m2n-3*e2n):/(m2n:^4)*vare^2)
mu4e=colsum(fe4)/R2

fa4=u2bar:^4:/m2n-fe4:/(m2n:^3)-3*(m2n-e2n):/(m2n:^4)*vare^2-6:/(m2n:^2)*vare*vara
mu4a=colsum(fa4)/R2

//calcalating the V-C matrix using H and g for the nonconcentrated MLE
kx=cols(Xstar)
k0=J(kx,1,0)
ki=I(kx)
Ups=J(kx+3,kx+3,0)
Gam=J(kx+3,kx+3,0)

for (m=2;m<=maxm;m++){
Rm=colsum(mr:==m) //number of groups with size equals to m
if (Rm>0){
Dm=select(Xstar,mn:==m)
Dm2=Dm'*Dm
Bm=select(Xbar,mn:==m)
Bm2=Bm'*Bm
Cm=colsum(Bm)/m

phim=(1/((m-1+lambda)*vare),-m/((1-lambda)*(vare+m*vara)),1/((m-1+lambda)*vare)*beta',-m/((1-lambda)*(vare+m*vara))*beta'\ ///
	 -1/(2*vare^2), -m/(2*(vare+m*vara)^2),k0',k0'\ ///
	 0,-m^2/(2*(vare+m*vara)^2),k0',k0'\ ///
	 k0,k0,-1/vare*ki,-m/(vare+m*vara)*ki)
	 	 
PsiGm=blockdiag(2*(m-1)*vare^2*Rm/N,blockdiag(2*(vara+vare/m)^2*Rm/N,blockdiag(vare*Dm2/N,(vara+vare/m)*(Bm2/m)/N)))

Psi11=Rm/N*blockdiag(2*(m-1)*vare^2,2*(vara+vare/m)^2+(mu4a-3*vara^2))+Rm/N*(mu4e-3*vare^2)*((m-1)^2/m,(m-1)/m^2\(m-1)/m^2,1/m^3)
Psi21=(k0,k0\(m-1)/m*mu3e*Cm'/N,(mu3a+mu3e/m^2)*Cm'/N)
Psi22=blockdiag(vare*Dm2/N,(vara+vare/m)*(Bm2/m)/N)

Psim=(Psi11,Psi21'\Psi21,Psi22) 
Upsm=phim*Psim*phim'
Gamm=phim*PsiGm*phim'

Ups=Ups+Upsm
Gam=Gam+Gamm
}
}
	
V2=invsym(Gam)/N
V3=makesymmetric(N*(V2*Ups*V2)) //caution: mata precision problem detected.

st_rclear()
st_matrix("r(V1)",V1)
st_matrix("r(V2)",V2)
st_matrix("r(V3)",V3)

st_matrix("r(b)",(lambda,vare,vara,beta'))
st_matrix("r(V)",V3)
st_numscalar("r(N)", N)      
st_numscalar("r(R)", R)
st_matrix("r(init)",init)
st_numscalar("r(lnf)", lnf)      
st_numscalar("r(lnf0)", lnf0)
st_numscalar("r(mu3e)",mu3e)
st_numscalar("r(mu4e)",mu4e) 
st_numscalar("r(mu3a)",mu3a) 
st_numscalar("r(mu4a)",mu4a)         
}
end	
