{smcl}
{* *! version 1.2.1  017aug2018}{...}
{viewerjumpto "Syntax" "limpe_re##syntax"}{...}
{viewerjumpto "Description" "limpe_re##description"}{...}
{viewerjumpto "Options" "limpe_re##options"}{...}
{viewerjumpto "Examples" "limpe_re##examples"}{...}
{title:Title}

{phang}
{bf:limpe_re} {hline 2} Estimate the linear-in-means peer effects model with random group effects.

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:limpe_re}
[{varlist}]
{ifin}
[{cmd:,} {it:options}]

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{p2coldent: * {opth groupid(varname)}} group identifier variable {p_end}
{synopt:{opth avevar(varlist)}} specify variables whose leave-out means are added as covariates. {p_end}
{synopt:{opt nocons:tant}} suppresses the constant term in the model. {p_end}
{synopt:{opt from()}} specify starting point of the optimization {p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}
* {cmd:groupid()} is required.{p_end}
{p 4 6 2}
{cmd:by} is allowed; see {manhelp by D}.{p_end}


{marker description}{...}
{title:Description}

{pstd}
{cmd:limpe_re} performs the quasi-maximum likelihood estimation for a linear-in-means peer effects model with random group effects, proposed by Kuersteiner, Prucha and Zeng (2021). Standard errors are calculated using the procedure desribed in the paper. Factor variables are currently not allowed in independent variables. You should generate dummies before adding them to the model. The program runs slowly when the sample size is large. It is therefore recommended to test your code on a subset of the groups before applying it to the whole sample. This version is created on Jul 20th,2021. 

{marker options}{...}
{title:Options}

{dlgtab:Main}

{phang}
{opt groupid} specifies a variable that contains a unique
        identifier for each group. The variable can be numeric or character. id() is required.

{phang}
{opt avevar(varlist)} specify a set of variables whose leave-out means (or spatial lags) are added as covariates. These are exogenous peer effects. 

{phang}
{opt noconstant} suppresses the constant term in the model. 

{phang}
{opt from()} specifies the starting point of optimization. The default is {cmd: from(lambda=0 stde=0.5 stda=0.5)}, so starding point is set to lambda=0, stde=0.5*std(y), stda=0.5*std(y), std(y) is the standard deviation of the dependent variable, stde specifies the starting value of the standard error of individual error terms, stda specifies the starting value of the standard error of random group effects. 

{marker examples}{...}
{title:Examples}
{phang}The example below estimates class level peer effects using the kindergarten sample of Project STAR. Control variables include poor, black, girl, the leave-out-means of these three variables, school fixed effects, and class type fixed effects. The starting point of optimization is set to lambda=0, stde=0.3*std(y), stda=0.2*std(y), where y is mathnorm (normalized math score). 

{phang2}{cmd:. use starpanel,clear}{p_end}
{phang2}{cmd:. tab schid,gen(Dsch)}{p_end}
{phang2}{cmd:. tab classtype,gen(Dtype)}{p_end}	
{phang2}{cmd:. limpe_re mathnorm poor black girl tchms  Dsch* Dtype* if grade==0, groupid(class) avevar(poor black girl) from(lambda=0 sige=0.3 siga=0.2) }{p_end}

{title:Saved results}

{pstd}
{cmd:limpe_re} saves the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(R)}}number of groups{p_end}
{synopt:{cmd:e(lnf0)}}value of the log likelihood function at the starting point{p_end}
{synopt:{cmd:e(lnf)}}value of the log likelihood at the optimum {p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrix}{p_end}
{synopt:{cmd:e(init)}}vector of starting point of optimization{p_end}
{synopt:{cmd:e(b)}} parameter vector {p_end}
{synopt:{cmd:e(V)}} variance-covariance matrix {p_end}

{title:Authors}

{phang}Guido Kuersteiner{p_end}
{phang}Department of Economics{p_end}
{phang}University of Maryland{p_end}
{phang}College Park, MD{p_end}
{phang}{browse "mailto:gkuerste@umd.edu":gkuerste@umd.edu}

{phang}Ingmar Prucha{p_end}
{phang}Department of Economics{p_end}
{phang}University of Maryland{p_end}
{phang}College Park, MD{p_end}
{phang}{browse "mailto:prucha@umd.edu":prucha@umd.edu}

{phang}Ying Zeng{p_end}
{phang}School of Economics and Wang Yanan Institute for Studies in Economics {p_end}
{phang}Xiamen University {p_end}
{phang}Xiamen, China{p_end}
{phang}{browse "mailto:zengying17@gmail.com":zengying17@gmail.com}{p_end}


{title:References}

{phang}Kuersteiner, G. M., Prucha, I. R., & Zeng, Y. (2021). Efficient Peer Effects Estimators with Random Group Effects. ArXiv:2105.04330 [Econ, Stat]. http://arxiv.org/abs/2105.04330
