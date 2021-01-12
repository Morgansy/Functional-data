#source files
source("LMCcopula.R")

set.seed(123)
s1=runif(10,-3,3)
s2=runif(10,-3,3)
#s1,s2 - Eucledian coordinates
#DistM - matrix of Eucledian distances 
Dist=DistM(s1,s2);

#########################################
### LMC copula model parameterization ###
#########################################

#LMC model for the Gaussian process
#zz1=par[10]*z0+sqrt(1-par[10]^2)*z1;
#zz2=par[11]*z0+sqrt(1-par[11]^2)*z2;
#z0,z1,z2-three independent Gaussian processes

#exponential factor model 
#W1=zz1+par[1]*E0u+par[3]*E1u-par[5]*E0l-par[4]*E1l; 
#W2=zz2+par[2]*E0u+   b2u*E2u-par[6]*E0l-   b2l*E2l;
#E0u,E1u,E0l,E1l,E2u,E2l - independent exponential factors

#covariance function for z0: exp(-par[7]*Dist^par[12])
#covariance function for z1: exp(-par[8]*Dist^par[13])
#covariance function for z2: exp(-par[9]*Dist^par[14])

####################################################
### data simulation and MLE for LMC copula model ###
####################################################

par=c(1.5,0.8,1.0,0.5,0.8,1.5, 0.6,0.7,0.8, 0.7,0.7, 1.4,1.5,1.6)

#simulate N=200 replicates of the process with a vector of parameters 'par'
#Dist - matrix of distances


#simulate from the model with Pareto factors instead of exponential factors
#model with a very strong tail dependence
#pr = T for Pareto factors
#pr = F for exponential factors
simdata=simexp3var(N=200,Dist,par,b2u=0.5,b2l=1,df=200,pr=T);

#uniform ranks for LMC copula model 
#uniform ranks should be used to get MLE
umat1=simdata$umat1
umat2=simdata$umat2

#umat1-variable1, umat2-variable2
#umat1, umat2 are matrices; rows are replicates; columns are locations
#Dist - matrix of distances
#LMC copula model MLE
#estimate misspecified model
 
#assuming exponential factors
MLE1=exp3varmle(umat1,umat2,Dist)



#par[10] reached the upper bound. We fix it (i.fix[10]=T)
i.fix=rep(F,14); i.fix[10]=T;
MLE1=exp3varmle(umat1,umat2,Dist,start=MLE1$estimate,i.fix=i.fix)

#now par[13] reached the lower bound. We fix it (i.fix[13]=T)
i.fix[13]=T;
MLE1=exp3varmle(umat1,umat2,Dist,start=MLE1$estimate,i.fix=i.fix)

#now par[8] reached the upper bound. We fix it (i.fix[8]=T)
i.fix[8]=T;
MLE1=exp3varmle(umat1,umat2,Dist,start=MLE1$estimate,i.fix=i.fix)


#output:
#> MLE1$minimum
#[1] -1766.205

#> MLE1$estimate
# [1] 0.58473533 0.41860523 0.57149010 0.06989667 0.44020262 0.95451069 0.66964047
# [8] 4.99999495 0.70920617 0.99899959 0.45766464 1.46663994 0.01000178 1.53966112



#simulate from the misspecified model
 
#Gaussian random field + exponential linear common factors
simdata0=simexp3var(N=10000,Dist,MLE1$estimate,b2u=0,b2l=0,df=200,pr=F);

umat10=simdata0$umat1

umat20=simdata0$umat2


#GoF for the misspecified model
 
#with exponential factors
#Spearman's rho
#cr1: -0.04/0.04

#cr2: -0.06/0.06

#cr12: -0.06/0.07

#fit: quite good
#lower tail
#lt1: 0.05/0.10

#lt2: 0.05/0.07

#lt12: 0.07/0.11
#fit: quite good
#upper tail
#ut1: 0.04/0.10

#ut2: 0.04/0.08

#ut12: 0.08/0.12

#fit: quite good


#now estimate parameters of the Pareto model 
#assuming Gaussian random field (par[1:6]=0)
MLE2=LMCvarmle(umat1,umat2,Dist,start=st.LMC,l.b=lb.LMC,u.b=ub.LMC,pl=1)

#output:
#> MLE2$minimum
#[1] -1635.735
#> MLE2$estimate
#[1] 0.2191684 0.3718317 0.6575322 0.6543466 0.8806295 0.8211629 0.8335379 1.9989999

#simulate from the misspecified Gaussian model
simdata0n=simexp3var(N=10000,Dist,c(rep(0,6),MLE2$estimate),b2u=0,b2l=0,df=200,pr=F);

umat10n=simdata0n$umat1

umat20n=simdata0n$umat2



#GoF for the misspecified Gaussian model
 
#Spearman's rho
#cr1: -0.02/0.07
#cr
2: -0.02/0.06
#cr12: -0.05/0.06
#fit: quite good
#lower tail
#lt1: 0.07/0.10
#lt2: 0.23/0.23
#lt12:0.34/0.34
#fit: bad! 
#dependence in the lower tail is underestimated! 
#upper tail
#ut1: 0.21/0.21
#ut2: 0.13/0.15
#ut12:0.03/0.08
#fit: quite bad. 
#dependence in the upper tail for var1 is underestimated!
 

#####################################################################
### predicted quantiles (on U(0,1) scale) for the estimated model ### 
#####################################################################

#for a new location (s1n,s2n)
s1n=0
s2n=0
#new matrix of distances, including the new location
pDist=DistM(c(s1n,s1),c(s2n,s2));
#conditional on the observed values u1, u2 (we select the first observation)
u1=umat1[1,]
u2=umat2[1,]

#par - estimated vector of copula parameters
par=MLE1$estimate
#nq - number of quadrature points for numerical integration
#nq = 25 gives very accurate results in most cases
#predicted 5% quantiles for both variables
q05=exp3qcondcdf(q=rep(0.05,2),u1,u2,pDist,par,nq=25)
#predicted medians for both variables
q50=exp3qcondcdf(q=rep(0.50,2),u1,u2,pDist,par,nq=25)
#predicted 95% quantiles for both variables
q95=exp3qcondcdf(q=rep(0.95,2),u1,u2,pDist,par,nq=25)

#> q05
#[1] 0.2590814 0.2893488
#> q50
#[1] 0.6378508 0.6653722
#> q95
#[1] 0.8976255 0.9302118



