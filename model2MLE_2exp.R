#library(mycopula)
library(mvtnorm)
library(cubature)

source("aux_functions.R")
source("mle_deriv_2exp.R")

#spatial Dist matrix of distances on a uniform grid
getDist = function(l.b,u.b,eps){
z0=seq(l.b,u.b,eps); 
d = length(z0);
Dist = matrix(0,d^2,d^2);
j1=rep(1:d,d);
j2=rep(1:d,each=d);
for(ii1 in 1:d^2){ 
  for(ii2 in 1:d^2){ 
    i1=floor((ii1-1)/d)+1; i2=floor((ii2-1)/d)+1; 
    Dist[ii1,ii2] = sqrt((z0[j1[ii1]]-z0[j1[ii2]])^2+(z0[j2[ii1]]-z0[j2[ii2]])^2);
  }
}
Dist
}

#spatial process generation; exponential covariance model (par - th); 
#common factor model 
#l.b, u.b - lower and upper bounds for X and Y coordinates of points to be generated
#eps - the grid size
genSpProcess = function(l.b,u.b,eps,rfm1,rfm2,th,par1,par2,cos=F,clsc=F){
z0=seq(l.b,u.b,eps); 
d = length(z0);
Dist = getDist(l.b,u.b,eps);
sigma = exp(-th*Dist);
if(cos == T) sigma = sigma*cos(3*Dist);
out = simspdata(1,sigma,rfm1,rfm2,par1,par2);
out=out$zdat; 
out=(rank(out)-.5)/length(out);
out=qnorm(out); 
out=matrix(out,d,d)
if(clsc==T) image.plot(z0,z0,out,xlab="",ylab="",xaxp=c(l.b,u.b,2),yaxp=c(l.b,u.b,2),col=heat.colors(12))
if(clsc==F) image(z0,z0,out,xlab="",ylab="",xaxp=c(l.b,u.b,2),yaxp=c(l.b,u.b,2),col=heat.colors(12))
out
}


#simulate data from the 1 factor spatial copula model
simspdata = function(n,sigma,rfm1,rfm2,par1,par2,df=300){
d = ncol(sigma);
if(df>=200) tem1 = rmvnorm(n,mean=rep(0,d),sigma=sigma);
if(df<200) tem1 = rmvt(n,sigma=sigma,df=df);
tem2 = rfm1(n,par1);
tem3 = rfm2(n,par2);
out = tem1 + tem2 - tem3;
udat = apply(out,2,rank);
udat = (udat-0.5)/n;
list(udat=udat,zdat=out)
}

#simulate data from the 1 factor spatial copula model+Student-t marginals
simspdata2marg = function(n,sigma,rfm1,rfm2,tmean,tsd,tdf,par1,par2,df=300){
d = ncol(sigma);
if(df>=200) tem1 = rmvnorm(n,mean=rep(0,d),sigma=sigma);
if(df<200) tem1 = rmvt(n,sigma=sigma,df=df);
tem2 = rfm1(n,par1);
tem3 = rfm2(n,par2);
out = tem1 + tem2 - tem3;
udat=exp1cdf(out,c(par1,par2));
zdat = tsd*qt(udat,tdf)+tmean;
#udat = apply(out,2,rank);
#udat = (udat-0.5)/n;
list(zdat=zdat,udat=udat)
}

#inverse margingal cdf 
exp1qcdf = function(q,par){
#qcdf=invfunc(q,exp1cdf,par,l.b=-10,u.b=10,tol=1e-12)
qcdf=invsec(q,exp1cdf,par,l.b=qnorm(q)-5e-1, u.b=qnorm(q),tol=1e-13)
qcdf
}

#doesn't work with adaptIntegrate......
#conditional cdf p|u 
exp1condcdf = function(p,u,pDist,par,nq=35){
Dist = pDist[-1,-1];
#cdf0=exp(-exp2lik0(u,Dist,par,N=6)); 
cdf0=exp(-llikf(u,Dist,par,grad=F)); 
#cdf1=function(x){out=exp2lik0(c(x,u),pDist,par,N=6); out= exp(-out); out}
cdf1=function(x){out=llikf(c(x,u),pDist,par,grad=F); out= exp(-out); out}
gl=gausslegendre(nq); wl=p*gl$weights; nl=p*gl$nodes;
cdf11=0;
for(iq in 1:nq) cdf11=cdf11+wl[iq]*cdf1(nl[iq]);
pcondcdf=cdf11/cdf0; 
pcondcdf
}


#predict u0 given u1,...,ud in the double exponential common factor model
#predicted quantiles given observed values Z(s1)=u1, ... , Z(sd)=ud 
exp1qcondcdf = function(q,u,pDist,par){
tem.f = function(p,par.f){return(exp1condcdf(p,u,pDist,par))}
#qcondcdf=invsec(q,tem.f,0,l.b=.3,u.b=0.7,tol=1e-6)
qcondcdf=invfunc(q,tem.f,0,l.b=1e-4,u.b=1-1e-4,tol=1e-6)
qcondcdf
}



#loglikelihood + derivatives wrt copula pars
#this is not optimal; replaced by llikf()
exp2lik = function(u,Dist,par,N){ 

if(min(Dist)<0){print("error in Dist: negative distances"); return(NULL)}
diag(Dist) = 0;

if(is.vector(u)) u = matrix(u,nrow=1)
n = nrow(u);
d = ncol(u);
th=par[3]; par=par[1:2];
sigma=exp(-th*Dist);
dsig1=-Dist*sigma;

siginv=solve(sigma);  
dets=det(sigma);    

llik = 0; lgrad = rep(0,2); lth = 0;
for(i in 1:n){ 
  z=exp1qcdf(u[i,],par); 
  cdf1  = exp1cdf(z,par,T);
  fpdf1 = fpdf(z,par,T);
  pdf1  = exp1pdf(z,par,T);
  lden  = sum(log(pdf1[,1])); 
  
  der1 = -cdf1[,2:3]/pdf1[,1];
  lden1 = pdf1[,4]*der1+pdf1[,2:3];
  lden1 = apply(lden1/pdf1[,1],2,sum);
  d.tem = function(x){        
    fV=fpdf(x,par,T);
    sigx=as.vector(siginv%*%(z-x)); 
    cf=(2*pi)^(d/2)*sqrt(dets);
    expd = exp(-t(z-x)%*%sigx/2);
    dmvn = expd/cf;
    out1 = dmvn*as.vector(fV[,1]);
    tem1 = as.vector(fV[,2:3])+apply(sigx*cdf1[,2:3]/pdf1[,1],2,sum)*as.vector(fV[,1]);
    out2 = tem1*dmvn; 
    out3 = .5*out1*as.vector(t(sigx)%*%dsig1%*%sigx); 
    c(out1,out2,out3)
  } 
  adint = adaptIntegrate(d.tem,lowerLimit=-N, upperLimit=N,fDim=4,tol=1e-5)
  intg = adint$integral; 
  llik = llik + lden -log(intg[1]); 
  lgrad = lgrad + lden1 - intg[2:3]/intg[1]; 
  lth = lth - intg[4]/intg[1]; 
}  
lth = lth + n*sum(diag(siginv%*%dsig1))/2; 
list(llik=llik,lgrad=c(lgrad,lth))
}



#exp2lik, likelihood only (for interpolation to get conditional u0|u1,...,ud)
exp2lik0 = function(u,Dist,par,N){ 

if(min(Dist)<0){print("error in Dist: negative distances"); return(NULL)}
diag(Dist) = 0;


if(is.vector(u)) u = matrix(u,nrow=1)
n = nrow(u);
d = ncol(u);
th=par[3]; par=par[1:2];
sigma=exp(-th*Dist);
dsig1=-Dist*sigma;

siginv=solve(sigma); 
dets=det(sigma);    

llik = 0; 
for(i in 1:n){ 
  z=exp1qcdf(u[i,],par); 
  pdf1  = exp1pdf(z,par,F);
  lden  = sum(log(pdf1)); 
  
  d.tem = function(x){        
    fV=fpdf(x,par,F);
    sigx=as.vector(siginv%*%(z-x)); 
    cf=(2*pi)^(d/2)*sqrt(dets);
    expd = exp(-t(z-x)%*%sigx/2);
    dmvn = expd/cf;
    out1 = dmvn*as.vector(fV);
    out1
  } 
  adint = adaptIntegrate(d.tem,lowerLimit=-N, upperLimit=N,fDim=1,tol=1e-5)
  intg = adint$integral; 
  llik = llik + lden -log(intg); 
}  
llik
}

#joint cdf in the double exponential common factor model
#very slow! doesn't work well with dim > 2
exp2cdf = function(u,sigma,par,N){
d = length(u);
z=pow1qcdf(u,par); 
f.tem = function(x){ return( pmvnorm(upper=z-x,mean=rep(0,d),sigma=sigma)*fpdf(x,par) ) }
adint = adaptIntegrate(f.tem,lowerLimit=-N, upperLimit=N,tol=1e-5)
lik=adint$integral
lik
}


#MLE in a selected common factor model; all parameters
#with gradient supplied, quite stable and fast. The algorithm usually reaches the maximum 
#sometimes fails because of poor estimation of the hessian matrix (iterates within tolerance).   
exp2mle = function(u,Dist,start,l.b=-Inf,u.b=Inf,N=12,pl=0){ 
llik=function(par){ 
  adj=0; 
  if(min(par - l.b) < 0 || max(par - u.b) > 0) adj = 1e5; 
  par[par<l.b] = l.b[par<l.b]; par[par>u.b]=u.b[par>u.b];
  out = exp2lik(u,Dist,par,N);  
  lik=out$llik; attr(lik,"gradient")=out$lgrad; lik+adj}
out = nlm(llik,p=start,print.level=pl,check.analyticals =F);
out
}

#fast mle with uniform ranks
#u is an n by d matrix of ranks or data with U(0,1) marginals
#n is number of replicates, d is number of locations
#Dist is a d by d matrix of distances
#start is vector of starting values (th1, th2, th, al)
#the model is E_U - E_L + Z(s)
#E_U ~ Exp(rate=th1), E_L ~ Exp(rate=th2)
#Z(s) has powered-exponential covariance function rho(d)=exp(-th*d^al)
exp2bmle = function(u,Dist,start,l.b=-Inf,u.b=Inf,pl=0){ 
llik=function(par){ 
  adj=0; 
  if(min(par - l.b) < 0 || max(par - u.b) > 0) { return(1e6) }
  #adj = 1e5; 
  #par[par<l.b] = l.b[par<l.b]; par[par>u.b]=u.b[par>u.b];
  #lik = llikf(u,Dist,par) + adj;  
  out = llikf(u,Dist,par,T);
  lik=out[1]; #attr(lik,"gradient")=out[2:4]; lik+adj}
  attr(lik,"gradient")=out[2:5]; lik+adj}
out = nlm(llik,p=start,print.level=pl,check.analyticals =F);
out
}

#slow full mle (no gradient) for model with t-marginals 
#(3 more pars to estimate: tmean, tsd, tdf) 
exp2cmle = function(tu,Dist,start,l.b=-Inf,u.b=Inf,pl=0){ 
llik=function(par){ 
  adj=0; 
  if(min(par - l.b) < 0 || max(par - u.b) > 0) { return(1e6) }
  tu0=(tu-par[1])/par[2]
  u=pt(tu0,df=par[3]);
  out = llikf(u,Dist,par[4:7],F); 
  lik=out[1]-sum(log(dt(tu0,df=par[3])))+length(tu)*log(par[2]); 
  lik }
out = nlm(llik,p=start,print.level=pl);
out
}

#Student-t density derivatives wrt z and df
dt1=function(z,df){
out=dt(z,df=df);
out1=-(1+1/df)*z*out/(1+z^2/df);
out2=(dt(z,df=df+1e-4)-dt(z,df=df-1e-4))/2e-4;
cbind(as.vector(out),as.vector(out1),as.vector(out2))
}

#fast full mle (with gradient) for model with t-marginals 
#(3 more pars to estimate: tmean, tsd, tdf)
#same as exp2bmle but we assume the process has t-marginals
#with mean tmean and sd=tsd, df = tdf
#assume the same marginal distribution for each location 
exp2dmle = function(tu,Dist,start,l.b=-Inf,u.b=Inf,pl=0){ 
llik=function(par){ 
  if(min(par - l.b) < 0 || max(par - u.b) > 0) { return(1e6) }
  tu0=(tu-par[1])/par[2]; 
  u=pt(tu0,df=par[3]);
  u1=(pt(tu0,df=par[3]+1e-4)-pt(tu0,df=par[3]-1e-4))/2e-4;
  grad0=dt1(tu0,df=par[3]);     
  du=grad0[,1]; 
  du1=grad0[,2]; ddf=grad0[,3];
  z=exp1qcdf(u,par[4:5]); zpdf=exp1pdf(z,par[4:5],gradth=F)
  out = llikf(z,Dist,par[4:7],T,zsc=T);  
  out1= llikf(z-1e-7*du/zpdf/par[2],Dist,par[4:7],F,zsc=T);
  out2= llikf(z-1e-7*du*tu0/zpdf/par[2],Dist,par[4:7],F,zsc=T);
  out3= llikf(z+1e-7*u1/zpdf,Dist,par[4:7],F,zsc=T);
  #lik=out[1]-sum(log(dt(tu0,df=par[3])))+length(tu)*log(par[2]); 
  lik=out[1]-sum(log(du))+length(tu)*log(par[2]); 
  lgrad=rep(0,7); lgrad[4:7]=out[2:5]; lgrad[1:3]=(c(out1,out2,out3)-out[1])/1e-7; 
  lgrad[1]=lgrad[1]+sum(du1/du)/par[2];
  lgrad[2]=lgrad[2]+sum(du1*tu0/du)/par[2]+length(tu)/par[2];
  lgrad[3]=lgrad[3]-sum(ddf/du);
  attr(lik,"gradient")=lgrad; lik  }
out = nlm(llik,p=start,print.level=pl,check.analyticals =F);
out
}