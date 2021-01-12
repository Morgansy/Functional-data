source("marginal_cpdf.R");
library(pbivnorm);
library(mvtnorm);
library(statmod);



#joint pdf for the model:
#W1=Z1+a1u*V0u-a1l*V0l+b1u*V1u-b1l*V1l (var1)
#W2=Z2+a2u*v0u-a2l*V0l                 (var2)
#w1,w2-vectors of length n (n locations)
#sigma is a (2*n) by (2*n) covariance matrix for Z=c(Z1,Z2)
#first/second n by n entries correspond to the correlation structure of Z1/Z2
#entries [1:n, (n+1):(2*d)] and [(n+1):(2*d), 1:n] - cross correlations
#the function is a one dimensional integral; but it is very fast
exppdf2var=function(w1,w2,a1u,a2u,b1u,b1l,a1l,a2l,sigma,N=11,nq=35){
#gl=gausslegendre(nq);
gl=gauss.quad.prob(nq);
nl=gl$nodes;
wl=gl$weights;
nl=N*nl;
wl=N*wl;
w=c(w1,w2);
d=length(w1);
invs=solve(sigma);
dets=det(sigma);
s11=sum(invs[1:d,1:d]);
s12=sum(invs[1:d,(d+1):(2*d)]);
s22=sum(invs[(d+1):(2*d),(d+1):(2*d)]);
sw=invs%*%w;
s1=sum(sw[1:d]);
s2=sum(sw[(d+1):(2*d)]);
c1=a1u*s1+a2u*s2-1;  
c2=-a1l*s1-a2l*s2-1; 
c11=a1u^2*s11+a2u^2*s22+2*a1u*a2u*s12;
c22=a1l^2*s11+a2l^2*s22+2*a1l*a2l*s12;
c12=a1l*a1u*s11+a2l*a2u*s22+(a1l*a2u+a1u*a2l)*s12;
dl=c11*c22-c12^2;  dl=max(dl,1e-10); #WARNING!
h1=-(c1*c22+c2*c12)/dl;
h2=-(c1*c12+c2*c11)/dl;
const1=-t(w)%*%invs%*%w-c1*h1-c2*h2;  
dlta=a1u*a2l-a1l*a2u; sdla=sign(dlta); 
kdlt=(a2u+a2l)/dlta;
sdlt=(s11*s22-s12^2);
cdlt=sdlt*dlta^2;
dconst=(2*pi)^(d-1)*sqrt(dl*dets);
pn1=h1*sqrt(dl/c22);
pn2=h2*sqrt(dl/c11);
rrho=c12/sqrt(c11*c22);  
out1a=pbivnorm(-pn1-nl*a2l*sdla*sqrt(sdlt/c22),-pn2-nl*a2u*sdla*sqrt(sdlt/c11),rrho);     
out2a=pbivnorm(-pn1+nl*a2l*sdla*sqrt(sdlt/c22),-pn2+nl*a2u*sdla*sqrt(sdlt/c11),rrho);  
if(min(out1a)< 1e-310 || min(out2a)<1e-310 ||is.finite(log(out1a))==FALSE || is.finite(log(out2a))==FALSE){ 
out1a[out1a<1e-310]=1e-310; out2a[out2a<1e-310]=1e-310; }

out1a=exp(log(out1a)+nl*kdlt-nl/b1u+const1/2-log(dconst));

out2a=exp(log(out2a)-nl*kdlt-nl/b1l+const1/2-log(dconst));  
if(is.finite(max(out1a))==FALSE || is.finite(max(out2a))==FALSE){ 
print("critical warning!!! small arguments rounded to zero in log() function");  return(NULL); }
out1=sum(wl*(out1a+out2a))/(b1u+b1l); 
out1
}




#joint copula density for the above model and the gradient (see exppdf2var())
#if abs(a1u/a2u-a1l/a2l) is very small, the function returns zero
#if the value is << 1e-300, the function returns zero
expcoppdf2var=function(u1,u2,a1u,a2u,b1u,b1l,a1l,a2l,cf1,cf2,sigma,N,nq,grad=F,
s0=NULL,s1=NULL,s2=NULL,sc1=NULL,sc2=NULL,sc0a=NULL,sc1a=NULL,sc2a=NULL,eps=1e-4,isw=F){

tol.lvl=0.08
#tol.lvl=0.02
if(abs(a1u-b1u)<1e-3) b1u=b1u+2e-3
if(abs(a1l-b1l)<1e-3) b1l=b1l+2e-3
d.check=a1u*a2l-a1l*a2u; 
 
if(d.check > 0 & d.check < tol.lvl) {add1=(tol.lvl-d.check)/(a1u+a2l);  a1u=a1u+add1; a2l=a2l+add1; } 
if(d.check <= 0 & d.check > -tol.lvl) {add1=(tol.lvl+d.check)/(a2u+a1l);  a1l=a1l+add1; a2u=a2u+add1;  } 

pr1=1/c(a1l,b1l,a1u,b1u);
pr2=1/c(a2u,a2l);
if(isw==F){ w1=exp2vqcdf(u1,pr1); w2=exp1qcdf(u2,pr2); }
if(isw==T){ w1=u1; w2=u2; grad=F; }
den=exppdf2var(w1,w2,a1u,a2u,b1u,b1l,a1l,a2l,sigma,N,nq); 
cpdf1=exp2vcdf(w1,1/a1l,1/b1l,1/a1u,1/b1u,der=grad);
cpdf2=exp1pdf(w2,1/c(a2u,a2l),gradth=grad);
cdf2=exp1cdf(w2,1/c(a2u,a2l),gradth=grad);
if(grad==F) cp12=prod(cpdf1*cpdf2);
if(grad==T) cp12=prod(cpdf1[,10]*cpdf2[,1]); 
copden=den/cp12;
if(grad==T){
da1u=cpdf1[,10]*a1u^2;
db1u=cpdf1[,10]*b1u^2;
da1l=cpdf1[,10]*a1l^2;
db1l=cpdf1[,10]*b1l^2;
da2u=cpdf2[,1]*a2u^2;
da2l=cpdf2[,1]*a2l^2;
inv.a1u=cpdf1[,2]/da1u;
inv.b1u=cpdf1[,3]/db1u;
inv.a1l=cpdf1[,4]/da1l;
inv.b1l=cpdf1[,5]/db1l;
inv.a2u=cdf2[,2]/da2u;
inv.a2l=cdf2[,3]/da2l;
w1.a1u=w1+eps*inv.a1u;
w1.b1u=w1+eps*inv.b1u;
w1.a1l=w1+eps*inv.a1l;
w1.b1l=w1+eps*inv.b1l;
w2.a2u=w2+eps*inv.a2u;
w2.a2l=w2+eps*inv.a2l;


den0=exppdf2var(w1,w2,a1u,a2u,b1u,b1l,a1l,a2l,s0,N,nq); 
den1=exppdf2var(w1,w2,a1u,a2u,b1u,b1l,a1l,a2l,s1,N,nq); 
den2=exppdf2var(w1,w2,a1u,a2u,b1u,b1l,a1l,a2l,s2,N,nq); 
denc1=exppdf2var(w1,w2,a1u,a2u,b1u,b1l,a1l,a2l,sc1,N,nq); 
denc2=exppdf2var(w1,w2,a1u,a2u,b1u,b1l,a1l,a2l,sc2,N,nq); 
denc0a=exppdf2var(w1,w2,a1u,a2u,b1u,b1l,a1l,a2l,sc0a,N,nq); 
denc1a=exppdf2var(w1,w2,a1u,a2u,b1u,b1l,a1l,a2l,sc1a,N,nq); 
denc2a=exppdf2var(w1,w2,a1u,a2u,b1u,b1l,a1l,a2l,sc2a,N,nq); 


den.0=(den0-den)/eps; 
den.1=(den1-den)/eps;
den.2=(den2-den)/eps;
den.c1=(denc1-den)/eps; 
den.c2=(denc2-den)/eps;
den.c0a=(denc0a-den)/eps; 
den.c1a=(denc1a-den)/eps; 
den.c2a=(denc2a-den)/eps;


copden.a1u=exppdf2var(w1.a1u,w2,a1u+eps,a2u,b1u,b1l,a1l,a2l,sigma,N,nq);
copden.a2u=exppdf2var(w1,w2.a2u,a1u,a2u+eps,b1u,b1l,a1l,a2l,sigma,N,nq); 
copden.b1u=exppdf2var(w1.b1u,w2,a1u,a2u,b1u+eps,b1l,a1l,a2l,sigma,N,nq); 
copden.b1l=exppdf2var(w1.b1l,w2,a1u,a2u,b1u,b1l+eps,a1l,a2l,sigma,N,nq); 
copden.a1l=exppdf2var(w1.a1l,w2,a1u,a2u,b1u,b1l,a1l+eps,a2l,sigma,N,nq); 
copden.a2l=exppdf2var(w1,w2.a2l,a1u,a2u,b1u,b1l,a1l,a2l+eps,sigma,N,nq); 

den.1.a1l=(log(copden.a1l)-log(den))/eps; 
den.1.b1l=(log(copden.b1l)-log(den))/eps;
den.1.a1u=(log(copden.a1u)-log(den))/eps;
den.1.b1u=(log(copden.b1u)-log(den))/eps;
den.1.a2u=(log(copden.a2u)-log(den))/eps;
den.1.a2l=(log(copden.a2l)-log(den))/eps;

ff1=cpdf1[,11]/cpdf1[,10];
ff2=cpdf2[,4]/cpdf2[,1];
der.a1u = den.1.a1u+sum(cpdf1[,6]/da1u-ff1*inv.a1u);  
der.a2u = den.1.a2u+sum(cpdf2[,2]/da2u-ff2*inv.a2u);   
der.b1u = den.1.b1u+sum(cpdf1[,7]/db1u-ff1*inv.b1u);
der.b1l = den.1.b1l+sum(cpdf1[,9]/db1l-ff1*inv.b1l);
der.a1l = den.1.a1l+sum(cpdf1[,8]/da1l-ff1*inv.a1l);
der.a2l = den.1.a2l+sum(cpdf2[,3]/da2l-ff2*inv.a2l); 
copden=c(copden,der.a1u,der.a2u,der.b1u,der.b1l,der.a1l,der.a2l,den.0/cp12,
den.1/cp12,den.2/cp12,den.c1/cp12,den.c2/cp12,den.c0a/cp12,den.c1a/cp12,den.c2a/cp12);
}
copden
}

#negative loglikelihood function and gradient (if grad=T) for LMC copula model
#eps - difference used to compute derivatives numerically
#N - upper bound used for numerical integration on the real line (N=12 works well in most cases)
#nq - number of quadrarue points for numerical integration (nq=21 gives accurate results)
#this is a faster method for computing gradient than a naive approach  
exp2varlik=function(umat1,umat2,par,Dist,N,nq,grad=T,eps=1e-7,isw=F){
if(min(Dist)<0){print("error in Dist: negative distances"); return(NULL)}
diag(Dist) = 0;

if(is.vector(umat1)) umat1 = matrix(umat1,nrow=1)
if(is.vector(umat2)) umat2 = matrix(umat2,nrow=1)
n = nrow(umat1);
d = ncol(umat1);
dd= ncol(Dist);
n2 = nrow(umat2);
d2 = ncol(umat2);
if(n!=n2 ||d!=d2) {print("dimensions of umat1 and umat2 don't match"); return(NULL)}
if(dd!=d) {print("umat and Dist are not compatible"); return(NULL)}

a1u=par[1];      a2u=par[2];      b1u=par[3];      b1l=par[4];      a1l=par[5]; 
a2l=par[6];      th0=par[7];      th1=par[8];      th2=par[9];  
cf1=par[10];     cf2=par[11];     cf0a=par[12];    cf1a=par[13];    cf2a=par[14];


cr0=exp(-th0*Dist^cf0a);
cr1=exp(-th1*Dist^cf1a);
cr2=exp(-th2*Dist^cf2a);
ii1=(1:d); ii2=((d+1):(2*d));

sigma=matrix(0,ncol=(2*d),nrow=(2*d));
sigma[ii1,ii1]=cf1*cf1*cr0+(1-cf1*cf1)*cr1;
sigma[ii1,ii2]=cf1*cf2*cr0;
sigma[ii2,ii1]=cf1*cf2*cr0;
sigma[ii2,ii2]=cf2*cf2*cr0+(1-cf2*cf2)*cr2;

#we calculate derivatives numerically
#we avoid multiple computations of the inverse cdfs to make the code faster
if(grad==T){

th0e=th0+eps; th1e=th1+eps; th2e=th2+eps;
cf0ae=cf0a+eps; cf1ae=cf1a+eps; cf2ae=cf2a+eps;

cr0e=exp(-th0e*Dist^cf0a);
cr1e=exp(-th1e*Dist^cf1a);
cr2e=exp(-th2e*Dist^cf2a);
cr0ae=exp(-th0*Dist^cf0ae);
cr1ae=exp(-th1*Dist^cf1ae);
cr2ae=exp(-th2*Dist^cf2ae);
s0=matrix(0,ncol=(2*d),nrow=(2*d));
s1=sigma;
s2=sigma;
sc1=sigma;
sc2=sigma;
sc0a=sigma;
sc1a=sigma;
sc2a=sigma;
cf1e=cf1+eps;
cf2e=cf2+eps;

s0[ii1,ii1]=cf1*cf1*cr0e+(1-cf1*cf1)*cr1;
s1[ii1,ii1]=cf1*cf1*cr0+(1-cf1*cf1)*cr1e;
sc1[ii1,ii1]=cf1e*cf1e*cr0+(1-cf1e*cf1e)*cr1;
sc1a[ii1,ii1]=cf1*cf1*cr0+(1-cf1*cf1)*cr1ae;
sc0a[ii1,ii1]=cf1*cf1*cr0ae+(1-cf1*cf1)*cr1;
s0[ii1,ii2]=cf1*cf2*cr0e;
sc0a[ii1,ii2]=cf1*cf2*cr0ae;
sc1[ii1,ii2]=cf1e*cf2*cr0;
sc2[ii1,ii2]=cf2e*cf1*cr0;
s0[ii2,ii1]=cf1*cf2*cr0e;
sc0a[ii2,ii1]=cf1*cf2*cr0ae;
sc1[ii2,ii1]=cf1e*cf2*cr0;
sc2[ii2,ii1]=cf2e*cf1*cr0;
s0[ii2,ii2]=cf2*cf2*cr0e+(1-cf2*cf2)*cr2;
s2[ii2,ii2]=cf2*cf2*cr0+(1-cf2*cf2)*cr2e;
sc2[ii2,ii2]=cf2e*cf2e*cr0+(1-cf2e*cf2e)*cr2;
sc2a[ii2,ii2]=cf2*cf2*cr0+(1-cf2*cf2)*cr2ae;
sc0a[ii2,ii2]=cf2*cf2*cr0ae+(1-cf2*cf2)*cr2;
}


llik = 0; 
if(grad==T) llik=rep(0,15);
for(i in 1:n){  
i.llik=expcoppdf2var(umat1[i,],umat2[i,],a1u,a2u,b1u,b1l,a1l,a2l,cf1,cf2,sigma,N,nq,grad,s0,s1,s2,sc1,sc2,sc0a,sc1a,sc2a,eps,isw=isw); 
if(abs(i.llik[1])<1e-305) i.llik[1]=1e-305; 
if(grad==F) llik = llik - log(i.llik); 
if(grad==T) {llik[1] = llik[1] - log(i.llik[1]); llik[2:7] = llik[2:7] - i.llik[2:7]; llik[8:15] = llik[8:15] - i.llik[8:15]/i.llik[1];}
}
llik
}


#vector of lower bounds
lb.EXP=c(rep(0.05,6),rep(0.01,3),rep(0.01,2),rep(0.01,3)); 
#vector of upper bounds
ub.EXP=c(rep(2,6),rep(5,3),rep(.999,2),rep(2,3));
#vector of starting values
st.EXP=c(1.5,1,1,1,1,1,0.5,0.5,0.5,0.5,0.5,1,1,1)
#MLE for LMC copula model
exp3varmle = function(umat1,umat2,Dist,start=st.EXP,l.b=lb.EXP,u.b=ub.EXP,pl=2,N=12,nq=21,i.fix=NULL){ 
llik=function(par){ 
  if(min((par - l.b)) < 0 || max((par - u.b)) > 0) return(1e6) 
  
  out = exp2varlik(umat1,umat2,par,Dist,N,nq,grad=T,eps=1e-7); 
  lik = out[1]; grout=out[2:15]; grout[i.fix]=0; attr(lik,"gradient")=grout; 
  lik}
out = nlm(llik,p=start,print.level=pl,check.analyticals =F);
out
}

#conditional cdf p|u 
exp3condcdf = function(p,u1,u2,pDist,par,nq=35,N=12,isw=F){
Dist = pDist[-1,-1]; 
pr1=1/par[c(5,4,1,3)]; pr2=1/par[c(2,6)];
if(isw==F) {w1=exp2vqcdf(u1,pr1); w2=exp1qcdf(u2,pr2);}
if(isw==T) {w1=u1; w2=u2;}
cdf0=exp(-exp2varlik(w1,w2,par,Dist,N,nq=21,grad=F,eps=1e-7,isw=T))
cdf1=function(x){ x1=exp2vqcdf(x[1],pr1);  x2=exp1qcdf(x[2],pr2); out=exp2varlik(c(x1,w1),c(x2,w2),par,pDist,N,nq=21,grad=F,eps=1e-7,isw=T); out= exp(-out); out}
p1=p[1]; p2=p[2];
gl=gauss.quad.prob(nq); wl1=p1*gl$weights; nl1=p1*gl$nodes; wl2=p2*gl$weights; nl2=p2*gl$nodes;
cdf11=0; 
for(iq1 in 1:nq) {
for(iq2 in 1:nq) {cdf11=cdf11+wl1[iq1]*wl2[iq2]*cdf1(c(nl1[iq1],nl2[iq2])); }}
pcondcdf=cdf11/cdf0; 
pcondcdf=max(1.1e-5,pcondcdf);
pcondcdf=min(1-1.1e-5,pcondcdf);
pcondcdf
}


#conditional cdf p|u 
mvn3condcdf = function(p,u1,u2,pDist,par,nq=25){
Dist = pDist[-1,-1]; 
cdf0=exp(-mvn3varlik(u1,u2,par,Dist))
cdf1=function(u){ out=mvn3varlik(c(u[1],u1),c(u[2],u2),par,pDist); out= exp(-out); out}
p1=p[1]; p2=p[2];
gl=gauss.quad.prob(nq); wl1=p1*gl$weights; nl1=p1*gl$nodes; wl2=p2*gl$weights; nl2=p2*gl$nodes;
cdf11=0; 
for(iq1 in 1:nq) {
for(iq2 in 1:nq) {cdf11=cdf11+wl1[iq1]*wl2[iq2]*cdf1(c(nl1[iq1],nl2[iq2])); }}
pcondcdf=cdf11/cdf0; 
pcondcdf=max(1.1e-4,pcondcdf);
pcondcdf=min(1-1.1e-4,pcondcdf);
pcondcdf
}


#conditional q-quantile given u1,...,ud in LMC copula model
exp3qcondcdf = function(q,u1,u2,pDist,par,nq=21){
pr1=1/par[c(5,4,1,3)]; pr2=1/par[c(2,6)];
w1=exp2vqcdf(u1,pr1); w2=exp1qcdf(u2,pr2);
tem.f1 = function(p,par.f){return(exp3condcdf(c(1,p),w1,w2,pDist,par,nq,isw=T))}
tem.f2 = function(p,par.f){return(exp3condcdf(c(p,1),w1,w2,pDist,par,nq,isw=T))}
qcondcdf1=invfunc(q[1],tem.f1,0,l.b=1e-4,u.b=1-1e-4,tol=1e-4)
qcondcdf2=invfunc(q[2],tem.f2,0,l.b=1e-4,u.b=1-1e-4,tol=1e-4)
c(qcondcdf2,qcondcdf1)
}


#conditional q-quantile given u1,...,ud in bivariate normal LMC model
mvn3qcondcdf = function(q,u1,u2,pDist,par,nq=21){
tem.f1 = function(p,par.f){return(mvn3condcdf(c(1,p),u1,u2,pDist,par,nq))}
tem.f2 = function(p,par.f){return(mvn3condcdf(c(p,1),u1,u2,pDist,par,nq))}
qcondcdf1=invfunc(q[1],tem.f1,0,l.b=1e-3,u.b=1-1e-3,tol=1e-4)
qcondcdf2=invfunc(q[2],tem.f2,0,l.b=1e-3,u.b=1-1e-3,tol=1e-4)
c(qcondcdf2,qcondcdf1)
}

##################################################################################################################################
#simulate data from a bivariate LMC copula model
#df < 100 - bivariate Student-t model + common factors
#df > 100 - bivariate normal model + common factors
#pr = T - Pareto common factors
#pr = F - exponential common factors
simexp3var = function(N,Dist,par,b2u=0,b2l=0,df=200,pr=F){

a1u=par[1];      a2u=par[2];      b1u=par[3];      b1l=par[4];      a1l=par[5]; 
a2l=par[6];      th0=par[7];      th1=par[8];      th2=par[9];  
cf1=par[10];     cf2=par[11];     cf0a = par[12];  cf1a=par[13];    cf2a=par[14];


cr0=exp(-th0*Dist^cf0a);
cr1=exp(-th1*Dist^cf1a);
cr2=exp(-th2*Dist^cf2a);
d=ncol(Dist);
if(pr==F) {E0u=rexp(N,1); E0l=rexp(N,1); E1u=rexp(N,1); E1l=rexp(N,1); E2u=rexp(N,1); E2l=rexp(N,1);}
if(pr==T) {E0u=rpow(N,c(1,4)); E0l=rpow(N,c(1,4)); E1u=rpow(N,c(1,4)); E1l=rpow(N,c(1,4)); E2u=rpow(N,c(1,4)); E2l=rpow(N,c(1,4));}

if(df > 100){
z0=rmvnorm(N,mean=rep(0,d),sigma=cr0);
z1=rmvnorm(N,mean=rep(0,d),sigma=cr1);
z2=rmvnorm(N,mean=rep(0,d),sigma=cr2); }
if(df <= 100){
z0=rmvt(N,sigma=cr0,df=df);
z1=rmvt(N,sigma=cr1,df=df);
z2=rmvt(N,sigma=cr2,df=df); }

zz1=cf1*z0+sqrt(1-cf1*cf1)*z1;
zz2=cf2*z0+sqrt(1-cf2*cf2)*z2;

W1=zz1+a1u*E0u+b1u*E1u-a1l*E0l-b1l*E1l;
W2=zz2+a2u*E0u+b2u*E2u-a2l*E0l-b2l*E2l;

umat1=(apply(W1,2,rank)-.5)/N
umat2=(apply(W2,2,rank)-.5)/N

list(W1=W1,W2=W2,umat1=umat1,umat2=umat2)
}

##################################################################################################################################


#bivariate normal LMC model; negative log-likelihood for the joint normal distribution 
LMCvarlik=function(zmat1,zmat2,par,Dist){
if(min(Dist)<0){print("error in Dist: negative distances"); return(NULL)}
diag(Dist) = 0;

if(is.vector(zmat1)) zmat1 = matrix(zmat1,nrow=1)
if(is.vector(zmat2)) zmat2 = matrix(zmat2,nrow=1)
n = nrow(zmat1);
d = ncol(zmat1);
dd= ncol(Dist);
n2 = nrow(zmat2);
d2 = ncol(zmat2);
if(n!=n2 ||d!=d2) {print("dimensions of umat1 and umat2 don't match"); return(NULL)}
if(dd!=d) {print("umat and Dist are not compatible"); return(NULL)}

th0=par[1]; th1=par[2]; th2=par[3];
pw0=par[6]; pw1=par[7]; pw2=par[8];
cf1=par[4]; cf2=par[5];


cr0=exp(-th0*Dist^pw0);
cr1=exp(-th1*Dist^pw1);
cr2=exp(-th2*Dist^pw2);
ii1=(1:d); ii2=((d+1):(2*d));

sigma=matrix(0,ncol=(2*d),nrow=(2*d));
sigma[ii1,ii1]=cf1*cf1*cr0+(1-cf1*cf1)*cr1;
sigma[ii1,ii2]=cf1*cf2*cr0;
sigma[ii2,ii1]=cf1*cf2*cr0;
sigma[ii2,ii2]=cf2*cf2*cr0+(1-cf2*cf2)*cr2;


zmat12=cbind(zmat1,zmat2); 

llik = 0; 
invs=solve(sigma); dets=det(sigma); 
for(i in 1:n){  
zmat12i=zmat12[i,]; llik=llik + .5*zmat12i%*%invs%*%zmat12i;
}
llik = llik + .5*n*log(dets) + n*d*log(2*pi);
llik
}

# bivariate normal LMC model; negative log-likelihood for the joint normal copula
mvn3varlik=function(umat1,umat2,par,Dist){
zmat1=qnorm(umat1); zmat2=qnorm(umat2);
llik=LMCvarlik(zmat1,zmat2,par,Dist)+sum(log(dnorm(zmat1)))+sum(log(dnorm(zmat2)));
llik
}

#vector of starting values
st.LMC=c(0.90,0.90,0.90,0.50,0.50,0.30,0.30,0.30); 
#vector of lower bounds
lb.LMC=c(0.01,0.01,0.01,0.000,0.000,0.001,0.001,0.001);
#vector of upper bounds
ub.LMC=c(5.00,5.00,5.00,1.000,1.000,1.999,1.999,1.999);
#bivariate normal LMC model MLE
LMCvarmle1 = function(umat1,umat2,Dist,start=st.LMC,l.b=lb.LMC,u.b=ub.LMC,pl=2){ 
llik=function(par){ 
  ee=0.001;
  ilb=(par-l.b<  ee & par > l.b);  iub=(par-u.b> -ee & par < u.b);
  if(min((par - l.b)[!(ilb | iub)]) < 0 || max((par - u.b)[!(ilb | iub)]) > 0) { return(1e6); } 
  par[par<l.b] = l.b[par<l.b]; par[par>u.b]=u.b[par>u.b]; 
  
  print(par); 
  lik = mvn3varlik(umat1,umat2,par,Dist); 
  emat=diag(rep(1e-4,8));
  glik=rep(0,8);
  for(i in 1:8){
  gg1=mvn3varlik(umat1,umat2,par+emat[i,],Dist);
  gg2=mvn3varlik(umat1,umat2,par-emat[i,],Dist);
  glik[i]=(gg1-gg2)/2e-4
  }
  glik[(ilb | iub)]=0; attr(lik,"gradient")=glik;
  lik}
out = nlm(llik,p=start,print.level=pl,iterlim = 300, check.analyticals=F);
out
}

#bivariate normal LMC model MLE
#same as LMCvarmle() 
#if some parameters reach boundaries (l.b or u.b) it restarts LMCvarmle() with these parameters fixed
LMCvarmle = function(zmat1,zmat2,Dist,start=st.LMC,l.b=lb.LMC,u.b=ub.LMC,pl=1){
  out=LMCvarmle1(zmat1,zmat2,Dist,start=st.LMC,l.b=lb.LMC,u.b=ub.LMC,pl=pl);
  if(mean(abs(out$gradient)) > 5){ cat("one of the covariance matrix parameters reached boundaries\n\n");
  out=LMCvarmle1(zmat1,zmat2,Dist,start=out$estimate,l.b=lb.LMC,u.b=ub.LMC,pl=pl);}
  if(mean(abs(out$gradient)) > 5){ cat("two of the covariance matrix parameters reached boundaries\n\n");
  out=LMCvarmle1(zmat1,zmat2,Dist,start=out$estimate,l.b=lb.LMC,u.b=ub.LMC,pl=pl);}
  out
}

