#function inverse
#bisection method
#func() - function of one argument to be inverted
#th - function parameter
#q - argument of the inverse function (can be a vector if func() supports a vector argument)
invfunc = function(q,func,th,l.b=1e-12,u.b=1-1e-12,tol=1e-8){ 
lq = length(q);
v0 = rep(l.b,lq);
v1 = rep(u.b,lq);
out = rep(0, lq);
iconv = rep(F, lq); 
func0 = function(v,q,th){return(func(v,th)-q)};
f0 = func0(v0,q,th); f1 = func0(v1,q,th);
l0 = (f0 > 0 & f1 > 0);
l1 = (f0 < 0 & f1 < 0);
iconv[l0 | l1] = T; 
tol0 = rep(1,lq); tol0[iconv] = 0; 
while(prod(iconv) < 0.5) {  
  v0a = v0[!iconv]; v1a = v1[!iconv]; 
  v2a = (v0a+v1a)/2; qa = q[!iconv]; 
  f0a = f0[!iconv]; f1a = f1[!iconv];
  f2a = func0(v2a,qa,th);                      
  i1 = (f0a*f2a < 0); i2 = (f1a*f2a < 0);
  v1a[i1] = v2a[i1]; v1[!iconv] = v1a; 
  v0a[i2] = v2a[i2]; v0[!iconv] = v0a;
  
  f0[!iconv][i1] = f0a[i1];
  f0[!iconv][i2] = f2a[i2];
  f1[!iconv][i1] = f2a[i1];
  f1[!iconv][i2] = f1a[i2];
  
  #tol0[!iconv] = abs(func0((v0a+v1a)/2,qa,th)); 
  tol0[!iconv] = abs(f2a); 
  iconv = (abs(tol0) < tol)     
}

out=(v0+v1)/2
out[l0] = l.b; out[l1] = u.b
out
}

#function inverse
#secant method
#func() - function of one argument to be inverted
#th - function parameter
#q - argument of the inverse function (can be a vector if func() supports a vector argument)
invsec = function(q,func,th,l.b=1e-12,u.b=1-1e-12,tol=1e-8){
lq = length(q);
v0 = l.b; v1 = u.b;
if(length(l.b)==1) v0=rep(v0,lq);
if(length(u.b)==1) v1=rep(v1,lq);
out = rep(0, lq);
iconv = rep(F, lq);
func0 = function(v,q,th){return(func(v,th)-q)};
f0 = func0(v0,q,th); f1 = func0(v1,q,th);
tol0 = rep(1,lq); tol0[iconv] = 0;  
while(prod(iconv) < 0.5) {  
  v0a = v0[!iconv]; v1a = v1[!iconv]; 
  qa = q[!iconv]; 
  f0a = f0[!iconv]; f1a = f1[!iconv]; 
  v2a = (v0a*f1a-v1a*f0a)/(f1a-f0a);  
  f2a = func0(v2a,qa,th); 
  v0a = v1a; v1a = v2a;         
  v1[!iconv] = v1a;  v0[!iconv] = v0a;

  f0[!iconv] = f1a;
  f1[!iconv] = f2a;
  
  tol0[!iconv] = abs(f2a); 
  iconv = (abs(tol0) < tol)
}

out=(v0+v1)/2
out
}


#Pareto cdf
powcdf=function(z,par){
th=par[1]; b=par[2];
if(th<0 || b<0) {print("error: negative par"); return(NULL)}
cdf=(z>th)*(1-(abs(z/th))^(-b));
cdf
}

#Pareto inverse cdf
powinvcdf=function(q,par){
lq=(1-q)^(-1/par[2]);
qcdf=par[1]*lq;
qcdf
}

#Pareto radnom number generation
rpow = function(n,par){
u = runif(n); 
z = powinvcdf(u,par);
z
}



#marginal cdf+ derivatives wrt pars+ wrt argument
#warning! pars are reciprocals of the scale parameters
#th1l=1/a1l, etc.
#if der==F, only pdf is calculated
#if der==T, cdf, its derivatives wrt pars + pdf (with no derivatives)
exp2vcdf = function(z,th1l,th2l,th1u,th2u,der=F){

dz=dnorm(z);
p1u=pnorm(z-th1u);
p2u=pnorm(z-th2u);
p1l=pnorm(-z-th1l);
p2l=pnorm(-z-th2l);
ex1u=exp(th1u^2/2-th1u*z);
ex2u=exp(th2u^2/2-th2u*z);
ex1l=exp(th1l^2/2+th1l*z);
ex2l=exp(th2l^2/2+th2l*z);

e1u=ex1u*p1u;
e2u=ex2u*p2u;
e1l=ex1l*p1l;
e2l=ex2l*p2l;

ee1u=e1u*(th1u-z)-dz;
ee2u=e2u*(th2u-z)-dz;
ee1l=e1l*(th1l+z)-dz;
ee2l=e2l*(th2l+z)-dz;


t1u=th1l*th2l*th2u/(th1l+th1u)/(th2l+th1u)/(th1u-th2u);
t2u=th1l*th2l*th1u/(th2l+th2u)/(th1l+th2u)/(th1u-th2u);
t1l=th2l*th1u*th2u/(th1l+th1u)/(th1l+th2u)/(th1l-th2l);
t2l=th1l*th1u*th2u/(th2l+th2u)/(th2l+th1u)/(th1l-th2l);

outz=-t1u*e1u*th1u+t2u*e2u*th2u-(t1l*e1l*th1l-t2l*e2l*th2l);
outz=outz+dz+(t1u-t2u+t1l-t2l)*dz;

if(der==T){

#let's try to add deriv-s of pdf wrt par-s
e1uz=dz-th1u*e1u;
e2uz=dz-th2u*e2u;
e1lz=-dz+th1l*e1l;
e2lz=-dz+th2l*e2l;

outzz=-t1u*e1uz*th1u+t2u*e2uz*th2u-(t1l*e1lz*th1l-t2l*e2lz*th2l)-z*dz*(1+t1u-t2u+t1l-t2l);

#zdz=z*dz;
#ee1uz=-e1u+zdz+(dz-th1u*e1u)*(th1u-z);
#ee2uz=-e2u+zdz+(dz-th2u*e2u)*(th2u-z);
#ee1lz= e1l+zdz+(-dz+th1l*e1l)*(th1l+z);
#ee2lz= e2l+zdz+(-dz+th2l*e2l)*(th2l+z);
ee1uz=-e1u*(1+th1u^2)+dz*th1u+z*th1u*e1u;
ee2uz=-e2u*(1+th2u^2)+dz*th2u+z*th2u*e2u;
ee1lz= e1l*(1+th1l^2)-dz*th1l+z*th1l*e1l;
ee2lz= e2l*(1+th2l^2)-dz*th2l+z*th2l*e2l;

t1u1u=-t1u*(1/(th1l+th1u)+1/(th2l+th1u)+1/(th1u-th2u));
t1u2u=t1u*th1u/(th1u-th2u)/th2u;
t1u1l=t1u*th1u/(th1l+th1u)/th1l;
t1u2l=t1u*th1u/(th2l+th1u)/th2l;

t2u1u=-t2u*th2u/(th1u-th2u)/th1u;
t2u2u=-t2u*(1/(th2l+th2u)+1/(th1l+th2u)-1/(th1u-th2u));
t2u1l=t2u*th2u/(th1l+th2u)/th1l;
t2u2l=t2u*th2u/(th2l+th2u)/th2l;

t1l1u=t1l*th1l/(th1l+th1u)/th1u;
t1l2u=t1l*th1l/(th1l+th2u)/th2u;
t1l1l=-t1l*(1/(th1l+th1u)+1/(th1l+th2u)+1/(th1l-th2l));
t1l2l=t1l*th1l/(th1l-th2l)/th2l;

t2l1u=t2l*th2l/(th2l+th1u)/th1u;
t2l2u=t2l*th2l/(th2l+th2u)/th2u;
t2l1l=-t2l*th2l/(th1l-th2l)/th1l;
t2l2l=-t2l*(1/(th2l+th2u)+1/(th2l+th1u)-1/(th1l-th2l));

out=pnorm(z)+t1u*e1u-t2u*e2u-(t1l*e1l-t2l*e2l);
out1u=t1u1u*e1u+t1u*ee1u-t2u1u*e2u-(t1l1u*e1l-t2l1u*e2l);
out2u=t1u2u*e1u-(t2u2u*e2u+t2u*ee2u)-(t1l2u*e1l-t2l2u*e2l);
out1l=t1u1l*e1u-t2u1l*e2u-(t1l1l*e1l+t1l*ee1l-t2l1l*e2l);
out2l=t1u2l*e1u-t2u2l*e2u-(t1l2l*e1l-(t2l2l*e2l+t2l*ee2l));

out1uz=t1u1u*e1uz+t1u*ee1uz-t2u1u*e2uz-(t1l1u*e1lz-t2l1u*e2lz);
out2uz=t1u2u*e1uz-(t2u2u*e2uz+t2u*ee2uz)-(t1l2u*e1lz-t2l2u*e2lz);
out1lz=t1u1l*e1uz-t2u1l*e2uz-(t1l1l*e1lz+t1l*ee1lz-t2l1l*e2lz);
out2lz=t1u2l*e1uz-t2u2l*e2uz-(t1l2l*e1lz-(t2l2l*e2lz+t2l*ee2lz));

outz = cbind(out,out1u,out2u,out1l,out2l,out1uz,out2uz,out1lz,out2lz,outz,outzz)
}
outz
}

#marginal cdf without derivatives
#used to find inverse marginal cdf
exp2vcdf0 = function(z,par){

th1l=par[1]; th2l=par[2]; th1u=par[3]; th2u=par[4];
dz=dnorm(z);
p1u=pnorm(z-th1u);
p2u=pnorm(z-th2u);
p1l=pnorm(-z-th1l);
p2l=pnorm(-z-th2l);
ex1u=exp(th1u^2/2-th1u*z);
ex2u=exp(th2u^2/2-th2u*z);
ex1l=exp(th1l^2/2+th1l*z);
ex2l=exp(th2l^2/2+th2l*z);

e1u=ex1u*p1u;
e2u=ex2u*p2u;
e1l=ex1l*p1l;
e2l=ex2l*p2l;

t1u=th1l*th2l*th2u/(th1l+th1u)/(th2l+th1u)/(th1u-th2u);
t2u=th1l*th2l*th1u/(th2l+th2u)/(th1l+th2u)/(th1u-th2u);
t1l=th2l*th1u*th2u/(th1l+th1u)/(th1l+th2u)/(th1l-th2l);
t2l=th1l*th1u*th2u/(th2l+th2u)/(th2l+th1u)/(th1l-th2l);

out=pnorm(z)+t1u*e1u-t2u*e2u-(t1l*e1l-t2l*e2l);

out
}

#inverse margingal cdf 
exp2vqcdf = function(q,par){
#qcdf=invfunc(q,exp2vcdf0,par,l.b=-10,u.b=10,tol=1e-12)
qcdf=invsec(q,exp2vcdf0,par,l.b=qnorm(q)-5e-1, u.b=qnorm(q),tol=1e-10) 
qcdf
}

#marginal cdf and its derivatives wrt th1,th2
exp1cdf = function(z,par,gradth=F){
th1=par[1]; th2=par[2];
t1=pnorm(z-th1);
e1=exp(th1^2/2-th1*z);
t2=pnorm(-z-th2);
e2=exp(th2^2/2+th2*z);
te1=t1*e1;
te2=t2*e2;
m1=te1*th2;
m2=te2*th1;
cdf=pnorm(z)-(m1-m2)/(th1+th2);
#grad.th=NULL;
if(gradth==T){
tem1=-(th1-z)*te1+dnorm(z)+(te1+te2)/(th1+th2);
tem2= (th2+z)*te2-dnorm(z)-(te1+te2)/(th1+th2);
grad.th1=tem1*th2/(th1+th2);
grad.th2=tem2*th1/(th1+th2);
grad.th=cbind(grad.th1,grad.th2);
cdf=cbind(cdf,grad.th)
}
cdf
}

#marginal pdf and its derivatives wrt arg, th1,th2
exp1pdf = function(z,par,gradth=F){
th1=par[1]; th2=par[2];
t1=pnorm(z-th1);
e1=exp(th1^2/2-th1*z);
t2=pnorm(-z-th2);
e2=exp(th2^2/2+th2*z);
te1=t1*e1;
te2=t2*e2;
pdf=th1*th2*(te1+te2)/(th1+th2);
#grad.th=NULL;
if(gradth==T){
tem1a = (th1-z)*te1-dnorm(z);
tem1b = (z+th2)*te2-dnorm(z);
tem0 = (te1+te2)/(th1+th2)^2;
grad.th1 = th2^2*tem0+th1*th2*tem1a/(th1+th2);
grad.th2 = th1^2*tem0+th1*th2*tem1b/(th1+th2);
grad.u = (th2*te2-th1*te1)*th1*th2/(th1+th2);
grad.th=cbind(grad.th1,grad.th2,grad.u);
pdf=cbind(pdf,grad.th)
}
pdf
}

#inverse margingal cdf 
exp1qcdf = function(q,par){
qcdf=invfunc(q,exp1cdf,par,l.b=-10,u.b=10,tol=1e-12)
qcdf
}

#matrix of Eucledian distances
#s1 is a vector of X coordinates
#s2 is a vector of Y coordinates
DistM=function(s1,s2){
d=length(s1);
out=matrix(0,ncol=d,nrow=d);
for(i1 in 1:d){
for(i2 in i1:d){ out[i1,i2]=sqrt((s1[i1]-s1[i2])^2+(s2[i1]-s2[i2])^2); 
                 out[i2,i1]=out[i1,i2]; }}
out
}