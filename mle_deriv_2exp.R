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

#double exponential density (for r.v. V0=Exp(th1)-Exp(th2))
fpdf = function(z,par,gradth=F){
th1=par[1]; th2=par[2];
t1=exp(-th1*z*(z>0))
t2=exp(th2*z*(z<0))
t12=t1*t2
out = th1*th2*t12/(th1+th2)
#grad.th=NULL;
if(gradth==T){
tem1=th2^2*(t12)/(th1+th2)^2-z*(z>0)*out
tem2=th1^2*(t12)/(th1+th2)^2+z*(z<0)*out
grad.th=cbind(tem1,tem2)
out = cbind(out,grad.th)
}
out
}



llikf = function(u,Dist,par,grad=F,zsc=F){

if(min(Dist)<0){print("error in Dist: negative distances"); return(NULL)}
diag(Dist) = 0;

if(is.vector(u)) u = matrix(u,nrow=1)
n = nrow(u);
d = ncol(u);
th =par[3]; al=par[4];
par=par[1:2];
th1=par[1]; th2=par[2];

sigma=exp(-th*Dist^al);

dsig1=-Dist^al*sigma;
Dist1=Dist; diag(Dist1)=1;
dsig2=-th*Dist^al*log(Dist1)*sigma;
siginv=solve(sigma); 
invs1=siginv%*%dsig1;
invs2=siginv%*%dsig2;
dets=det(sigma);

c32=sum(siginv);
c31=sqrt(c32);
c3.3 =-sum(invs1%*%siginv); 
c3.4 =-sum(invs2%*%siginv);


llik=0; lgrad = 0;
for(i in 1:n){

if(zsc==F) z=exp1qcdf(u[i,],par); #changed!
if(zsc==T) z=u[i,]; #changed!
pdf1  = exp1pdf(z,par);
lden  = sum(log(pdf1));

sigz=siginv%*%z; 
c1=t(z)%*%sigz;
c2=sum(sigz);
fpar1=(c2+th2)/c31;
fpar2=(c2-th1)/c31;
pw1=fpar1*fpar1;
pw2=fpar2*fpar2;
ex1=exp(.5*pw1); 
ex2=exp(.5*pw2);
px1=pnorm(-fpar1);
px2=pnorm(fpar2);
eliki=ex1*px1+ex2*px2;
lliki=log(eliki); 
llik = llik+lden+.5*c1-lliki;

if(grad==T){
cdf1=exp1cdf(z,par,T);
pdf1=exp1pdf(z,par,T);
zth=-cdf1[,2:3]/pdf1[,1];

sigz1=siginv%*%zth;
ssz1 =apply(sigz1,2,sum);
c1.12=2*t(zth)%*%sigz;
c1.3 =-t(sigz)%*%dsig1%*%sigz;
c1.4 =-t(sigz)%*%dsig2%*%sigz;
c1s.12=(ssz1+c(0,1))/c31;
c2s.12=(ssz1+c(-1,0))/c31;
tem.3 =-sum(invs1%*%sigz)/c31;
tem.4 =-sum(invs2%*%sigz)/c31;
c33=c31*c32;
c1s.3=tem.3-c3.3*(c2+th2)/2/c33
c2s.3=tem.3-c3.3*(c2-th1)/2/c33
c1s.4=tem.4-c3.4*(c2+th2)/2/c33
c2s.4=tem.4-c3.4*(c2-th1)/2/c33
lnum=px1*ex1*fpar1*c(c1s.12,c1s.3)+px2*ex2*fpar2*c(c2s.12,c2s.3)+(c(c2s.12,c2s.3)-c(c1s.12,c1s.3))/sqrt(2*pi);
lnum=px1*ex1*fpar1*c(c1s.12,c1s.3,c1s.4)+px2*ex2*fpar2*c(c2s.12,c2s.3,c2s.4)+(c(c2s.12,c2s.3,c2s.4)-c(c1s.12,c1s.3,c1s.4))/sqrt(2*pi);
#lgrad=lgrad+.5*c(c1.12,c1.3)-lnum/eliki;
lgrad=lgrad+.5*c(c1.12,c1.3,c1.4)-lnum/eliki;
lgrad[1:2]=lgrad[1:2]+apply((pdf1[,4]*zth+pdf1[,2:3])/pdf1[,1],2,sum);
}

}
llik = llik+.5*n*((d-1)*log(2*pi)+log(dets)+log(c32))+n*(log(th1+th2)-log(th1)-log(th2));
if(grad==T){
lgrad[1:2] = lgrad[1:2] + n*(1/sum(par)-1/par);
lgrad[3] = lgrad[3] +.5*n*sum(diag(invs1)) +.5*n*c3.3/c32;
lgrad[4] = lgrad[4] +.5*n*sum(diag(invs2)) +.5*n*c3.4/c32;##
llik = c(llik,lgrad);
}

llik
}