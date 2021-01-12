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
###lth= length(th); #####some changes
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
  ###tha=th; if(lth > 1) tha=th[!iconv]; #####some changes
  f0a = f0[!iconv]; f1a = f1[!iconv]; 
  v2a = (v0a*f1a-v1a*f0a)/(f1a-f0a);  
  f2a = func0(v2a,qa,th); 
  ###f2a = func0(v2a,qa,tha); #####some changes 
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


glint = function(func,nq,l.b,u.b){
out = 0;
dl = u.b-l.b;
gl = gausslegendre(nq);
zl = gl$nodes;
wl = gl$weights;
for(i in 1:nq) out = out + wl[i]*func(dl*zl[i]+l.b);
dl*out
}

glint.vec = function(func,nq,l.b,u.b){
out = 0;
dl = u.b-l.b;
gl = gausslegendre(nq);
zl = gl$nodes;
wl = gl$weights;
fout = func(dl*zl+l.b); 
if(is.vector(fout)) fout = matrix(fout,ncol=1)
out = apply(wl*fout,2,sum);
dl*out
}

glint.vecp = function(func,fpar,nq,l.b,u.b){
out = 0;
dl = u.b-l.b;
gl = gausslegendre(nq);
zl = gl$nodes;
wl = gl$weights;
fout = func(dl*zl+l.b,fpar); 
if(is.vector(fout)) fout = matrix(fout,ncol=1)
out = apply(wl*fout,2,sum);
dl*out
}