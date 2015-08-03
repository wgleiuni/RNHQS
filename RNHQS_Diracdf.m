function out = RNHQS_Diracdf(n,t,a,b,sigma,r,ri,f,w,phi,k)
Sigma=sigma*(r+1i*ri*sin(w*t+phi));
if n==1
    out=Sigma*a+k.*b;
elseif n==2
    out=k.*a-Sigma*b;
end
out=out/1i;
end