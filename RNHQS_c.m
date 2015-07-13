function k = RNHQS_c(n,t,c,gama,v,A,omega,f,phi)
s=-A*(sin(omega*t)+f*sin(2*omega*t+phi));
if n==1
    k=0.5*(1i*gama+s)*c(1)+v*c(2);
elseif n==2
    k=v*c(1)-0.5*(1i*gama+s)*c(2);
end
k=k/1i;
end

