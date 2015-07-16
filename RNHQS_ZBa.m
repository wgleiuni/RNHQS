function k = RNHQS_ZBa(t,H,a,sigma,r,ri,f,w,phi)
N=length(diag(H));
tH=kron(eye(N/2),[1,0;0,-1])*(r*sigma*(sin(w*t)+f*sin(2*w*t+pi*phi))+ri*sigma*sin(w*t)*1i)/1i;
H=H+sparse(tH);
k=H*a;
end