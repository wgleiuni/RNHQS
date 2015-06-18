function dE=RNHQS_V(E,t)
phi_=0;
p=3;
mu_=0.07;
alpha_=0.014;
ws=3.2;
wx=0.3;
f=0.25;
w=0.2168;
k=1;

N=101;
x=linspace(-2,2,N)*ws;
dx=x(end)-x(end-1);

V0=-p*(exp(-((x+ws/1)/wx).^6)+exp(-((x-ws/1)/wx).^6));
V1=-p*mu_*(exp(-((x+ws/1)/wx).^6)-exp(-((x-ws/1)/wx).^6));
VI=-p*alpha_*(exp(-((x+ws/1)/wx).^6)-exp(-((x-ws/1)/wx).^6));
F=sin(w*t)+f*sin(2*w*t+phi_);
V=V0+V1.*F+1i*VI;
E=conj(E');
dE=(V-2/dx^2*(-0.5)).*E+1/dx^2*(-0.5)*[E(2:end),0]+1/dx^2*(-0.5)*[0,E(1:end-1)];
% V=diag(V+1/k*1/dx^2);
% V(1:end-1,2:end)=V(1:end-1,2:end)+diag(ones(N-1,1)*(-1/(2*k)*1/dx^2));
% V(2:end,1:end-1)=V(2:end,1:end-1)+diag(ones(N-1,1)*(-1/(2*k)*1/dx^2));
% dE=V*E/1i;
dE=conj(dE')/1i;