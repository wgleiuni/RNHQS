h=0.0001;
% E=ones(1001,1);
% E(1)=0;
% E(end)=0;

%%
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

V0=-p*(exp(-((x+ws/2)/wx).^6)+exp(-((x-ws/2)/wx).^6));
V=diag(V0+1/k*1/dx^2);
V(1:end-1,2:end)=V(1:end-1,2:end)+diag(ones(N-1,1)*(-1/(2*k)*1/dx^2));
V(2:end,1:end-1)=V(2:end,1:end-1)+diag(ones(N-1,1)*(-1/(2*k)*1/dx^2));
[vec,ee]=eig(V);
ee=diag(ee);
E=(vec(:,1)+vec(:,2))/sqrt(2);
%%
step=100000;
tE=zeros(101,step);
for i=1:step
    t=i*h;
    k1=RNHQS_V(E,t);
    k2=RNHQS_V(E+h/2*k1,t+h/2);
    k3=RNHQS_V(E+h/2*k2,t+h/2);
    k4=RNHQS_V(E+h*k3,t+h);
    E=E+h/6*(k1+2*k2+2*k3+k4);
    tE(:,i)=E;
end