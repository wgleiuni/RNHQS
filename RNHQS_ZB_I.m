function I=RNHQS_ZB_I(T,sigma,sigmai,w,kappa,k)
%t=linspace(0,2*pi/w,T);
dt=2*pi/w/T;
t=(0:T-1)*dt;
s2=(sigma*(1+1i*sigmai*(sin(w*t)))).^2;
k2=kappa^2*k.^2;
tI=zeros(T,length(k));
for i=1:length(k)
    tI(:,i)=s2+k2(i);
end
tI=sqrt(tI)*dt;
I=zeros(T,length(k));
I(1,:)=tI(1,:);
for i=2:T
    I(i,:)=sum(tI(1:i,:));
end
end