gama=0.001;
omega=2;
A=0.0*omega;
f=0.25;
v=0.01;
phi=0;

c(1)=1/sqrt(2);
c(2)=-1/sqrt(2);

numdt=100000;
h=0.002;
k=zeros(4,1);
l=k;
t=0;
C=zeros(numdt,2);
for n=1:numdt
    for i=1:4
        if i==1
            k(i)=RNHQS_c(1,t,c,gama,v,A,omega,f,phi);
            l(i)=RNHQS_c(2,t,c,gama,v,A,omega,f,phi);
        elseif i==2 || i==3
            k(i)=RNHQS_c(1,t+h/2,c+h/2*[k(i-1),l(i-1)],gama,v,A,omega,f,phi);
            l(i)=RNHQS_c(2,t+h/2,c+h/2*[k(i-1),l(i-1)],gama,v,A,omega,f,phi);
        elseif i==4
            k(i)=RNHQS_c(1,t+h,c+h*[k(i-1),l(i-1)],gama,v,A,omega,f,phi);
            l(i)=RNHQS_c(2,t+h,c+h*[k(i-1),l(i-1)],gama,v,A,omega,f,phi);
        end
    end
    c=c+[k(1)+2*k(2)+2*k(3)+k(4),l(1)+2*l(2)+2*l(3)+l(4)]*h/6;
    C(n,:)=c;
end
figure
set(gcf,'position',[2000 400 560 420])
subplot(3,1,1)
plot(abs(C(:,1)))
subplot(3,1,2)
plot(abs(C(:,2)))
subplot(3,1,3)
plot(sqrt(abs(C(:,1)).^2+abs(C(:,2)).^2))