kappa=1;
sigma=2.1;
w0=105e-6;
d=16e-6;
w0=w0/d;
kstep=500;
K=linspace(-24/w0,24/w0,kstep);
Gk2=exp(-K.^2*w0^2/16);

%%
numdt=20000;
h=2*pi/1/1000;
k=zeros(4,kstep);
l=zeros(4,kstep);
t=0;
A=zeros(kstep,numdt);
B=zeros(kstep,numdt);
r=1.0;ri=0.5;f=0.0;w=20;phi=0.0;
a=ones(1,kstep);
b=a;
for n=1:numdt
    for i=1:4
        if i==1
            k(i,:)=RNHQS_Diracdf(1,t,a,b,sigma,r,ri,f,w,phi,K);
            l(i,:)=RNHQS_Diracdf(2,t,a,b,sigma,r,ri,f,w,phi,K);
        elseif i==2 || i==3
            k(i,:)=RNHQS_Diracdf(1,t+h/2,a+h/2*k(i-1,:),b+h/2*l(i-1,:),sigma,r,ri,f,w,phi,K);
            l(i,:)=RNHQS_Diracdf(2,t+h/2,a+h/2*k(i-1,:),b+h/2*l(i-1,:),sigma,r,ri,f,w,phi,K);
        elseif i==4
            k(i,:)=RNHQS_Diracdf(1,t+h,a+h*k(i-1,:),b+h*l(i-1,:),sigma,r,ri,f,w,phi,K);
            l(i,:)=RNHQS_Diracdf(2,t+h,a+h*k(i-1,:),b+h*l(i-1,:),sigma,r,ri,f,w,phi,K);
        end
    end
    a=a+(k(1,:)+2*k(2,:)+2*k(3,:)+k(4,:))*h/6;
    b=b+(l(1,:)+2*l(2,:)+2*l(3,:)+l(4,:))*h/6;
    A(:,n)=a;
    B(:,n)=b;
    t=t+h;
end
A2=abs(A).^2;
B2=abs(B).^2;
modu=Gk2.^2*(A2+B2);
%%
figure
set(gcf,'position',[2000 400 560 840],'color','w')
subplot(4,1,1)
plot((1:numdt)*h,modu)
subplot(4,1,2)

subplot(4,1,3)

subplot(4,1,4)