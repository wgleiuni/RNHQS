i=10;
eval(['load EF',num2str(i),'.txt']);
EF=eval(['EF',num2str(i)]);
xLen=1001;
Len=floor(length(EF)/xLen);
p1=mod(length(EF),xLen);
if p1==0
    EF=reshape(EF,[xLen,Len]);
else
    EF(Length+1:(Len+1)*xLen)=0;
    Len=Len+1;
end
[X,Y]=meshgrid((1:Len)*0.0005*1000*0.2168/(2*pi),linspace(-10,10,xLen));
figure;mesh(X,Y,(EF));
view(0,90)

%%
i=12;
eval(['load EF',num2str(i),'.txt']);
eval(['load XP',num2str(i),'.txt']);
EF=eval(['EF',num2str(i)]);
XP=eval(['XP',num2str(i)]);
Nw=10;
xLen=100*Nw+0;
Len=floor(length(EF)/xLen);
p1=mod(length(EF),xLen);
if p1==0
    EF=reshape(EF,[xLen,Len]);
else
    EF(Len+1:(Len+1)*xLen)=0;
    Len=Len+1;
end
[X,Y]=meshgrid((1:Len)*0.0005*10000*0.2168/(2*pi),linspace(-Nw,Nw,xLen));
figure
subplot(2,1,1)
mesh(X,Y,(EF));
view(0,90)
subplot(2,1,2)
plot(XP(:,1))
set(gcf,'position',[2000 500 560 420])
%%
i=40;
eval(['load data/data9/EF',num2str(i),'.mat']);
EF=eval(['EF',num2str(i)]);
eval(['clear EF',num2str(i)]);
xLen=1000;
Len=floor(length(EF)/xLen);
p1=mod(length(EF),xLen);
if p1==0
    EF=reshape(EF,[xLen,Len]);
else
    EF(Length+1:(Len+1)*xLen)=0;
    Len=Len+1;
end
[X,Y]=meshgrid((1:Len)*0.0005*1000*0.2168/(2*pi),linspace(-10,10,xLen));
%[X,Y]=meshgrid((1:(Len-5000))*0.0005*1000*0.2168/(2*pi),linspace(-10,10,xLen));
%figure;mesh(X,Y,EF(:,1:end-5000));
figure;mesh(X,Y,EF);
view(0,90)

%%
for i=1:112
    eval(['load data/data9/EF',num2str(i),'.mat']);
    EF=eval(['EF',num2str(i)]);
    mEF=max(max(EF));
    disp(mEF);
    eval(['clear EF',num2str(i)])
end

%%
pf=zeros(112,1);
for i=1:112
    eval(['load data/data23/XP',num2str(i),'.mat']);
    xp=eval(['XP',num2str(i)]);
    temp=xp(:,4);
    pf(i)=max(temp);
    eval(['clear XP',num2str(i)])
end
if 1==1
    PF=[PF,pf];
end

%% disorder
i=12;
Rmean=load(['Rmean',num2str(i),'.txt']);
figure
set(gcf,'position',[2000 400 560 630])
subplot(5,1,1)
plot(sqrt(Rmean(:,1).^2+Rmean(:,2).^2));
subplot(5,1,2)
plot(sqrt(Rmean(:,3).^2+Rmean(:,4).^2));
subplot(5,1,3)
plot(sqrt(Rmean(:,1).^2+Rmean(:,2).^2+Rmean(:,3).^2+Rmean(:,4).^2));
subplot(5,1,4)
plot(sqrt(Rmean(:,1).^2+Rmean(:,2).^2)./sqrt(Rmean(:,1).^2+Rmean(:,2).^2+Rmean(:,3).^2+Rmean(:,4).^2));
subplot(5,1,5)
plot(sqrt(Rmean(:,3).^2+Rmean(:,4).^2)./sqrt(Rmean(:,1).^2+Rmean(:,2).^2+Rmean(:,3).^2+Rmean(:,4).^2));
%%
a2=0.5*gama; a1=0;
M=[a2,a1,0,v;-a1,a2,-v,0;0,v,-a2,-a1;-v,0,a1,-a2];
eig(M)

%%
lambda=633e-9;
w0=20e-6;
x=linspace(-16e-6*7*2,16e-6*7*2,1000);
temp=1./sqrt(1+(lambda*x/(pi*w0)).^2);
G=temp.*exp(-(x/w0).^2.*temp.^2);
figure;plot(G)

d=1e-6;
E=G.*exp(1i*pi*x/(2*d));
figure;
plot(real(E))
hold on
plot(imag(E))

%%
syms si
tH=[si,-1,0,-1;-1,-si,-1,0;0,-1,si,-1;-1,0,-1,-si];
tE=eig(tH);

%%
z=linspace(0,20,1000);
sigr=2.1;
sigi=ri*sigr*sin(w*z)*1i;
meanz=1./(2*(sigr+sigi)).*sin(2*(sigr+0*sigi).*z);
figure;plot(abs(meanz))