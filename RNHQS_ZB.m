sigma=2.1;
N=80;
%a=zeros(N,1);

%%
lambda=633e-9;
w0=105e-6;
d=16e-6;
x=(-(N/2-1)-0.5:1:(N/2-1)+0.5)*d;
temp=1./sqrt(1+(lambda*x/(pi*w0^2)).^2);
G=temp.*exp(-(x/w0).^2.*temp.^2);
%figure;plot(G)

E=G.*exp(1i*pi*x/(2*d));
% figure;
% plot(real(E))
% hold on
% plot(imag(E),'r')
a=conj(E');
a=a/norm(a);
%%

numdt=10000;
h=0.005;
k=zeros(N,4);
t=0;
A=zeros(N,numdt);
r=0.0;ri=0.5;f=0.25;w=0.5;phi=0.0;
H=kron(eye(N/2),[1,0;0,-1])*(sigma*1+0*sigma*0.01*1i);
H(1:end-1,2:end)=H(1:end-1,2:end)-1*eye(N-1);
H(2:end,1:end-1)=H(2:end,1:end-1)-1*eye(N-1);
H(1,end)=-1;
H(end,1)=-1;
H=sparse(H./1i);
for n=1:numdt
    for i=1:4
        if i==1
            k(:,i)=RNHQS_ZBa(t,H,a,sigma,r,ri,f,w,phi);
        elseif i==2 || i==3
            k(:,i)=RNHQS_ZBa(t+h/2,H,a+h/2*k(:,i-1),sigma,r,ri,f,w,phi);
        elseif i==4
            k(:,i)=RNHQS_ZBa(t+h,H,a+h*k(:,i-1),sigma,r,ri,f,w,phi);
        end
    end
    a=a+(k(:,1)+2*k(:,2)+2*k(:,3)+k(:,4))*h/6;
    A(:,n)=a;
    t=t+h;
end
%%
figure
set(gcf,'position',[2000 400 560 840],'color','w')
resu=abs(A).^2;
subplot(4,1,1:2)
mesh(resu)
view(90,90)
subplot(4,1,3)
plot(sum(resu,1))
box on
subplot(4,1,4)
meanx=(-(N/2-1)-0.5:1:(N/2-1)+0.5)*resu./sum(resu,1);
plot(meanx)
%ylim([-1 3])
box on

% subplot(3,1,1)
% plot(abs(C(:,1)))
% subplot(3,1,2)
% plot(abs(C(:,2)))
% subplot(3,1,3)
% plot(sqrt(abs(C(:,1)).^2+abs(C(:,2)).^2))