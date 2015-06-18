EF=load('EF3.txt');
xLen=401;
Len=length(EF)/xLen;
EF=reshape(EF,[xLen,Len]);
[X,Y]=meshgrid((1:Len)*0.0001*10000*0.2168/(2*pi),linspace(-2,2,xLen));
figure;mesh(X,Y,EF);
view(0,90)