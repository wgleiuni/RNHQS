EF=load('EF3.txt');
xLen=101;
Len=length(EF)/xLen;
EF=reshape(EF,[xLen,Len]);
[X,Y]=meshgrid((1:Len)*0.0005*1000*0.2168/(2*pi),linspace(-1,1,xLen));
figure;mesh(X,Y,(EF));
view(0,90)