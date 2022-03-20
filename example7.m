phi=0.2;
L=10000;
A=200000;
mu=1;
Bw=1;
ct=1e-6;

deltaX=[2000 3000 1500 3500];
ki=[10 100 50 20];
Pinitial=1000;
Pb=2000;


kinterblock=(deltaX(1:end-1)+deltaX(2:end))./...
    (deltaX(1:end-1)./ki(1:end-1)+deltaX(2:end)./ki(2:end));

Tinterblock=2*kinterblock*A./(mu*Bw*(deltaX(1:end-1)+deltaX(2:end)));

T1=ki(1)*A/(mu*Bw*deltaX(1));

%unutma t1/2 yok. t1/2 2t1 dır. T ler 3/2 den başlıyor. aynı şekilde t9/2 de
%yok enson 7/2 var. toplam 3 ara t var

T=(diag([Tinterblock(1)+2*T1, Tinterblock(1:end-1)+Tinterblock(2:end),Tinterblock(end)])+ ...
    diag(-Tinterblock,1)+diag(-Tinterblock,-1))*6.33e-3;

B=diag(A*deltaX*ct*phi/Bw);

P0=ones(4,1)*Pinitial;
Q=zeros(4,1);
Q(1)=2*Pb*T1*6.33e-3;

deltaT=1; 
Pn=P0;

Pimp=zeros(10,4);
for i=1:10
    A=T+B/deltaT;
    b=(B/deltaT)*Pn+Q;
    Pnp1=A\b;
    Pn=Pnp1;
    Pimp(i,:)=Pnp1;
end

disp(Pimp)

