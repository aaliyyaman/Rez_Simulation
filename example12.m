phi=0.2;
ki=50;
L=10000;
Alan=200000;
mu=1;
Bw=1;
ct=1e-6;

deltax=5000;
deltay=5000;
deltaz=10;

Pinitial=1000;
Pb=0;

th=ki*deltay*deltaz/(mu*Bw*deltax)*6.33e-3;
tv=ki*deltay*deltax/(mu*Bw*deltaz)*6.33e-3;

B=diag(ones(1,8)*deltay*deltaz*deltax*ct*phi/Bw);

pn=ones(8,1)*Pinitial;
Q=zeros(8,1);

z=[5 5 5 5 15 15 15 15]';


T=[2*th+tv -th -th 0 -tv 0 0 0
   -th 2*th+tv 0 -th 0 -tv 0 0
   -th 0 2*th+tv -th 0 0 -tv 0
   0 -th -th 2*th+tv 0 0 0 -tv
   -tv 0 0 0 2*th+tv -th -th 0
   0 -tv 0 0 -th 2*th+tv 0 -th
   0 0 -tv 0 -th 0 2*th+tv -th
   0 0 0 -tv 0 -th -th 2*th+tv]; % md*ft/cp den ft^3/psi/day çerirmen lazım
G=0.433*T*z;

deltat=0.0001;
Pimp=zeros(10,8);
for i=1:10
    b=B/deltat*pn+Q+G;
    A=T+B/deltat;
    pimp=A\b;
    pn=pimp;
    Pimp(i,:)=pimp;
end

disp(Pimp)