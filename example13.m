phi=0.2;
kha=50;
khb=40;
kv=5;
L=10000;
Alan=200000;
mu=1;
Bw=1;
ct=1e-6;
rw=0.25;

deltax=5000;
deltay=5000;
deltaz=10;

Pinitial=1000;
Pb=0;
pwf=800;

tha=kha*deltay*deltaz/(mu*Bw*deltax)*6.33e-3;
thb=khb*deltay*deltaz/(mu*Bw*deltax)*6.33e-3;
tv=kv*deltay*deltax/(mu*Bw*deltaz)*6.33e-3;

B=diag(ones(1,8)*deltay*deltaz*deltax*ct*phi/Bw);

req=0.2*deltax;
j=2*pi*kha*deltaz/(mu*Bw*log(req/rw))*6.33e-3;
J=diag([0 0 0 j 0 0 0 0]);

pn=ones(8,1)*Pinitial;
Q=zeros(8,1);
Q(1)=1000;
Q(4)=j*pwf;

z=[5 5 5 5 15 15 15 15]';


T=[2*tha+tv -tha -tha 0 -tv 0 0 0
   -tha 2*tha+tv 0 -tha 0 -tv 0 0
   -tha 0 2*tha+tv -tha 0 0 -tv 0
   0 -tha -tha 2*tha+tv 0 0 0 -tv
   -tv 0 0 0 2*thb+tv -thb -thb 0
   0 -tv 0 0 -thb 2*thb+tv 0 -thb
   0 0 -tv 0 -thb 0 2*thb+tv -thb
   0 0 0 -tv 0 -thb -thb 2*thb+tv]; % md*ft/cp den ft^3/psi/day çerirmen lazım
G=0.433*T*z;

deltat=1;
Pimp=zeros(10,8);
for i=1:10
    b=B/deltat*pn+Q+G;
    A=T+J+B/deltat;
    pimp=A\b;
    pn=pimp;
    Pimp(i,:)=pimp;
end

disp(Pimp)