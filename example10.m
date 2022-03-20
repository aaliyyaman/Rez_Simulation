phi=0.2;
L=10000;
Alan=200000;
mu=1;
Bw=1;
ct=1e-6;
rw=0.25;
pwf=800;

deltaX=3333;
ki=50;
Pinitial=1000;
Pb=2000;
h=20;

T=ki*h/(mu*Bw)*6.33e-3;
B=diag(ones(1,9)*deltaX^2*h*ct*phi/Bw);

req=0.2*deltaX;
j=2*pi*ki*h/(mu*Bw*log(req/rw))*6.33e-3;
J=diag([zeros(1,8) j]);



pn=ones(9,1)*Pinitial;
Q=zeros(9,1);
Q(3)=2*T*Pb;
Q(5)=1000;
Q(6)=2*T*Pb;
Q(9)=2*T*Pb+j*pwf;

T=T*[2 -1 0 -1 0 0 0 0 0;
     -1 3 -1 0 -1 0 0 0 0; ...
     0 -1 4 0 0 -1 0 0 0; ...
     -1 0 0 3 -1 0 -1 0 0; ...
     0 -1 0 -1 4 -1 0 -1 0;...
     0 0 -1 0 -1 5 0 0 -1; ...
     0 0 0 -1 0 0 2 -1 0;
     0 0 0 0 -1 0 -1 3 -1;
     0 0 0 0 0 -1 0 -1 4]; % md*ft/cp den ft^3/psi/day çerirmen lazım

Pimp=zeros(10,9);
for i=1:10
    b=B*pn+Q;
    A=T+B+J;
    pimp=A\b;
    pn=pimp;
    Pimp(i,:)=pimp;
end

disp(Pimp)