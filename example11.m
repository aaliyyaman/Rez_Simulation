phi=0.2;
L=10000;
mu=1;
Bw=1;
ct=1e-6;
rw=0.25;
pwf=800;

deltaX=3333;
deltaY=3333;
h=20;
kx=50;ky=50;kz=5;
Pinitial=1000;
Pb=2000;


T=kx*h/(mu*Bw)*6.33e-3;
B=diag(ones(1,9)*deltaX^2*h*ct*phi/Bw);

req=0.28*(sqrt(sqrt(kz/ky)*deltaY^2+sqrt(ky/kz)*h^2)/(nthroot(kz/ky,4)+nthroot(ky/kz,4)));
j=2*pi*sqrt(ky*kz)*deltaX/(mu*Bw*(log(req/rw)-0.75))*6.33e-3;
J=diag([zeros(1,3) j j j zeros(1,3)]);



pn=ones(9,1)*Pinitial;
Q=zeros(9,1);
Q(3)=2*T*Pb;
Q(4)=j*pwf;
Q(5)=j*pwf;
Q(6)=2*T*Pb+j*pwf;
Q(9)=2*T*Pb;

T=T*[2 -1 0 -1 0 0 0 0 0;
     -1 3 -1 0 -1 0 0 0 0;
     0 -1 4 0 0 -1 0 0 0;
     -1 0 0 3 -1 0 -1 0 0;
     0 -1 0 -1 4 -1 0 -1 0;
     0 0 -1 0 -1 5 0 0 -1;
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