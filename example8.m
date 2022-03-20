phi=0.2;
L=10000;
Alan=200000;
mu=1;
Bw=1;
ct=1e-6;
rw=0.25;
pwf=800;

deltaX=2500;
ki=50;
Pinitial=1000;
Pb=2000;

T=ki*Alan/(mu*Bw*deltaX)*6.33e-3;
B=diag(ones(1,4)*Alan*deltaX*ct*phi/Bw);

h=Alan/deltaX;
req=0.2*deltaX;
j=2*pi*ki*h/(mu*Bw*log(req/rw))*6.33e-3;
J=diag([0 0 j 0]);



pn=ones(4,1)*Pinitial;
Q=zeros(4,1);
Q(4)=2*T*Pb;
Q(3)=j*pwf;
Q(1)=1000;

T=T*(diag([1,2,2,3])+diag([-1,-1,-1],-1)+diag([-1,-1,-1],1)); % md*ft/cp den ft^3/psi/day çerirmen lazım

Pimp=zeros(10,4);
for i=1:10
    b=B*pn+Q;
    A=T+B+J;
    pimp=A\b;
    pn=pimp;
    Pimp(i,:)=pimp;
end

disp(Pimp)