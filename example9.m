phi=0.2;
L=10000;
Alan=200000;
mu=1;
Bw=1;
ct=1e-6;

deltaX=2500;
ki=50;
Pinitial=1000;
Pb=1800;

T=ki*Alan/(mu*Bw*deltaX)*6.33e-3;
B=diag(ones(1,4)*Alan*deltaX*ct*phi/Bw);

pn=ones(4,1)*Pinitial;
Q=zeros(4,1);
Q(1)=2*T*Pb;
deltaX=linspace(deltaX/2,7*deltaX/2,4);
z=sin(60*pi/180)*deltaX;
z=z';

T=T*(diag([3,2,2,1])+diag([-1,-1,-1],-1)+diag([-1,-1,-1],1)); % md*ft/cp den ft^3/psi/day çerirmen lazım
G=-0.433*T*z;

Pimp=zeros(10,4);
for i=1:10
    b=B*pn+Q+G;
    A=T+B;
    pimp=A\b;
    pn=pimp;
    Pimp(i,:)=pimp;
end

disp(Pimp)