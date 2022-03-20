phi=0.2;
L=10000;
Alan=200000;
mu=1;
Bw=1;
ct=1e-6;

deltaX=2500;
ki=50;
Pinitial=1000;
Pb=2000;

T=ki*Alan/(mu*Bw*deltaX)*6.33e-3;
B=diag(ones(1,4)*Alan*deltaX*ct*phi/Bw);

pn=ones(4,1)*Pinitial;
Q=zeros(4,1);
Q(1)=2*T*Pb;

T=T*(diag([3,2,2,1])+diag([-1,-1,-1],-1)+diag([-1,-1,-1],1)); % md*ft/cp den ft^3/psi/day çerirmen lazım

Pexp=zeros(10,4);
Pimp=zeros(10,4);
Pcr_nic=zeros(10,4);
for i=1:10
    A=(1-1/2)*T+B;
    b=(B-1/2*T)*pn+Q;
    pyeni=A\b;
    pn=pyeni;
    Pcr_nic(i,:)=pyeni;
    
end
