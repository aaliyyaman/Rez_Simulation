phi=0.2;
L=10000;
Alan=200000;
mu=1;
Bw=1;
ct=1e-6;

deltaX=2500;
ki=[10 100 50 20];
Pinitial=1000;
Pb=2000;

ki=[ki(1) 2*(1./ki(1:end-1)+1./ki(2:end)).^(-1)];
T=ki*Alan/(mu*Bw*deltaX)*6.33e-3;
B=diag(ones(1,4)*Alan*deltaX*ct*phi/Bw);

pn=ones(4,1)*Pinitial;
Q=zeros(4,1);
Q(1)=2*T(1)*Pb;

T=diag([2*T(1)+T(2),T(2:end-1)+T(3:end),T(end)])+diag(-T(2:end),-1)+diag(-T(2:end),1); 

Pimp=zeros(10,4);
for i=1:10
    b=B*pn+Q;
    A=T+B;
    pimp=A\b;
    pn=pimp;
    Pimp(i,:)=pimp;
end

disp(Pimp)