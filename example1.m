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

n=ki/(mu*phi*ct)*1/deltaX^2*6.33e-3;
A=diag([-3,-2,-2,-1])+diag([1,1,1],-1)+diag([1,1,1],1);
pn=ones(4,1)*Pinitial;
Q=zeros(4,1);
Q(1)=2*Pb*n;

PnP1=zeros(10,4);
for i=1:10
    Pnp1=pn+n*A*pn+Q;
    pn=Pnp1;
    PnP1(i,:)=Pnp1;
end
disp(PnP1)