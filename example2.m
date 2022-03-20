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
A=diag([1+3*n,1+2*n,1+2*n,1+n])+diag([-n,-n,-n],-1)+diag([-n,-n,-n],1);
pn=ones(4,1)*Pinitial;
Q=zeros(4,1);
Q(1)=2*Pb*n;

PnP1=zeros(10,4);

for i=1:10
    b=pn+Q;
    yeni=A\b;
    pn=yeni;
    PnP1(i,:)=yeni;
end
disp(PnP1)