function [Pimp,Sexp,B,T]=example14
por=0.2;
k=100;
visc=1;
c=1e-5;
fvf=1;
L=1000;
deltax=333;
Area=10000;
Pinitial=1000;
swi=0.2;
sw=0.2; % for initial
swr=0.2;
Qw=zeros(3,1);
Qw(1)=426.5;
Qo=zeros(3,1);
Qo(1)=-426.5;

deltat=0.01;

ct=(1-swi)*c+swi*c;
B=diag(ones(1,3)*(Area*L*por*ct));
pn=ones(3,1)*Pinitial;
Q=Qw+Qo;
sn=ones(3,1)*sw;
d12=fvf*deltat/(Area*L*por);

% function [pimp,Sexp]=Pres(B,pn,Q,deltat,sn,d12,Qw)

Pimp=zeros(10,3);
Sexp=zeros(10,3);

for i=1:10
    [T,Tw]=Trans(sw,swi,swr,k,Area,visc,fvf,deltax);
    b=B/deltat*pn+Q;
    A=T+B/deltat;
    pimp=A\b;
    pn=pimp;
    Pimp(i,:)=pimp;
    
    
    sexp=sn+d12*(Tw*pn-Qw);
    sn=sexp;
    Sexp(i,:)=sexp;
end

end



function [T,Tw]=Trans(sw,swi,swr,k,Area,visc,fvf,deltax)
[krw,kro,~]=RelPerm(sw,swi,swr);
tw=k*Area/(visc*fvf*deltax)*krw*(6.33e-3);
to=k*Area/(visc*fvf*deltax)*kro*(6.33e-3);

Tw=tw*(diag([1,2,1])+diag([-1,-1],-1)+diag([-1,-1],1));
To=to*(diag([1,2,1])+diag([-1,-1],-1)+diag([-1,-1],1));

T=Tw+To;
end

function [krw,kro,sw]=RelPerm(sw,swi,swr)
sw=(sw-swi)/(1-swi-swr);
krw=0.2*sw^3;
kro=(1-sw)^3;
end

