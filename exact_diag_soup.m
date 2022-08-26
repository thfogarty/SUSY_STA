function [U1,E1,H1]=exact_diag_soup(N,x,V)

Ngx=length(x);
dx=abs(x(2)-x(1));

diagV=spdiags([V.'],0:0,Ngx,Ngx);

% HOO is the second derivative operator with the minus sign
e=ones(Ngx,1);
%H00=spdiags([e -2*e e],-1:1,Ngx,Ngx);
H00=spdiags([(1./90)*e -(3./20)*e (3./2)*e -(49./18)*e (3./2)*e -(3./20)*e (1./90)*e],-3:3,Ngx,Ngx);

H00=-H00/dx^2;

% find sp eigenstates
OPTS1.disp=0;
Nstates=N;

H=H00./2+diagV;
[u1,GF]=eigs(H,Nstates,'SM',OPTS1);
[G,Ind]=sort(diag(GF));


u=zeros(length(x),N);
for ii=1:N
    u(:,ii)=u1(:,Ind(ii))/sqrt(sum(abs(u1(:,Ind(ii))).^2)*dx);
end

H1=H;
clear u1;
E1=G(1:N);
clear V diagV H00 H GF G Ind
% for rr=1:1:N
%     if u(length(x)/2,rr)<0
%         u(:,rr)=-u(:,rr);
%     end
% end

U1=u(:,1:N);

clear u;
  