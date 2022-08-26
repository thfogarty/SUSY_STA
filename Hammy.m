function H = Hammy(x,V)

Ngx=length(x);
dx=abs(x(2)-x(1));

diagV=spdiags([V],0:0,Ngx,Ngx);

% HOO is the second derivative operator with the minus sign
e=ones(Ngx,1);
%H00=spdiags([e -2*e e],-1:1,Ngx,Ngx);
H00=spdiags([(1./90)*e -(3./20)*e (3./2)*e -(49./18)*e (3./2)*e -(3./20)*e (1./90)*e],-3:3,Ngx,Ngx);

H00=-H00/dx^2;

H=H00./2+diagV;
