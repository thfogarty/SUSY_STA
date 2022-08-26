function [inf_pot,SUSY_pot,superpot] = SUSY_inf_pot(L,x,order)

Nx = length(x);
dx = max(x)*2/Nx;

inf_pot = zeros(Nx,1);
inf_pot(abs(x) >= L/2) = 100000; 


% 
% figure(2) 
% plot(x,Orig_Trap) 
% ylim([-20 100])
e = ones(Nx,1); % Basis
Lap2 =spdiags([(1./90)*e -(3./20)*e (3./2)*e -(49./18)*e (3./2)*e -(3./20)*e (1./90)*e],-3:3,Nx,Nx)/(dx^2);
H_kin = -0.5.*Lap2;

H_pot = spdiags(inf_pot,0:0,Nx,Nx);
Ham = H_kin + H_pot;
[gs,E_gs] = eigs(Ham,order,'sm');
E_gs = diag(E_gs);


%gs = gs.*sqrt(1./(dx.*sum(abs(gs).^2)));

superpot = (order-1).*(pi/sqrt(2.*L.^2)).*tan(x.*pi./L); 

%superpot = -gradient(gs,dx)./(gs*sqrt(2)); 

SUSY_pot = superpot.^2 + gradient(superpot,dx)/sqrt(2) + E_gs(order-1);
SUSY_pot(abs(x) >= (L/2-dx)) = max(SUSY_pot(abs(x) <= L/2));
SUSY_pot(SUSY_pot > 15000) = 15000;

superpot((x) >= (L/2-dx)) =Inf; 
superpot((x) <= (-L/2+dx)) =-Inf; 















