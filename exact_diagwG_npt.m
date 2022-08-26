function [U1,E1,x,dx,H1]=exact_diagwG_npt(N,sigma_init,h_init,d_init,Xmax,Ngx,w,Vh,npt,x,dx)






omegax=w;    %omegax in Hz

mass=1; %Rb mass in kg
hbar=1;   %hbar in Js
xunit=1;        %spatial unit is in micron
ax0=sqrt(hbar./(2.*mass.*omegax)); %in meters
%[x,dx,pz,dpz]=fftdef(Xmax,Ngx); 


%No of single particle states used to diagonalise V 
Nparticles=N;
% HO potential

%Vh(1,:)=(0.5)*omegax.^2*mass*(x.^2);%harmonic potential

%Centre Gaussian
if sigma_init>0
V_init(1,:)=1./(sigma_init.*sqrt(2.*pi)).*h_init*exp((-x.^2)./(2*(sigma_init^2)));
else

% Delta fn
yap=zeros(1,length(x));
cap=find(abs(x)<1e-10);
yap(cap)=1;

V_init(1,:)=h_init*yap./dx;
end
V(1,:)=Vh;%+V_init;        

% V=V';
%x=x';


%Below code diagonalises the potential V to find both the single particle
%eigenstates and the eigenvalues;
%Depending on how many states you take initially you can build N good ones.
%This depends on the size of your x grid and finess of your mesh.


diagV=spdiags([V.'],0:0,Ngx,Ngx);

% HOO is the second derivative operator with the minus sign
%e=ones(Ngx,1);
%H00=spdiags([e -2*e e],-1:1,Ngx,Ngx);
%H00=spdiags([(1./90)*e -(3./20)*e (3./2)*e -(49./18)*e (3./2)*e -(3./20)*e (1./90)*e],-3:3,Ngx,Ngx);


%  W=ufdwt(1,npt,2);
%            %coeffs=round(W(round((npt./2)+1/2),:),12);
%            
%           
%            noff=round((npt./2)-1/2);
% 
%            coeffsX(1:(noff+1))=W(round((npt./2)+1/2),1:(noff+1));
%            coeffsX(noff+2:npt)=coeffsX(noff:-1:1);
%            
%            
%             e=ones(Ngx,1);
%             H00=spdiags([coeffsX].*e,-noff:noff,Ngx,Ngx);

e=ones(Ngx,1);
%H00=spdiags([e -2*e e],-1:1,Ngx,Ngx);
H00=spdiags([(1./90)*e -(3./20)*e (3./2)*e -(49./18)*e (3./2)*e -(3./20)*e (1./90)*e],-3:3,Ngx,Ngx);


H00=-H00/dx^2;

% find sp eigenstates
OPTS1.disp=0;
Nstates=Nparticles;

H=H00./2+diagV;
[u1,GF]=eigs(H,Nstates,'SM',OPTS1);
[G,Ind]=sort(diag(GF));


u=zeros(length(x),Nparticles);
for ii=1:Nparticles
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
  