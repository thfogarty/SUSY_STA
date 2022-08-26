%% Fidelity and QSL for isospectral case with no cost considered


clear

Maxtime=[0.05:0.05:5.1];%[0.25:0.25:10];%(2.35);
dt=5e-3; %time step

Li = pi/sqrt(2);
Lf = 2*Li;
xmax = max(Lf/2+1);
Nx = 2001;
[x,dx,~,~] = fftdef(xmax,Nx);

for cc=1:length(Maxtime)

    maxtime=Maxtime(cc);    
    time=0:dt:maxtime;
    t_length(cc) = length(time); 
    maxtimeX=maxtime.*1;

    dtX=dt.*0.01;
    %timeX=0:dtX:maxtimeX;
    timeX=time+dtX;

    Hbar = 1;
    tf=maxtime;
    opts.disp = 0; 



    Lt=Li+(Lf-Li).*((time.^3./(maxtime.^3)).*(1+3.*(1-time./maxtime)+6.*(1-time./maxtime).^2)); %time dependent ramp K(t) - 6th order polynomial


    LtX=Li+(Lf-Li).*(timeX.^3./(maxtimeX.^3)).*(1+3.*(1-timeX./maxtimeX)+6.*(1-timeX./maxtimeX).^2); %time dependent ramp K(t+dtX)

    [box,not_box] = SUSY_inf_pot(Lt(1),x,4);
    H4 = Hammy(x,not_box);
    [Vs,~]=eigs(H4,1,'sm',opts);
    %inst0 = Ds;
    wfmi4 = normalize(Vs(:,1),1,dx);
    wfmK4=wfmi4;
    
    [box,not_box] = SUSY_inf_pot(Lt(1),x,3);
    H3 = Hammy(x,not_box);
    [Vs,~]=eigs(H3,2,'sm',opts);
    %inst0 = Ds;
    wfmi3 = normalize(Vs(:,2),1,dx);
    wfmK3=wfmi3;
    
    [box,not_box] = SUSY_inf_pot(Lt(1),x,2);
    H2 = Hammy(x,not_box);
    [Vs,~]=eigs(H2,3,'sm',opts);
    %inst0 = Ds;
    wfmi2 = normalize(Vs(:,3),1,dx);
    wfmK2=wfmi2;
  tic  
for tt = 1:length(time)

    dL = (Lt(tt) - LtX(tt))/dtX; 
    L = Lt(tt); 


    [~,not_box1] = SUSY_inf_pot(L,x,4);
    H4 = Hammy(x,not_box1);
    [~,Ds]=eigs(H4,1,'sm',opts);
    E4 = Ds(1,1);
    

    [~,not_box2] = SUSY_inf_pot(L,x,3);
    H3 = Hammy(x,not_box2);
    [~,Ds]=eigs(H3,2,'sm',opts);
    E3 = Ds(2,2);
%     
    [~,not_box3] = SUSY_inf_pot(L,x,2);
    H2 = Hammy(x,not_box3);
    [~,Ds]=eigs(H2,3,'sm',opts);
    E2 = Ds(3,3);
    
    wfmK2 = evolve(wfmK2,dt,H2); 
    wfmK3 = evolve(wfmK3,dt,H3);
    wfmK4 = evolve(wfmK4,dt,H4);
    
    Etf2(cc,tt)=real(dx.*wfmK2'*H2*wfmK2);
    Etf3(cc,tt)=real(dx.*wfmK3'*H3*wfmK3);
    Etf4(cc,tt)=real(dx.*wfmK4'*H4*wfmK4);
%     
%     

    
end
toc

    [~,not_box1] = SUSY_inf_pot(L,x,4);
    H4 = Hammy(x,not_box1);
    [Vs,~]=eigs(H4,1,'sm',opts);
     wfmf4 = normalize(Vs,1,dx);
     
    [~,not_box2] = SUSY_inf_pot(L,x,3);
    H3 = Hammy(x,not_box2);
    [Vs,~]=eigs(H3,2,'sm',opts);
    wfmf3 = normalize(Vs(:,2),1,dx);
%     
    [~,not_box3] = SUSY_inf_pot(L,x,2);
    H2 = Hammy(x,not_box3);
    [Vs,~]=eigs(H2,3,'sm',opts);
    wfmf2 = normalize(Vs(:,3),1,dx);
    
% 
% wfmf2 = sqrt(2/(15*L))*(4*cos(4*pi.*x./L)+sin(4.*pi.*x/L).*tan(pi.*x./L));
% wfmf2(abs(x)>=L/2) = 0; 
% 
% wfmf3 = sqrt(2/(45*L))*(-12*cos(2.*x.*pi./L) -3*cos(4.*x.*pi./L) -9 ).*tan(pi.*x./L);
% wfmf3(abs(x)>=L/2) = 0; 
% 
% wfmf4 = sqrt(4/(70*L))*(-cos(4.*x.*pi./L)-4*cos(2.*x.*pi./L)-3);
% wfmf4(abs(x)>=L/2) = 0; 

fid2(cc)=abs(dx.*wfmK2'*wfmf2).^2;
fid3(cc)=abs(dx.*wfmK3'*wfmf3).^2;
fid4(cc)=abs(dx.*wfmK4'*wfmf4).^2;

m_E2(cc) = mean(Etf2(cc,1:t_length(cc)));
m_E3(cc) = mean(Etf3(cc,1:t_length(cc)));
m_E4(cc) = mean(Etf4(cc,1:t_length(cc)));

   [box,not_box] = SUSY_inf_pot(Lf,x,4);
    H4 = Hammy(x,not_box);
    [Vs,~]=eigs(H4,1,'sm',opts);
    %inst0 = Ds;
    wfmf4 = normalize(Vs(:,1),1,dx);

    
    [box,not_box] = SUSY_inf_pot(Lf,x,3);
    H1 = Hammy(x,not_box);
    [Vs,~]=eigs(H1,2,'sm',opts);
    %inst0 = Ds;
    wfmf3 = normalize(Vs(:,2),1,dx);

    
    [box,not_box] = SUSY_inf_pot(Lf,x,2);
    H1 = Hammy(x,not_box);
    [Vs,~]=eigs(H1,3,'sm',opts);
    %inst0 = Ds;
    wfmf2 = normalize(Vs(:,3),1,dx);

    
    Bures_angle2 = acos(abs(dx.*wfmi2'*wfmf2));
    Bures_angle3 = acos(abs(dx.*wfmi3'*wfmf3));
    Bures_angle4 = acos(abs(dx.*wfmi4'*wfmf4));

    t_QSL2(cc) = sin(Bures_angle2).^2./(2.*m_E2(cc));
    t_QSL3(cc) = sin(Bures_angle3).^2./(2.*m_E3(cc));
    t_QSL4(cc) = sin(Bures_angle4).^2./(2.*m_E4(cc));



end

 


