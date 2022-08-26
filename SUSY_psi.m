function psit = SUSY_psi(Lt,x,order)
    
    dx = x(2) - x(1); 
    opts.disp = 0;
        
    [~,not_box,~] = SUSY_inf_pot(Lt,x,order);
    H0 = Hammy(x,not_box);
    [Vs,Ds]=eigs(H0,1,'sm',opts);
    
    
    psit=normalize(Vs(:,1),1,dx); %instantaneous eigenstates at time t
