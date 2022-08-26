function wave = evolve(wave,dt,H0)
        Ns = length(wave);
        
        AC=eye(Ns)+1i.*dt./2.*(H0);
        BC=eye(Ns)-1i.*dt./2.*(H0); 
        
        wave = AC\(BC*wave);