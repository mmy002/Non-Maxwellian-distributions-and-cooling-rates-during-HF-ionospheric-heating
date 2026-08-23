function cross_sect = sigma_N2_vib_05(E,E_old,xs_old)

% Convert from [eV] to [J] 
eVtoJ = 1.602176634e-19; 

% Interpolation 
cross_sect = exp(interp1(E_old,log(xs_old),E,'pchip')); 

% Set to zero below excitation energy and above sample data 
for i = 1:length(E)
    if E(i)*(1/eVtoJ) < E_old*(1/eVtoJ)
        cross_sect(i) = 0; 
    end 
    if E(i)*(1/eVtoJ) > E_old*(1/eVtoJ)
        cross_sect(i) = 0; 
    end 
end 

end 