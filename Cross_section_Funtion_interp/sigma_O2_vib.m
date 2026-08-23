function cross_sect = sigma_O2_vib(E,E_old,xs_old)
% Calling: cross_sect = sigma_O2_vib(E)
% Input: E, energy [J]
% Output: cross_sect, cross sections [m^2] per energy [J] for vibrational excitation of O2 v:0->1,2,3,4 osv.  

% Interpolation 
%cross_sect = exp(interp1(E_old,log(xs_old),E,'pchip')); 
cross_sect = interp1(E_old,xs_old,E,'pchip'); 

% Set to zero below excitation energy and above sample data 
for i = 1:length(E)
    if E(i) < E_old(1)
        cross_sect(i) = 0; 
    end 
    if E(i) > E_old(end)
        cross_sect(i) = 0; 
    end 
end 

end 