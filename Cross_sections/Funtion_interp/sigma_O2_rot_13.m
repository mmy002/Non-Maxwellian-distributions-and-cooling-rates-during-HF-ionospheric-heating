function cross_sect = sigma_O2_rot_13(E,E_old,xs_old)

% Calling: cross_sect = sigma_N2_rot_02(E)
% Input: E, energy [J], for interpolation 
% Output: cross_sect, cross sections for rotational excitation (0->6) of O2 [m^2] per energy [J] 

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

% % Convert from cm to m
% cm_to_m = 0.01; 
% 
% % Convert from [eV] to [J] 
% eVtoJ = 1.602176634e-19;     
% 
% % From Itakawa 1989 fig.5.1. Units cm^2 per eV.
% rot_O2_1_3 = load('O2_rot_1_3.mat'); 
% E_old = (10.^(rot_O2_1_3.O2_rot_1_3(:,1)))*eVtoJ; 
% xs_old = (10.^(rot_O2_1_3.O2_rot_1_3(:,2)))*(cm_to_m)^2;
% 
% % Find values in cross sections energy vector (eV) closest to E energy vector 
% % (E might be out of range of cross section energy vector)
% [~,closestIndex_1] = min(abs(E-E_old (1))); 
% [~,closestIndex_end] = min(abs(E-E_old (end))); 
% E_interp = E(closestIndex_1:closestIndex_end); 
% 
% %----------------------------------------------------------------------------------------------------------------
% % INTERPOLATION 
% %----------------------------------------------------------------------------------------------------------------
% sigma_interp = interp1(E_old,xs_old,E_interp,'pchip'); 
% cross_sect = [zeros((closestIndex_1-1),1); sigma_interp; zeros(length(E)-closestIndex_end,1)]; 

end 