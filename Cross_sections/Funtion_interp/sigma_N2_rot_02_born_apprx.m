function cross_sect = sigma_N2_rot_02_born_apprx(E)

% Calling: cross_sect = sigma_N2_rot_02(E)
% Input: E, energy [J], for interpolation 
% Output: cross_sect, cross sections for rotational excitation (0->2) of N2 [m^2] per energy [J] 

% Convert from cm to m
cm_to_m = 0.01; 

% Convert from [eV] to [J] 
eVtoJ = 1.602176634e-19;     

% From Itakava 1986 fig. 5.1. Units cm^2 per eV. Born approximation. 
rot_N2_0_2 = load('N2_rot_0_2_born.mat'); 
E_old = (10.^(rot_N2_0_2.N2_rot_0_2(:,1)))*eVtoJ; 
xs_old = (10.^(rot_N2_0_2.N2_rot_0_2(:,2)))*(cm_to_m)^2; 

% Remove some values
E_old(22:end) = []; 
xs_old(22:end) = []; 

% Find values in cross sections energy vector (eV) closest to E energy vector 
% (E might be out of range of cross section energy vector)
[~,closestIndex_1] = min(abs(E-E_old(1))); 
[~,closestIndex_end] = min(abs(E-E_old(end))); 
E_interp = E(closestIndex_1:closestIndex_end); 

%----------------------------------------------------------------------------------------------------------------
% INTERPOLATION 
%----------------------------------------------------------------------------------------------------------------
sigma_interp = interp1(E_old,xs_old,E_interp,'pchip'); 
cross_sect = [zeros((closestIndex_1-1),1); sigma_interp; zeros(length(E)-closestIndex_end,1)]; 

end 