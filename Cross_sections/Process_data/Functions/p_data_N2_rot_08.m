function [E_old,xs_old] = p_data_N2_rot_08()

% Convert from cm to m
cm_to_m = 0.01; 

% Convert from [eV] to [J] 
eVtoJ = 1.602176634e-19;     

% From Itakava 1986 fig. 5.1 and 5.2. Units cm^2 per eV. 
rot_N2_0_8 = load('N2_rot_0_8.mat'); 
E_old = (10.^(rot_N2_0_8.N2_rot_0_8(:,1)))*eVtoJ; 
xs_old = (10.^(rot_N2_0_8.N2_rot_0_8(:,2)))*(cm_to_m)^2; 