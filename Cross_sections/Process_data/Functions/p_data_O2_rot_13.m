function [E_old,xs_old] = p_data_O2_rot_13()

% Convert from cm to m
cm_to_m = 0.01; 

% Convert from [eV] to [J] 
eVtoJ = 1.602176634e-19;    

% From Itakawa 1989 fig.5.1. Units cm^2 per eV.
rot_O2_1_3 = load('O2_rot_1_3.mat'); 
E_old = (10.^(rot_O2_1_3.O2_rot_1_3(:,1)))*eVtoJ; 
xs_old = (10.^(rot_O2_1_3.O2_rot_1_3(:,2)))*(cm_to_m)^2;