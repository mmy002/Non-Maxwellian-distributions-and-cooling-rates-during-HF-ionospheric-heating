function [E_old,xs_old] = p_data_N2_rot_15_new()

% Convert from [eV] to [J] 
eVtoJ = 1.602176634e-19;     

% From LXCat, www.lxcat.net Units m^2 per eV. 
rot_N2_1_5 = load('N2_rot_1_5_new.mat'); 
E_old = rot_N2_1_5.AA(:,1)*eVtoJ; 
xs_old = rot_N2_1_5.AA(:,2); 

% Remove first value, which is zero
E_old(1) = []; 
xs_old(1) = []; 

end 