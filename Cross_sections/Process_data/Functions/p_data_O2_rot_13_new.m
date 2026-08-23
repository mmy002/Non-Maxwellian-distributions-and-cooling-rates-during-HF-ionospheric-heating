function [E_old,xs_old] = p_data_O2_rot_13_new()

% Convert from [eV] to [J] 
eVtoJ = 1.602176634e-19;     

% From LXCat, www.lxcat.net Units m^2 per eV. 
rot_O2_1_3 = load('O2_rot_1_3_new.mat'); 
E_old = rot_O2_1_3.AA(:,1)*eVtoJ; 
xs_old = rot_O2_1_3.AA(:,2); 

% Remove first value, which is zero
E_old(1) = []; 
xs_old(1) = []; 

end 