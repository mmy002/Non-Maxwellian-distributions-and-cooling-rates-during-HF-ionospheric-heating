function [E_old,xs_old] = p_data_O2_rot_57_new()

% Convert from [eV] to [J] 
eVtoJ = 1.602176634e-19;     

% From LXCat, www.lxcat.net Units m^2 per eV. 
rot_O2_5_7 = load('O2_rot_5_7_new.mat'); 
E_old = rot_O2_5_7.AA(:,1)*eVtoJ; 
xs_old = rot_O2_5_7.AA(:,2); 

% Remove first value, which is zero
E_old(1) = []; 
xs_old(1) = []; 

end 