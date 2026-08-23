function [E_old,xs_old] = p_data_N2_rot_04()

% Convert from cm to m
cm_to_m = 0.01; 

% Convert from [eV] to [J] 
eVtoJ = 1.602176634e-19;     

% From Itakava 1986 fig. 5.1. Units cm^2 per eV. 
rot_N2_0_4_other_part = load('N2_rot_0_4.mat'); 
x_rot_N2_0_4_other_part = (10.^(rot_N2_0_4_other_part.N2_rot_0_4(:,1)))*eVtoJ; 
y_rot_N2_0_4_other_part = (10.^(rot_N2_0_4_other_part.N2_rot_0_4(:,2)))*(cm_to_m)^2; 

% From Itakava 1986 fig. 5.2. Units cm^2 per eV. 
rot_N2_0_4_res_part = load('N2_rot_0_4_res.mat'); 
x_rot_N2_0_4_res_part = (rot_N2_0_4_res_part.N2_rot_0_4_res(:,1))*eVtoJ; 
y_rot_N2_0_4_res_part = (rot_N2_0_4_res_part.N2_rot_0_4_res(:,2))*1e-16*(cm_to_m)^2; 

% Remove some values (makes the interpolation work better...)
x_rot_N2_0_4_res_part(1:4) = []; 
y_rot_N2_0_4_res_part(1:4) = []; 
x_rot_N2_0_4_other_part(17:end) = []; 
y_rot_N2_0_4_other_part(17:end) = []; 

% Combine to one data set 
E_old = [x_rot_N2_0_4_other_part; x_rot_N2_0_4_res_part]; 
xs_old = [y_rot_N2_0_4_other_part; y_rot_N2_0_4_res_part]; 

end 