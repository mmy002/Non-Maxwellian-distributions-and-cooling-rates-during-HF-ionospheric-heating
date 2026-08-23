function [E_old,xs_old] = p_data_O2_vib()

% Convert from cm to m
cm_to_m = 0.01; 

% Convert from [eV] to [J] 
eVtoJ = 1.602176634e-19;     

% Figure in 6.1 in Itakawa et al. 1989. v:0->1 at lower energies
O2_vib_itakawa_v_0_1 = load('cross_section_vib_O2_itikawa_et_al_v_0_1.mat'); 
E_old_le = (10.^(O2_vib_itakawa_v_0_1.cross_section_vib_O2_itikawa_et_al_v_0_1(:,1)))*eVtoJ; % Units: [J]
xs_old_le = (10.^(O2_vib_itakawa_v_0_1.cross_section_vib_O2_itikawa_et_al_v_0_1(:,2)))*(cm_to_m)^2; % Units: m^2: v:0->1 

% Figure in 6.1 in Itakawa et al. 1989. Sum v:0->1,2,3 osv. at higher energies 
O2_vib_itakawa_sum_v_0_up = load('cross_section_vib_O2_itikawa_et_al_sum_v_0_up.mat'); 
E_old_he = (10.^(O2_vib_itakawa_sum_v_0_up.cross_section_vib_O2_itikawa_et_al_sum_v_0_up(:,1)))*eVtoJ; % [J]
xs_old_he = (10.^(O2_vib_itakawa_sum_v_0_up.cross_section_vib_O2_itikawa_et_al_sum_v_0_up(:,2)))*(cm_to_m)^2; % Units: m^2:

% Split data 
E_split_O2_vib_1 = E_old_le(1:15); 
E_split_O2_vib_2 = E_old_le(16:30); 
E_split_O2_vib_3 = E_old_le(31:42); 
E_split_O2_vib_4 = E_old_le(43:53); 
E_split_O2_vib_5 = E_old_le(54:62); 
E_split_O2_vib_6 = E_old_le(63:71); 
E_split_O2_vib_7 = E_old_le(72:79); 
E_split_O2_vib_8 = E_old_le(80:86); 
E_split_O2_vib_9 = E_old_le(90:end);  
xs_split_O2_vib_1 = xs_old_le(1:15); 
xs_split_O2_vib_2 = xs_old_le(16:30); 
xs_split_O2_vib_3 = xs_old_le(31:42); 
xs_split_O2_vib_4 = xs_old_le(43:53); 
xs_split_O2_vib_5 = xs_old_le(54:62); 
xs_split_O2_vib_6 = xs_old_le(63:71); 
xs_split_O2_vib_7 = xs_old_le(72:79); 
xs_split_O2_vib_8 = xs_old_le(80:86); 
xs_split_O2_vib_9 = xs_old_le(90:end); 

dE = E_old_le(2) - E_old_le(1); % Delta E 
num_E = 30; % Number of zeros between split data 
xs_new = zeros(num_E,1);

% Create zeros between split data 
E_new1 = (linspace(E_old_le(15)+dE,E_old_le(16)-dE,num_E))'; 
E_new2 = (linspace(E_old_le(30)+dE,E_old_le(31)-dE,num_E))'; 
E_new3 = (linspace(E_old_le(42)+dE,E_old_le(43)-dE,num_E))'; 
E_new4 = (linspace(E_old_le(53)+dE,E_old_le(54)-dE,num_E))'; 
E_new5 = (linspace(E_old_le(62)+dE,E_old_le(63)-dE,num_E))'; 
E_new6 = (linspace(E_old_le(71)+dE,E_old_le(72)-dE,num_E))'; 
E_new7 = (linspace(E_old_le(79)+dE,E_old_le(80)-dE,num_E))'; 
E_new8 = (linspace(E_old_le(86)+dE,E_old_le(90)-dE,num_E))'; 
E_new9 = (linspace(E_old_le(95)+dE,E_old_he(1)-dE,num_E))'; 

% Combine data with zeros between the split data points (zero padding)
E_old = [E_split_O2_vib_1; E_new1; E_split_O2_vib_2; E_new2; E_split_O2_vib_3; E_new3; ...
   E_split_O2_vib_4; E_new4; E_split_O2_vib_5; E_new5; E_split_O2_vib_6; E_new6; ...
   E_split_O2_vib_7; E_new7; E_split_O2_vib_8; E_new8; E_split_O2_vib_9; ....
   E_new9; E_old_he]; 

xs_old = [xs_split_O2_vib_1; xs_new; xs_split_O2_vib_2; xs_new; xs_split_O2_vib_3; xs_new; ...
   xs_split_O2_vib_4; xs_new; xs_split_O2_vib_5; xs_new; xs_split_O2_vib_6; xs_new; ...
   xs_split_O2_vib_7; xs_new; xs_split_O2_vib_8; xs_new; xs_split_O2_vib_9; ....
   xs_new; xs_old_he]; 


end