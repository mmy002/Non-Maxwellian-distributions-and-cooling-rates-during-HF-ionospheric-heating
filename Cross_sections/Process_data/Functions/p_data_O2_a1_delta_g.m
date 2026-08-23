function [E_old,xs_old] = p_data_O2_a1_delta_g()

% Convert from cm to m
cm_to_m = 0.01; 

% Convert from [eV] to [J] 
eVtoJ = 1.602176634e-19;     

% From Itakava 1989 fig. 7.2. Units: cm^2 
cross_sections_itakawa_1989_fig_7_2 = load('O2_a1_delta_g.mat'); 
E_old = (10.^cross_sections_itakawa_1989_fig_7_2.O2_a1_delta_g(:,1))*eVtoJ ;
xs_old = (10.^cross_sections_itakawa_1989_fig_7_2.O2_a1_delta_g(:,2))*(cm_to_m)^2; 