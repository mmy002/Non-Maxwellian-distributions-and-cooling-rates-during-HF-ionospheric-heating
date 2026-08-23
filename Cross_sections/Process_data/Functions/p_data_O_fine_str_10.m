function [E_old,xs_old] = p_data_O_fine_str_10()

% Convert from cm to m
cm_to_m = 0.01; 

% Convert from [eV] to [J] 
eVtoJ = 1.602176634e-19;     

% From Itakava 1990 fig. 5.1. Units: cm^2 
cross_sections_itakawa_1990_fig_5_1 = load('O_fine_str_J_1_to_0.mat'); 
E_old = (10.^cross_sections_itakawa_1990_fig_5_1.O_fine_str_J_1_to_0(:,1))*eVtoJ ;
xs_old = (10.^cross_sections_itakawa_1990_fig_5_1.O_fine_str_J_1_to_0(:,2))*(cm_to_m)^2; 