function [E_old,xs_old] = p_data_O_3P_1D()

% Convert from cm to m
cm_to_m = 0.01; 

% Convert from [eV] to [J] 
eVtoJ = 1.602176634e-19;     

% From Itakawa 1990 fig. 5.2. Units: cm^2 
cross_sections_itakawa_1990_fig_5_2 = load('O_3P_1D.mat'); 
E_old = (10.^cross_sections_itakawa_1990_fig_5_2.O_3P_1D(:,1))*eVtoJ ;
xs_old = (10.^cross_sections_itakawa_1990_fig_5_2.O_3P_1D(:,2))*(cm_to_m)^2; 