function [E_old,xs_old] = p_data_O_el()

% Convert from cm to m
cm_to_m = 0.01; 

% Convert from [eV] to [J] 
eVtoJ = 1.602176634e-19;    

% Figure 4.1 in Itikawa 1990. Units: 10^-16 cm^2 per energy eV
O_el = load('O_el3.mat'); 
E_old = (O_el.O_el3(:,1))*eVtoJ;
xs_old = (O_el.O_el3(:,2))*1e-16*(cm_to_m)^2; 