function [E_old,xs_old] = p_data_O2_el()

% Convert from cm to m
cm_to_m = 0.01; 

% Convert from [eV] to [J] 
eVtoJ = 1.602176634e-19;    

% Figure 4.2 in Itikawa 1989. Units: 10^-16 cm^2 per energy eV
O2_el = load('O2_el.mat'); 
E_old = (10.^O2_el.O2_el(:,1))*eVtoJ;
xs_old = (O2_el.O2_el(:,2))*1e-16*(cm_to_m)^2; 