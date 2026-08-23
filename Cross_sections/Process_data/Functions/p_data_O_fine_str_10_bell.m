function [E_old,xs_old] = p_data_O_fine_str_10_bell()

% Convert from [eV] to [J] 
eVtoJ = 1.602176634e-19;  

% Converted from Bell et al. 1998 fig.2. Units: energy: eV, cross sections:m^2
cross_sect = load('O_fs_xs_1_0_bell_1998.mat'); 
E_old = (10.^cross_sect.O_fs_xs_1_0_bell_1998(:,1))*eVtoJ; 
xs_old = (10.^cross_sect.O_fs_xs_1_0_bell_1998(:,2)); 