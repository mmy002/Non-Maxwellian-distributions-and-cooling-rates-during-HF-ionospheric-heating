function [E_old,xs_old] = p_data_N2_vib_03()

% Convert from cm to m
cm_to_m = 0.01; 

% Convert from [eV] to [J] 
eVtoJ = 1.602176634e-19;   

% Figure 1 in Cambell et al. 2004 with recommended (solid line) vib. excitations. Units: 10^-17 cm^2 
N2_vib_cross_sections_v_0_3_cambell_et_al = load('N2_vib_cross_sections_v_0_3_cambell_et_al.mat'); 
E_old = (N2_vib_cross_sections_v_0_3_cambell_et_al.N2_vib_cross_sections_v_0_3_cambell_et_al(:,1))*eVtoJ; % [J]
xs_old = N2_vib_cross_sections_v_0_3_cambell_et_al.N2_vib_cross_sections_v_0_3_cambell_et_al(:,2)*1e-17*(cm_to_m)^2; % [m^2]

end 