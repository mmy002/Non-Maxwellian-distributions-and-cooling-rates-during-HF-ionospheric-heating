function [E_old,xs_old] = p_data_N2_vib_06()

% Convert from cm to m
cm_to_m = 0.01; 

% Convert from [eV] to [J] 
eVtoJ = 1.602176634e-19;   

% Figure 1 in Cambell et al. 2004 with recommended (solid line) vib. excitations. Units: 10^-17 cm^2 per energy eV 
N2_vib_v_0_6_cambell_et_al = load('N2_vib_cross_sections_v_0_6_cambell_et_al.mat'); 
E_old = (N2_vib_v_0_6_cambell_et_al.cross_section_vib_N2_resonance_v0_6_itikawa_et_al(:,1))*eVtoJ;
xs_old = (N2_vib_v_0_6_cambell_et_al.cross_section_vib_N2_resonance_v0_6_itikawa_et_al(:,2))*1e-17*(cm_to_m)^2; 

E_old(1) = []; 
xs_old(1) = []; 

end 