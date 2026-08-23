function [E_old,xs_old] = p_data_N2_vib_01()

% Convert from cm to m
cm_to_m = 0.01; 

% Convert from [eV] to [J] 
eVtoJ = 1.602176634e-19;        

% Ångstrøm to m
angstrom_to_m = 1e-10; 

%----------------------------------------------------------------------------------------------------------------
% RESONANCE REGION 
%----------------------------------------------------------------------------------------------------------------
% Figure 1 in Cambell et al. 2004 with recommended (solid line) vib.
% excitations from v = 0 to v =1. Units: 10^-17 cm^2 per energy eV
N2_vib_cross_sections_v_0_1_cambell_et_al = load('N2_vib_cross_sections_v_0_1_cambell_et_al.mat'); 
N2_vib_x_v_0_1_res_region_old = N2_vib_cross_sections_v_0_1_cambell_et_al.N2_vib_cross_sections_v_0_1_cambell_et_al(:,1);
N2_vib_y_v_0_1_res_region_old = N2_vib_cross_sections_v_0_1_cambell_et_al.N2_vib_cross_sections_v_0_1_cambell_et_al(:,2); 
% Remove negative values 
N2_vib_x_v_0_1_res_region_old(6) = []; 
N2_vib_y_v_0_1_res_region_old(6) = []; 
% Find duplicated values and remove them or x not increasing 
N2_vib_x_v_0_1_res_region_old(14) = []; 
N2_vib_y_v_0_1_res_region_old(14) = []; 
N2_vib_x_v_0_1_res_region_old(18) = []; 
N2_vib_y_v_0_1_res_region_old(18) = []; 
N2_vib_x_v_0_1_res_region_old(34) = []; 
N2_vib_y_v_0_1_res_region_old(34) = []; 
N2_vib_x_v_0_1_res_region_old(46) = []; 
N2_vib_y_v_0_1_res_region_old(46) = []; 
anyDuplicates = ~all(diff(sort(N2_vib_x_v_0_1_res_region_old(N2_vib_x_v_0_1_res_region_old ~= 0))));

% Remove low energy tail (it is not so accurate, the one from Figure 2 in Cambell et al. 2004 is much better)
N2_vib_x_v_0_1_res_region = N2_vib_x_v_0_1_res_region_old(12:end); 
N2_vib_y_v_0_1_res_region = N2_vib_y_v_0_1_res_region_old(12:end); 

%----------------------------------------------------------------------------------------------------------------
% LOW ENERGY TAIL REGION 
%----------------------------------------------------------------------------------------------------------------
% Figure 2 in Cambell et al. 2004 with recommended (big dotted line) vib.
% excitations from v=0->1. Units: angstrøm^2 per energy eV
N2_vib_cross_sections_v_0_1_cambell_et_al_tail = load('N2_vib_cross_sections_v_0_1_cambell_et_al_tail.mat'); 
N2_vib_x_v_0_1_low_E_tail_old = N2_vib_cross_sections_v_0_1_cambell_et_al_tail.N2_vib_cross_sections_v_0_1_cambell_et_al_tail(:,1); 
N2_vib_y_v_0_1_low_E_tail_old = 10.^(N2_vib_cross_sections_v_0_1_cambell_et_al_tail.N2_vib_cross_sections_v_0_1_cambell_et_al_tail(:,2)); 

N2_vib_x_v_0_1_low_E_tail = N2_vib_x_v_0_1_low_E_tail_old(1:28); 
N2_vib_y_v_0_1_low_E_tail = N2_vib_y_v_0_1_low_E_tail_old(1:28); 

%----------------------------------------------------------------------------------------------------------------
% HIGH ENERGY TAIL REGION 
%----------------------------------------------------------------------------------------------------------------
% Figure 6.1 in Itikawa et al. 1986. Units: cm^2 per energy eV
cross_section_vib_N2_v0_1_itikawa_et_al = load('cross_section_vib_N2_v0_1_itikawa_et_al.mat'); 
N2_vib_x_v_0_1_high_E_tail = 10.^(cross_section_vib_N2_v0_1_itikawa_et_al.cross_section_vib_N2_v0_1_itikawa_et_al(9:end,1)); % [eV]
N2_vib_y_v_0_1_high_E_tail = 10.^(cross_section_vib_N2_v0_1_itikawa_et_al.cross_section_vib_N2_v0_1_itikawa_et_al(9:end,2)); % [cm^2]

%----------------------------------------------------------------------------------------------------------------
% CONVERT TO [m^2]
%----------------------------------------------------------------------------------------------------------------
N2_vib_y_v_0_1_res_region_m_2 = N2_vib_y_v_0_1_res_region*1e-17*(cm_to_m)^2; % [m^2]
N2_vib_y_v_0_1_low_E_tail_m_2 = N2_vib_y_v_0_1_low_E_tail*(angstrom_to_m)^2; % [m^2]
N2_vib_y_v_0_1_high_E_tail_m_2 = N2_vib_y_v_0_1_high_E_tail*(cm_to_m)^2;     % [m^2]

%----------------------------------------------------------------------------------------------------------------
% COMBINE TO ONE DATA SET
%----------------------------------------------------------------------------------------------------------------
E_old = ([N2_vib_x_v_0_1_low_E_tail; N2_vib_x_v_0_1_res_region; N2_vib_x_v_0_1_high_E_tail])*eVtoJ; % Convert to [J]
xs_old = [N2_vib_y_v_0_1_low_E_tail_m_2; N2_vib_y_v_0_1_res_region_m_2; N2_vib_y_v_0_1_high_E_tail_m_2]; 

end 