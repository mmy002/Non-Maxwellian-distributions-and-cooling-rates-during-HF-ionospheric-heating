function out = get_data_withN2vib()
% Load data for cross sections
% Original one with less cross sections 
% Vib. O2 only 0->1
% Vib. N2 0->,1,2,3,4,5,6,7,8
% No rot. O2

% Process data 
[E_old_N2_el,         xs_old_N2_el]         = p_data_N2_el(); 
[E_old_O2_el,         xs_old_O2_el]         = p_data_O2_el(); 
[E_old_O_el,          xs_old_O_el]          = p_data_O_el(); 
[E_old_N2_vib_01,     xs_old_N2_vib_01]     = p_data_N2_vib_01(); 
[E_old_N2_vib_02,     xs_old_N2_vib_02]     = p_data_N2_vib_02(); 
[E_old_N2_vib_03,     xs_old_N2_vib_03]     = p_data_N2_vib_03(); 
[E_old_N2_vib_04,     xs_old_N2_vib_04]     = p_data_N2_vib_04(); 
[E_old_N2_vib_05,     xs_old_N2_vib_05]     = p_data_N2_vib_05(); 
[E_old_N2_vib_06,     xs_old_N2_vib_06]     = p_data_N2_vib_06(); 
[E_old_N2_vib_07,     xs_old_N2_vib_07]     = p_data_N2_vib_07(); 
[E_old_N2_vib_08,     xs_old_N2_vib_08]     = p_data_N2_vib_08();
[E_old_N2_rot_02,     xs_old_N2_rot_02]     = p_data_N2_rot_02(); 
[E_old_N2_rot_04,     xs_old_N2_rot_04]     = p_data_N2_rot_04(); 
[E_old_N2_rot_06,     xs_old_N2_rot_06]     = p_data_N2_rot_06(); 
[E_old_N2_rot_08,     xs_old_N2_rot_08]     = p_data_N2_rot_08(); 
[E_old_O2_a1_delta_g, xs_old_O2_a1_delta_g] = p_data_O2_a1_delta_g(); 
[E_old_O2_b1_delta_g, xs_old_O2_b1_delta_g] = p_data_O2_b1_delta_g(); 
[E_old_O_3P_1D,       xs_old_O_3P_1D]       = p_data_O_3P_1D(); 
[E_old_O_3P_1S,       xs_old_O_3P_1S]       = p_data_O_3P_1S(); 
[E_old_O_fine_str_10, xs_old_O_fine_str_10] = p_data_O_fine_str_10(); 
[E_old_O_fine_str_20, xs_old_O_fine_str_20] = p_data_O_fine_str_20(); 
[E_old_O_fine_str_21, xs_old_O_fine_str_21] = p_data_O_fine_str_21();
[E_old_O2_vib,        xs_old_O2_vib]        = p_data_O2_vib(); 

out = cell(1,23);
out{1}  = [E_old_N2_el         xs_old_N2_el]; 
out{2}  = [E_old_O2_el         xs_old_O2_el]; 
out{3}  = [E_old_O_el          xs_old_O_el]; 
out{4}  = [E_old_N2_vib_01     xs_old_N2_vib_01]; 
out{5}  = [E_old_N2_vib_02     xs_old_N2_vib_02]; 
out{6}  = [E_old_N2_vib_03     xs_old_N2_vib_03]; 
out{7}  = [E_old_N2_vib_04     xs_old_N2_vib_04]; 
out{8}  = [E_old_N2_vib_05     xs_old_N2_vib_05]; 
out{9}  = [E_old_N2_vib_06     xs_old_N2_vib_06]; 
out{10} = [E_old_N2_vib_07     xs_old_N2_vib_07]; 
out{11} = [E_old_N2_vib_08     xs_old_N2_vib_08]; 
out{12} = [E_old_N2_rot_02     xs_old_N2_rot_02];  
out{13} = [E_old_N2_rot_04     xs_old_N2_rot_04]; 
out{14} = [E_old_N2_rot_06     xs_old_N2_rot_06];
out{15} = [E_old_N2_rot_08     xs_old_N2_rot_08];
out{16} = [E_old_O2_a1_delta_g xs_old_O2_a1_delta_g];
out{17} = [E_old_O2_b1_delta_g xs_old_O2_b1_delta_g];
out{18} = [E_old_O_fine_str_10 xs_old_O_fine_str_10];
out{19} = [E_old_O_fine_str_20 xs_old_O_fine_str_20];
out{20} = [E_old_O_fine_str_21 xs_old_O_fine_str_21];
out{21} = [E_old_O_3P_1D       xs_old_O_3P_1D];
out{22} = [E_old_O_3P_1S       xs_old_O_3P_1S];
out{23} = [E_old_O2_vib        xs_old_O2_vib]; 
