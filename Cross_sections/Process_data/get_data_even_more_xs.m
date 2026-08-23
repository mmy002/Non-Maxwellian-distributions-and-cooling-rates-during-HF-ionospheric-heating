function out = get_data_even_more_xs()
% Load data for cross sections
% Now with vib. N2 0->9, 0->10, 1->2 and rot. O2 and more vib. O2 0->2,3,4,5,6,7
% With new cross sections from Bell for fs O 
% With new cross sections for O2 rot. and N2 rot. 

% Process data 
[E_old_N2_el, xs_old_N2_el] = p_data_N2_el();
[E_old_O2_el, xs_old_O2_el] = p_data_O2_el();
[E_old_O_el, xs_old_O_el] = p_data_O_el();
[E_old_N2_vib_01, xs_old_N2_vib_01]         = p_data_N2_vib_01();
[E_old_N2_vib_02, xs_old_N2_vib_02]         = p_data_N2_vib_02();
[E_old_N2_vib_03, xs_old_N2_vib_03]         = p_data_N2_vib_03();
[E_old_N2_vib_04, xs_old_N2_vib_04]         = p_data_N2_vib_04();
[E_old_N2_vib_05, xs_old_N2_vib_05]         = p_data_N2_vib_05();
[E_old_N2_vib_06, xs_old_N2_vib_06]         = p_data_N2_vib_06();
[E_old_N2_vib_07, xs_old_N2_vib_07]         = p_data_N2_vib_07();
[E_old_N2_vib_08, xs_old_N2_vib_08]         = p_data_N2_vib_08();
[E_old_N2_vib_09, xs_old_N2_vib_09]         = p_data_N2_vib_09();
[E_old_N2_vib_010,xs_old_N2_vib_010]        = p_data_N2_vib_010();
[E_old_N2_vib_12, xs_old_N2_vib_12]         = p_data_N2_vib_12();
[E_old_N2_rot_02, xs_old_N2_rot_02]         = p_data_N2_rot_02_new();
[E_old_N2_rot_04, xs_old_N2_rot_04]         = p_data_N2_rot_04_new();
[E_old_N2_rot_06, xs_old_N2_rot_06]         = p_data_N2_rot_06_new();
[E_old_N2_rot_08, xs_old_N2_rot_08]         = p_data_N2_rot_08();
[E_old_O2_a1_delta_g, xs_old_O2_a1_delta_g] = p_data_O2_a1_delta_g();
[E_old_O2_b1_delta_g, xs_old_O2_b1_delta_g] = p_data_O2_b1_delta_g();
[E_old_O_3P_1D, xs_old_O_3P_1D]             = p_data_O_3P_1D();
[E_old_O_3P_1S, xs_old_O_3P_1S]             = p_data_O_3P_1S();
[E_old_O_fine_str_10, xs_old_O_fine_str_10] = p_data_O_fine_str_10_bell();
[E_old_O_fine_str_20, xs_old_O_fine_str_20] = p_data_O_fine_str_20_bell();
[E_old_O_fine_str_21, xs_old_O_fine_str_21] = p_data_O_fine_str_21_bell();
[E_old_O2_vib_01, xs_old_O2_vib_01]         = p_data_O2_vib_01();
[E_old_O2_vib_02, xs_old_O2_vib_02]         = p_data_O2_vib_02();
[E_old_O2_vib_03, xs_old_O2_vib_03]         = p_data_O2_vib_03();
[E_old_O2_vib_04, xs_old_O2_vib_04]         = p_data_O2_vib_04();
[E_old_O2_vib_05, xs_old_O2_vib_05]         = p_data_O2_vib_05();
[E_old_O2_vib_06, xs_old_O2_vib_06]         = p_data_O2_vib_06();
[E_old_O2_vib_07, xs_old_O2_vib_07]         = p_data_O2_vib_07();
[E_old_O2_rot_13, xs_old_O2_rot_13]         = p_data_O2_rot_13_new();

% New for rot
[E_old_N2_rot_13, xs_old_N2_rot_13] = p_data_N2_rot_13_new();
[E_old_N2_rot_24, xs_old_N2_rot_24] = p_data_N2_rot_24_new();
[E_old_N2_rot_35, xs_old_N2_rot_35] = p_data_N2_rot_35_new();
[E_old_O2_rot_35, xs_old_O2_rot_35] = p_data_O2_rot_35_new();
[E_old_O2_rot_57, xs_old_O2_rot_57] = p_data_O2_rot_57_new();
[E_old_O2_rot_79, xs_old_O2_rot_79] = p_data_O2_rot_79_new();

[E_old_N2_rot_15, xs_old_N2_rot_15] = p_data_N2_rot_15_new();
[E_old_N2_rot_17, xs_old_N2_rot_17] = p_data_N2_rot_17_new();
[E_old_N2_rot_26, xs_old_N2_rot_26] = p_data_N2_rot_26_new();
[E_old_N2_rot_37, xs_old_N2_rot_37] = p_data_N2_rot_37_new();
[E_old_N2_rot_46, xs_old_N2_rot_46] = p_data_N2_rot_46_new();
[E_old_N2_rot_57, xs_old_N2_rot_57] = p_data_N2_rot_57_new();

[E_old_O2_rot_911,  xs_old_O2_rot_911]  = p_data_O2_rot_911_new();
[E_old_O2_rot_1113, xs_old_O2_rot_1113] = p_data_O2_rot_1113_new();
[E_old_O2_rot_1315, xs_old_O2_rot_1315] = p_data_O2_rot_1315_new();
[E_old_O2_rot_1517, xs_old_O2_rot_1517] = p_data_O2_rot_1517_new();
[E_old_O2_rot_1719, xs_old_O2_rot_1719] = p_data_O2_rot_1719_new();
[E_old_O2_rot_1921, xs_old_O2_rot_1921] = p_data_O2_rot_1921_new();


out = cell(1,51);
out{1}  = [E_old_N2_el xs_old_N2_el];
out{2}  = [E_old_O2_el xs_old_O2_el];
out{3}  = [E_old_O_el xs_old_O_el];
out{4}  = [E_old_N2_vib_01 xs_old_N2_vib_01];
out{5}  = [E_old_N2_vib_02 xs_old_N2_vib_02];
out{6}  = [E_old_N2_vib_03 xs_old_N2_vib_03];
out{7}  = [E_old_N2_vib_04 xs_old_N2_vib_04];
out{8}  = [E_old_N2_vib_05 xs_old_N2_vib_05];
out{9}  = [E_old_N2_vib_06 xs_old_N2_vib_06];
out{10} = [E_old_N2_vib_07 xs_old_N2_vib_07];
out{11} = [E_old_N2_vib_08 xs_old_N2_vib_08];
out{12} = [E_old_N2_vib_09 xs_old_N2_vib_09];
out{13} = [E_old_N2_vib_010 xs_old_N2_vib_010];
out{14} = [E_old_N2_vib_12 xs_old_N2_vib_12];
out{15} = [E_old_N2_rot_02 xs_old_N2_rot_02];
out{16} = [E_old_N2_rot_04 xs_old_N2_rot_04];
out{17} = [E_old_N2_rot_06 xs_old_N2_rot_06];
out{18} = [E_old_N2_rot_08 xs_old_N2_rot_08];
out{19} = [E_old_O2_a1_delta_g xs_old_O2_a1_delta_g];
out{20} = [E_old_O2_b1_delta_g xs_old_O2_b1_delta_g];
out{21} = [E_old_O_fine_str_10 xs_old_O_fine_str_10];
out{22} = [E_old_O_fine_str_20 xs_old_O_fine_str_20];
out{23} = [E_old_O_fine_str_21 xs_old_O_fine_str_21];
out{24} = [E_old_O_3P_1D xs_old_O_3P_1D];
out{25} = [E_old_O_3P_1S xs_old_O_3P_1S];
out{26} = [E_old_O2_vib_01 xs_old_O2_vib_01];
out{27} = [E_old_O2_vib_02 xs_old_O2_vib_02];
out{28} = [E_old_O2_vib_03 xs_old_O2_vib_03];
out{29} = [E_old_O2_vib_04 xs_old_O2_vib_04];
out{30} = [E_old_O2_vib_05 xs_old_O2_vib_05];
out{31} = [E_old_O2_vib_06 xs_old_O2_vib_06];
out{32} = [E_old_O2_vib_07 xs_old_O2_vib_07];
out{33} = [E_old_O2_rot_13 xs_old_O2_rot_13];

% NEW
out{34} = [E_old_N2_rot_13 xs_old_N2_rot_13];
out{35} = [E_old_N2_rot_24 xs_old_N2_rot_24];
out{36} = [E_old_N2_rot_35 xs_old_N2_rot_35];
out{37} = [E_old_O2_rot_35 xs_old_O2_rot_35];
out{38} = [E_old_O2_rot_57 xs_old_O2_rot_57];
out{39} = [E_old_O2_rot_79 xs_old_O2_rot_79];
out{40} = [E_old_N2_rot_15,xs_old_N2_rot_15];
out{41} = [E_old_N2_rot_17,xs_old_N2_rot_17];
out{42} = [E_old_N2_rot_26,xs_old_N2_rot_26];
out{43} = [E_old_N2_rot_37,xs_old_N2_rot_37];
out{44} = [E_old_N2_rot_46,xs_old_N2_rot_46];
out{45} = [E_old_N2_rot_57,xs_old_N2_rot_57];
out{46} = [E_old_O2_rot_911,xs_old_O2_rot_911];
out{47} = [E_old_O2_rot_1113,xs_old_O2_rot_1113];
out{48} = [E_old_O2_rot_1315,xs_old_O2_rot_1315];
out{49} = [E_old_O2_rot_1517,xs_old_O2_rot_1517];
out{50} = [E_old_O2_rot_1719,xs_old_O2_rot_1719];
out{51} = [E_old_O2_rot_1921,xs_old_O2_rot_1921];
