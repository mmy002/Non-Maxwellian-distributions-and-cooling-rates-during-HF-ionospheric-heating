function [L_non_max_vib,L_maxw_vib_2nd,L_maxw_vib_L,Q_non_max,Q_maxw_2nd,Q_maxw_L] = ...
    analytical_cooling_rates(v_bins_stubbe,f0_exp_stubbe,Te_2nd,Te_L,A,nNe,nN2,nO2,T_vib,xs_data)
% Function that compute cooling rates from analytical expression for
% vibrational excitation of N2 and O2
% INPUT: 
% v_bins: [bin edges x 1] Bin energy of speed [m/s], eneries = bin egdes -1
% f0_exp_stubbe [energies x 1] Distribution function [no unit]
% Te_2nd: [1 x 1] Electron temperature from second moment of distribution [K]
% Te_L: [1 x 1] Electron temperature from Maxwellian cooling rates with the same
% total energy as non-Maxwellian [K]
% A: [1 x 1] Normalization constant 
% n_e: [1 x 1] Electron density [m^-3]
% nN2: [1 x 1] Density of molecular nitrogen [m^-3]
% nO2: [1 x 1] Density of molecular oxygen [m^-3]
% T_vib: [1 x 1] Neutral temperature [K]
% x_data: cell[1 x 27] Load and process data for cross sections
% OUTPUT:
% L_non_max_vib: [2 x 1] (First row: N2, second row: O2) Non-Maxwellian cooling rates [eV/cm^3s^1] 
% L_maxw_vib_2nd: [2 x 1] (First row: N2, second row: O2) Maxwellian cooling rates at T_2nd [eV/cm^3s^1] 
% L_maxw_vib_L: [2 x 1] (First row: N2, second row: O2) Maxwellian cooling rates at T_L [eV/cm^3s^1]  
% Q_non_max: cell[1 x 2] (First column N2, second column: O2) Non-Maxwellian Q [eV cm^3 s^-1] 
% Q_maxw_2nd: cell[1 x 2] (First column N2, second column: O2) Maxwellian Q at T_2nd [eV cm^3 s^-1]  
% Q_maxw_L: cell[1 x 2] (First column N2, second column: O2) Maxwellian Q at T_L [eV cm^3 s^-1] 

%% Some constant 
me        = 9.10938291e-31;     % Electron mass [kg]
eVtoJ     = 1.602176634e-19;    % Convert from [eV] to [J] 
k_b       = 1.380649e-23;       % Boltzmann constant [J/K] 
m_to_cm   = 100;                % Convert [m] to [cm]
K_to_eV   = 8.6175e-5;          % Convert from K to eV
E1_N2_K   = 3353;               % [K]
E1_N2     = 3353*K_to_eV*eVtoJ; % Energy of first vibrational level of N2 [J] 
deltaE_N2 = 20.6*K_to_eV*eVtoJ; % [J]
E1_O2_K   = 2239;               % [K]
E1_O2     = 2239*K_to_eV*eVtoJ; % Energy of first vibrational level of O2 [J] 
deltaE_O2 = 17.3*K_to_eV*eVtoJ; % [J]

%% Convert densities from m^-3 to cm^-3
m_neg3_to_cm_neg3 = 1e-6; 
nNe_cm_neg3 = nNe*m_neg3_to_cm_neg3;  
nN2_cm_neg3 = nN2*m_neg3_to_cm_neg3; 
nO2_cm_neg3 = nO2*m_neg3_to_cm_neg3; 

%% Energy/speed bins with same resolution as Stubbe data, but higher maximum in the enery grid

% Resolution Stubbe
delta_v = v_bins_stubbe(2) - v_bins_stubbe(1); 

% New max edges 
max_E_2 = 12.5;

% New v-bins with same resolution as Stubbe 
v      = v_bins_stubbe(end):delta_v:sqrt((2*eVtoJ*max_E_2)/me); 
v_bins = [v_bins_stubbe; v(2:end)']; 

% Energy/speed in the middle fo the bin
v_bin_mid = v_bins(1:end-1)/2 + v_bins(2:end)/2; % [m/s]
E_bin_mid = (me*v_bin_mid.^2)/2; % [J]

%% Zero pad non-Maxwellian
f0_exp_stubbe = [f0_exp_stubbe; zeros(length(v_bins)-length(v_bins_stubbe),1)]; 

%% Cross sections 

E1_O2_eV     = 2239*K_to_eV; % Energy of first vibrational level of O2 [eV] 
deltaE_O2_eV = 17.3*K_to_eV; % [eV]
max_level_O2 = 7; 
E_v_O2       = zeros(max_level_O2,1); 
for v_O2 = 1:max_level_O2
  E_v_O2(v_O2) = v_O2*E1_O2_eV - v_O2*(v_O2-1)*deltaE_O2_eV; % eV
end 

% Excitation energy for inelastic collisions [eV] from Itakawa 1986, 1989, 1990
exE_N2vib0_1  = 0.2888; 
exE_N2vib0_2  = 0.5742; 
exE_N2vib0_3  = 0.8559; 
exE_N2vib0_4  = 1.1342; 
exE_N2vib0_5  = 1.4088; 
exE_N2vib0_6  = 1.6800; 
exE_N2vib0_7  = 1.9475; 
exE_N2vib0_8  = 2.2115;
exE_N2vib0_9  = 2.4718;
exE_N2vib0_10 = 2.7284;
exE_O2a1Dg    = 0.9770; 
exE_O2b1Sgp   = 1.6270; 
exE_Ofine_1_0 = 0.0085; 
exE_Ofine_2_0 = 0.0281; 
exE_Ofine_2_1 = 0.0196; 
exE_O1S       = 4.1900; 
exE_O1D       = 1.9670; 
exE_O2vib_01  = 0.193; 
exE_O2vib_02  = E_v_O2(2); 
exE_O2vib_03  = E_v_O2(3); 
exE_O2vib_04  = E_v_O2(4); 
exE_O2vib_05  = E_v_O2(5); 
exE_O2vib_06  = E_v_O2(6); 
exE_O2vib_07  = E_v_O2(7); 

% Excitation energy. Units [J]. 
ex_E = [exE_N2vib0_1; exE_N2vib0_2; exE_N2vib0_3; exE_N2vib0_4; exE_N2vib0_5; ...
	exE_N2vib0_6; exE_N2vib0_7; exE_N2vib0_8; exE_N2vib0_9; exE_N2vib0_10; ...
	exE_O2a1Dg; exE_O2b1Sgp; ...
	exE_Ofine_1_0; exE_Ofine_2_0; exE_Ofine_2_1; exE_O1S; exE_O1D; ...
	exE_O2vib_01; exE_O2vib_02; exE_O2vib_03; exE_O2vib_04; exE_O2vib_05; ...
	exE_O2vib_06; exE_O2vib_07]*eVtoJ; 

% Cross sections [m^2]. Energy in must be in [J]. 
xs_e_N2vib0_1  = sigma_N2_vib_01( E_bin_mid, xs_data{4}(:,1),  xs_data{4}(:,2)); 
xs_e_N2vib0_2  = sigma_N2_vib_02( E_bin_mid, xs_data{5}(:,1),  xs_data{5}(:,2));
xs_e_N2vib0_3  = sigma_N2_vib_03( E_bin_mid, xs_data{6}(:,1),  xs_data{6}(:,2)); 
xs_e_N2vib0_4  = sigma_N2_vib_04( E_bin_mid, xs_data{7}(:,1),  xs_data{7}(:,2));
xs_e_N2vib0_5  = sigma_N2_vib_05( E_bin_mid, xs_data{8}(:,1),  xs_data{8}(:,2)); 
xs_e_N2vib0_6  = sigma_N2_vib_06( E_bin_mid, xs_data{9}(:,1),  xs_data{9}(:,2));
xs_e_N2vib0_7  = sigma_N2_vib_07( E_bin_mid, xs_data{10}(:,1), xs_data{10}(:,2)); 
xs_e_N2vib0_8  = sigma_N2_vib_08( E_bin_mid, xs_data{11}(:,1), xs_data{11}(:,2)); 
xs_e_N2vib0_9  = sigma_N2_vib_09( E_bin_mid, xs_data{12}(:,1), xs_data{12}(:,2)); 
xs_e_N2vib0_10 = sigma_N2_vib_010(E_bin_mid, xs_data{13}(:,1), xs_data{13}(:,2));
xs_e_O2a1Dg = sigma_O2_a1_delta_g(E_bin_mid, xs_data{19}(:,1), xs_data{19}(:,2)); 
xs_e_O2b1Sgp = sigma_O2_b1_delta_g(E_bin_mid,xs_data{20}(:,1), xs_data{20}(:,2)); 
xs_e_Ofine_1_0 = sigma_O_fine_str_bell( E_bin_mid, xs_data{21}(:,1), xs_data{21}(:,2));  
xs_e_Ofine_2_0 = sigma_O_fine_str_bell( E_bin_mid, xs_data{22}(:,1), xs_data{22}(:,2)); 
xs_e_Ofine_2_1 = sigma_O_fine_str_bell( E_bin_mid, xs_data{23}(:,1), xs_data{23}(:,2));
xs_e_O3P1D     = sigma_O_3P_1D( E_bin_mid, xs_data{24}(:,1), xs_data{24}(:,2)); 
xs_e_O3P1S     = sigma_O_3P_1S( E_bin_mid, xs_data{25}(:,1), xs_data{25}(:,2)); 
xs_e_O2vib_01  = sigma_O2_vib_02( E_bin_mid, xs_data{26}(:,1), xs_data{26}(:,2)); 
xs_e_O2vib_02  = sigma_O2_vib_02( E_bin_mid, xs_data{27}(:,1), xs_data{27}(:,2)); 
xs_e_O2vib_03  = sigma_O2_vib_03( E_bin_mid, xs_data{28}(:,1), xs_data{28}(:,2)); 
xs_e_O2vib_04  = sigma_O2_vib_04( E_bin_mid, xs_data{29}(:,1), xs_data{29}(:,2)); 
xs_e_O2vib_05  = sigma_O2_vib_05( E_bin_mid, xs_data{30}(:,1), xs_data{30}(:,2)); 
xs_e_O2vib_06  = sigma_O2_vib_06( E_bin_mid, xs_data{31}(:,1), xs_data{31}(:,2)); 
xs_e_O2vib_07  = sigma_O2_vib_07( E_bin_mid, xs_data{32}(:,1), xs_data{32}(:,2)); 
xs_i = [xs_e_N2vib0_1, xs_e_N2vib0_2, xs_e_N2vib0_3, xs_e_N2vib0_4, xs_e_N2vib0_5, ...
	xs_e_N2vib0_6, xs_e_N2vib0_7, xs_e_N2vib0_8, xs_e_N2vib0_9, xs_e_N2vib0_10, ...
	xs_e_O2a1Dg, xs_e_O2b1Sgp, ...
	xs_e_Ofine_1_0, xs_e_Ofine_2_0, xs_e_Ofine_2_1, xs_e_O3P1D, xs_e_O3P1S, ...
	xs_e_O2vib_01, xs_e_O2vib_02, xs_e_O2vib_03, xs_e_O2vib_04, xs_e_O2vib_05, ...
	xs_e_O2vib_06, xs_e_O2vib_07];

% Set cross sections to zero below threshold for excitation 
for vv = 1:length(ex_E)
  [~,ind_trh] = min(abs(E_bin_mid-ex_E(vv))); % J
  xs_i(1:ind_trh,vv) = 0; 
end  

%% Q_0v for vibrational excitation of N2 from Pavlov 1998a equation 11 

xs_N2_vib = [xs_e_N2vib0_1, xs_e_N2vib0_2, xs_e_N2vib0_3, xs_e_N2vib0_4, xs_e_N2vib0_5, ...
	     xs_e_N2vib0_6, xs_e_N2vib0_7, xs_e_N2vib0_8, xs_e_N2vib0_9, xs_e_N2vib0_10];
max_level_N2 = 10;
E_v_N2       = zeros(max_level_N2,1); 
Q0_v_N2_non_max  = zeros(max_level_N2,1);
Q0_v_N2_maxw_2nd = zeros(max_level_N2,1);
Q0_v_N2_maxw_L   = zeros(max_level_N2,1);

for v_N2 = 1:max_level_N2
  
  % Energy of the v-th vibrational level  
  E_v_N2(v_N2) = v_N2*E1_N2 - v_N2*(v_N2-1)*deltaE_N2; 
  
  % Non-Maxwellian [eV cm^3 s^-1] 
  Q0_v_N2_non_max(v_N2) = E_v_N2(v_N2)*8*pi/me^2*A * ...
			  trapz(E_bin_mid,...
				xs_N2_vib(:,v_N2).*E_bin_mid.*f0_exp_stubbe) * ...
			  (1/eVtoJ)*m_to_cm^3*(1 - exp(v_N2*E1_N2_K*((1/Te_2nd) - (1/T_vib))));

  % Maxwellian: Computed with temperature Te_2nd [eV cm^3 s^-1] 
  Q0_v_N2_maxw_2nd(v_N2) = E_v_N2(v_N2)*(8/(pi*me)).^0.5*(k_b*Te_2nd)^(-3/2) * ...
			   trapz(E_bin_mid,...
				 xs_N2_vib(:,v_N2).*E_bin_mid.*exp(-E_bin_mid./(k_b*Te_2nd))) * ...
			   (1/eVtoJ)*m_to_cm^3* (1 - exp(v_N2*E1_N2_K*((1/Te_2nd) - (1/T_vib))));  

  % Maxwellian: Computed with temperature Te_L [eV cm^3 s^-1] 
  Q0_v_N2_maxw_L(v_N2) = E_v_N2(v_N2)*(8/(pi*me)).^0.5*(k_b*Te_L)^(-3/2) * ...
			 trapz(E_bin_mid,...
			       xs_N2_vib(:,v_N2).*E_bin_mid.*exp(-E_bin_mid./(k_b*Te_L))) * ...
			 (1/eVtoJ)*m_to_cm^3*(1 - exp(v_N2*E1_N2_K*((1/Te_L) - (1/T_vib))));  

end

% N2: Equation 11 from Pavlov (1998a) [eV cm^-3 s^-1]
L_vib_N2_non_max  = nNe_cm_neg3*nN2_cm_neg3*(1-exp(-E1_N2_K/T_vib))*sum(Q0_v_N2_non_max);
L_vib_N2_maxw_2nd = nNe_cm_neg3*nN2_cm_neg3*(1-exp(-E1_N2_K/T_vib))*sum(Q0_v_N2_maxw_2nd);
L_vib_N2_maxw_L   = nNe_cm_neg3*nN2_cm_neg3*(1-exp(-E1_N2_K/T_vib))*sum(Q0_v_N2_maxw_L);

%% Q_0v for vibrational excitation of O2 from from Pavlov 1998c equation 8 

xs_O2_vib = [xs_e_O2vib_01, xs_e_O2vib_02, xs_e_O2vib_03, xs_e_O2vib_04, ...
	     xs_e_O2vib_05, xs_e_O2vib_06, xs_e_O2vib_07];
max_level_O2 = 7; 
E_v_O2 = zeros(max_level_O2,1); 
Q0_v_O2_non_max  = zeros(max_level_O2,1);
Q0_v_O2_maxw_2nd = zeros(max_level_O2,1); 
Q0_v_O2_maxw_L   = zeros(max_level_O2,1);

for v_O2 = 1:max_level_O2
  
  % Energy of the v-th vibrational level    
  E_v_O2(v_O2) = v_O2*E1_O2 - v_O2*(v_O2-1)*deltaE_O2; 
  
  % Non-Maxwellian [eV cm^3 s^-1] 
  Q0_v_O2_non_max(v_O2) = E_v_O2(v_O2)*8*pi/me^2*A * ...
			  trapz(E_bin_mid,...
				xs_O2_vib(:,v_O2).*E_bin_mid.*f0_exp_stubbe) * ...
			  (1/eVtoJ)*m_to_cm^3*(1 - exp(v_O2*E1_O2_K*((1/Te_2nd) - (1/T_vib))));

  % Maxwellian: Computed with temperature  Te_2nd [eV cm^3 s^-1] 
  Q0_v_O2_maxw_2nd(v_O2) = E_v_O2(v_O2)*(8/(pi*me)).^0.5*(k_b*Te_2nd)^(-3/2) * ...
			   trapz(E_bin_mid,...
				 xs_O2_vib(:,v_O2).*E_bin_mid.*exp(-E_bin_mid./(k_b*Te_2nd))) * ...
			   (1/eVtoJ)*m_to_cm^3 * (1 - exp(v_O2*E1_O2_K*((1/Te_2nd) - (1/T_vib))));  

  % Maxwellian: Computed with temperature  Te_L[eV cm^3 s^-1]  
  Q0_v_O2_maxw_L(v_O2) = E_v_O2(v_O2)*(8/(pi*me)).^0.5*(k_b*Te_L)^(-3/2) * ...
			 trapz(E_bin_mid,...
			       xs_O2_vib(:,v_O2).*E_bin_mid.*exp(-E_bin_mid./(k_b*Te_L))) * ...
			 (1/eVtoJ)*m_to_cm^3 * (1 - exp(v_O2*E1_O2_K*((1/Te_L) - (1/T_vib)))); 
  
end

% O2: Equation 8 from Pavlov (1998c) [eV cm^-3 s^-1] 
L_vib_O2_non_max  = nNe_cm_neg3*nO2_cm_neg3*(1-exp(-E1_O2_K/T_vib))*sum(Q0_v_O2_non_max);
L_vib_O2_maxw_2nd = nNe_cm_neg3*nO2_cm_neg3*(1-exp(-E1_O2_K/T_vib))*sum(Q0_v_O2_maxw_2nd);
L_vib_O2_maxw_L   = nNe_cm_neg3*nO2_cm_neg3*(1-exp(-E1_O2_K/T_vib))*sum(Q0_v_O2_maxw_L);

%% Q [eV*cm^3s^-1]: 
Q_non_max       = cell(1,2); 
Q_non_max{1,1}  = Q0_v_N2_non_max;  
Q_non_max{1,2}  = Q0_v_O2_non_max; 
Q_maxw_2nd      = cell(1,2); 
Q_maxw_2nd{1,1} = Q0_v_N2_maxw_2nd;  
Q_maxw_2nd{1,2} = Q0_v_O2_maxw_2nd; 
Q_maxw_L        = cell(1,2); 
Q_maxw_L{1,1}   = Q0_v_N2_maxw_L;  
Q_maxw_L{1,2}   = Q0_v_O2_maxw_L; 

%% Electron cooling rate [eV/cm^3s^1]: 
L_non_max_vib  = [L_vib_N2_non_max; L_vib_O2_non_max]; 
L_maxw_vib_2nd = [L_vib_N2_maxw_2nd; L_vib_O2_maxw_2nd]; 
L_maxw_vib_L   = [L_vib_N2_maxw_L; L_vib_O2_maxw_L]; 

end
