function [L_rot_N2_SI,L_rot_O2_SI,L_el_SI] = L_rot_el_new(v_bins_stubbe,f0_exp_stubbe,A,xs_data,...
    nNe,nN2,nO2,nO,Te_non_max,Tn)
% Function to compute electron cooling rates for rotational excitation and elastic
% colisions 

%% Some constant 
me                = 9.10938291e-31;    % Electron mass [kg]
eVtoJ             = 1.602176634*1e-19; % Convert from eV to Joule 
k_b               = 1.380649e-23;      % Boltzmann constant [J/K] 
m_to_cm           = 100;               % Convert [m] to [cm]
cm_neg3_to_m_neg3 = 1e6;               % from 1/cm^3 to 1/m^3

%% Convert densities from m^-3 to cm^-3
m_neg3_to_cm_neg3 = 1e-6; 
nNe_cm_neg3       = nNe * m_neg3_to_cm_neg3;  
nN2_cm_neg3       = nN2 * m_neg3_to_cm_neg3; 
nO2_cm_neg3       = nO2 * m_neg3_to_cm_neg3; 
nO_cm_neg3        = nO  * m_neg3_to_cm_neg3;

%% Energy/speed in the middle fo the bin
v_bin_mid = v_bins_stubbe(1:end-1)/2 + v_bins_stubbe(2:end)/2; % [m/s]
E_bin_mid = (me*v_bin_mid.^2)/2; % [J]

%% Cross sections 
KtoJ    = 1.3806e-23; % Convert from [K] to [J]
g_to_kg = 1/1000;     % Convert from [g] to [kg]

m_N2 = 4.65e-23*g_to_kg; % Mass of N2 [kg]
m_O2 = 5.31e-23*g_to_kg; % Mass of O2 [kg]
m_O  = 2.65e-23*g_to_kg; % Mass of O [kg]

% Excitation energy [eV]
exE_N2elastic  = (me/m_N2)*(Te_non_max - Tn)*KtoJ*(1/eVtoJ);  %1e-6; ????????
exE_O2elastic  = (me/m_O2)*(Te_non_max - Tn)*KtoJ*(1/eVtoJ); %1e-6; ????????
exE_Oelastic   = (me/m_O)*(Te_non_max - Tn)*KtoJ*(1/eVtoJ); %1e-6; ????????
exE_N2rot0_2   = 0.0015;
exE_N2rot0_4   = 0.005; 
exE_N2rot0_6   = 0.0104;  
exE_N2rot0_8   = 8*2.4668e-4; 
exE_O2rot1_3   = 0.0025;
exE_N2rot1_3   = 0.018; 
exE_N2rot2_4   = 0.0035; 
exE_N2rot3_5   = 0.0045; 
exE_O2rot3_5   = 0.032; 
exE_O2rot5_7   = 0.046; 
exE_O2rot7_9   = 0.066; 
exE_N2rot1_5   = 0.0069;
exE_N2rot1_7   = 0.0134; 
exE_N2rot2_6   = 0.0089;  
exE_N2rot3_7   = 0.0109; 
exE_N2rot4_6   = 0.0054;
exE_N2rot5_7   = 0.0064; 
exE_O2rot9_11  = 0.0075; 
exE_O2rot11_13 = 0.0089; 
exE_O2rot13_15 = 0.0103; 
exE_O2rot15_17 = 0.0118; 
exE_O2rot17_19 = 0.0132; 
exE_O2rot19_21 = 0.0146; 

ex_E = [exE_N2elastic;
	exE_O2elastic;
	exE_Oelastic;
	exE_N2rot0_2;
	exE_N2rot0_4; 
	exE_N2rot0_6;
	exE_N2rot0_8;
	exE_O2rot1_3;
	exE_N2rot1_3;
	exE_N2rot2_4;
	exE_N2rot3_5;
	exE_O2rot3_5;
	exE_O2rot5_7;
	exE_O2rot7_9;
	exE_N2rot1_5;
	exE_N2rot1_7;
	exE_N2rot2_6;
	exE_N2rot3_7;
	exE_N2rot4_6;
	exE_N2rot5_7;
	exE_O2rot9_11;
	exE_O2rot11_13;
	exE_O2rot13_15;
	exE_O2rot15_17;
	exE_O2rot17_19; ...
	exE_O2rot19_21]*eVtoJ; %Units [J].

% Cross sections [m^2]. Energy input must be in [J]. 
xs_e_N2elastic  = sigma_N2_el(E_bin_mid,xs_data{1}(:,1),xs_data{1}(:,2));
xs_e_O2elastic  = sigma_O2_el(E_bin_mid,xs_data{2}(:,1),xs_data{2}(:,2));
xs_e_Oelastic   = sigma_O_el_2(E_bin_mid,xs_data{3}(:,1),xs_data{3}(:,2));
xs_e_N2rot0_2   = sigma_interp(E_bin_mid,xs_data{15}(:,1),xs_data{15}(:,2)); 
xs_e_N2rot0_4   = sigma_interp(E_bin_mid,xs_data{16}(:,1),xs_data{16}(:,2));
xs_e_N2rot0_6   = sigma_interp(E_bin_mid,xs_data{17}(:,1),xs_data{17}(:,2)); 
xs_e_N2rot0_8   = sigma_N2_rot_08(E_bin_mid,xs_data{18}(:,1),xs_data{18}(:,2)); 
xs_e_O2rot1_3   = sigma_interp(E_bin_mid,xs_data{33}(:,1),xs_data{33}(:,2));
%%% new
xs_e_N2rot1_3   = sigma_interp(E_bin_mid,xs_data{34}(:,1),xs_data{34}(:,2));
xs_e_N2rot2_4   = sigma_interp(E_bin_mid,xs_data{35}(:,1),xs_data{35}(:,2));
xs_e_N2rot3_5   = sigma_interp(E_bin_mid,xs_data{36}(:,1),xs_data{36}(:,2));
xs_e_O2rot3_5   = sigma_interp(E_bin_mid,xs_data{37}(:,1),xs_data{37}(:,2));
xs_e_O2rot5_7   = sigma_interp(E_bin_mid,xs_data{38}(:,1),xs_data{38}(:,2));
xs_e_O2rot7_9   = sigma_interp(E_bin_mid,xs_data{39}(:,1),xs_data{39}(:,2));
xs_e_N2rot1_5   = sigma_interp(E_bin_mid,xs_data{40}(:,1),xs_data{40}(:,2));
xs_e_N2rot1_7   = sigma_interp(E_bin_mid,xs_data{41}(:,1),xs_data{41}(:,2));
xs_e_N2rot2_6   = sigma_interp(E_bin_mid,xs_data{42}(:,1),xs_data{42}(:,2));
xs_e_N2rot3_7   = sigma_interp(E_bin_mid,xs_data{43}(:,1),xs_data{43}(:,2));
xs_e_N2rot4_6   = sigma_interp(E_bin_mid,xs_data{44}(:,1),xs_data{44}(:,2));
xs_e_N2rot5_7   = sigma_interp(E_bin_mid,xs_data{45}(:,1),xs_data{45}(:,2));
xs_e_O2rot9_11  = sigma_interp(E_bin_mid,xs_data{46}(:,1),xs_data{46}(:,2));
xs_e_O2rot11_13 = sigma_interp(E_bin_mid,xs_data{47}(:,1),xs_data{47}(:,2));
xs_e_O2rot13_15 = sigma_interp(E_bin_mid,xs_data{48}(:,1),xs_data{48}(:,2));
xs_e_O2rot15_17 = sigma_interp(E_bin_mid,xs_data{49}(:,1),xs_data{49}(:,2));
xs_e_O2rot17_19 = sigma_interp(E_bin_mid,xs_data{50}(:,1),xs_data{50}(:,2));
xs_e_O2rot19_21 = sigma_interp(E_bin_mid,xs_data{51}(:,1),xs_data{51}(:,2));

xs_i = [xs_e_N2elastic, xs_e_O2elastic, xs_e_Oelastic, ...
	xs_e_N2rot0_2, xs_e_N2rot0_4, xs_e_N2rot0_6,   xs_e_N2rot0_8, ...
	xs_e_O2rot1_3, ...
	xs_e_N2rot1_3, xs_e_N2rot2_4, xs_e_N2rot3_5, ...
	xs_e_O2rot3_5, xs_e_O2rot5_7, xs_e_O2rot7_9, ...
	xs_e_N2rot1_5, xs_e_N2rot1_7, ...
	xs_e_N2rot2_6, xs_e_N2rot3_7, xs_e_N2rot4_6, xs_e_N2rot5_7, ...
	xs_e_O2rot9_11, xs_e_O2rot11_13, xs_e_O2rot13_15, xs_e_O2rot15_17, xs_e_O2rot17_19, xs_e_O2rot19_21];

% Set cross sections to zero below threshold for excitation 
for m = 4:length(ex_E)
  [~,ind_trh] = min(abs(E_bin_mid-ex_E(m))); % J
  xs_i(1:ind_trh,m) = 0; 
end 

%% New stuff from Eq. 1 in Stubbe and Varum 
% g_i_N2 = [1 1 1 1 3 5 7]; % Statistical weight of the i-th state of N_2
% g_i_O2 = [3 7 11 15]; % Statistical weight of the i-th stateof O_2 
% g_j_N2 = [5 9 13 17 7 9 11]; % Statistical weight of the j-th state of N_2
% g_j_O2 = [7 11 15 19]; % Statistical weight of the j-th state of O_2
% E_i_N2 = [0 0 0 0 2.4668e-4 0.0015 0.0027]*eVtoJ; % Energy of i-th state of N_2 [J]
% E_i_O2 = [3.6e-4 2.16e-3 5.4e-3 1.008e-2]*eVtoJ;  % Energy of i-th state of O_2 [J]
% E_j_N2 = [0.0015 0.005 0.0104 8*2.4668e-4 0.0027 0.005 0.0071]*eVtoJ; % Energy of j-th state of O2 [J]
% E_j_O2 = [2.16e-3 5.4e-3 1.008e-2 1.62e-2]*eVtoJ; % Energy of j-th state of O2 [J]

% Statistical weight of the i-th state or j-th state
i_N2 = [0 0 0 0 1 2 3 1 1 2 3 4 5]; 
j_N2 = [2 4 6 8 3 4 5 5 7 6 7 6 7]; 

g_i_N2 = zeros(1,length(i_N2)); 
for iN2 = 1:length(i_N2)
  g_i_N2(iN2) = 2*i_N2(iN2) + 1; 
end
g_j_N2 = zeros(1,length(j_N2)); 
for jN2 = 1:length(j_N2)
  g_j_N2(jN2) = 2*j_N2(jN2) + 1; 
end 

i_O2 = [1 3 5 7 9 11 13 15 17 19]; 
j_O2 = [3 5 7 9 11 13 15 17 19 21]; 

g_i_O2 = zeros(1,length(i_O2)); 
for iO2 = 1:length(i_O2)
  g_i_O2(iO2) = 2*i_O2(iO2) + 1; 
end 
g_j_O2 = zeros(1,length(j_O2)); 
for jO2 = 1:length(j_O2)
  g_j_O2(jO2) = 2*j_O2(jO2) + 1; 
end 

% Energy of i-th state or j-th state
E_i_N2 = [0 0 0 0 2.4668e-4 0.0015 0.0027 2.4668e-4 2.4668e-4 0.0015 0.0027 0.005 0.0071]*eVtoJ; 
E_j_N2 = [0.0015 0.005 0.0104 8*2.4668e-4 0.0027 0.005 0.0071 0.0071 0.0136 0.0104 0.0136 0.0104 0.0136]*eVtoJ; 
E_i_O2 = [3.6e-4 2.16e-3 5.4e-3 1.008e-2 1.62e-2 2.376e-2 3.276e-2 4.32e-2 5.508e-2 6.84e-2]*eVtoJ;  
E_j_O2 = [2.16e-3 5.4e-3 1.008e-2 1.62e-2 2.376e-2 3.276e-2 4.32e-2 5.508e-2 6.84e-2 8.316e-2]*eVtoJ; 

%% N2 rot.
max_level_N2 = 13;
E_v_N2 = [exE_N2rot0_2, exE_N2rot0_4, exE_N2rot0_6, exE_N2rot0_8, ...
	  exE_N2rot1_3, exE_N2rot2_4, exE_N2rot3_5, ...
	  exE_N2rot1_5, exE_N2rot1_7, ...
	  exE_N2rot2_6, exE_N2rot3_7, ...
	  exE_N2rot4_6, exE_N2rot5_7]*eVtoJ; 
xs_N2_rot = [xs_e_N2rot0_2, xs_e_N2rot0_4, xs_e_N2rot0_6, xs_e_N2rot0_8, ...
	     xs_e_N2rot1_3, xs_e_N2rot2_4, xs_e_N2rot3_5, ...
	     xs_e_N2rot1_5, xs_e_N2rot1_7, ...
	     xs_e_N2rot2_6, xs_e_N2rot3_7, ...
	     xs_e_N2rot4_6, xs_e_N2rot5_7];
Q_N2_non_max = zeros(max_level_N2,1); 

% NEW
Z_N2 = zeros(max_level_N2,1);  
for v_N2 = 1:max_level_N2
  Z_N2(v_N2) = g_j_N2(v_N2)*exp(-E_j_N2(v_N2)/(k_b*Tn)); % Partition function 
end 
Z_N2 = sum(Z_N2); 
% end NEW

% Q [eV cm^3 s^-1] 
for v_N2 = 1:max_level_N2
  
  % Non-Maxwellian 
  dist_N2            = exp(-E_i_N2(v_N2)/(k_b*Tn)) * (exp(( (-E_v_N2(v_N2)/(k_b*Te_non_max*Tn)) * (Te_non_max-Tn) )-1)); 
  Q_0v_N2_SI_non_max = g_i_N2(v_N2) * dist_N2 * E_v_N2(v_N2) * 8*pi/me^2 * A * ...
		       trapz( E_bin_mid, ...
			      xs_N2_rot(:,v_N2) .* E_bin_mid .* f0_exp_stubbe ); % [J m^3 s^-1]
  Q_N2_non_max(v_N2) = Q_0v_N2_SI_non_max*(1/eVtoJ)*m_to_cm^3;  
  
end

% Cooling rates
L_rot_N2    = (1/Z_N2) * nNe_cm_neg3 * nN2_cm_neg3 * sum(Q_N2_non_max); % [eV cm^-3 s^-1]
L_rot_N2_SI = L_rot_N2 * eVtoJ*cm_neg3_to_m_neg3;                       % [J m^-3 s^-1] 


%% O2 rot. 
max_level_O2 = 10;
E_v_O2 = [exE_O2rot1_3, exE_O2rot3_5, exE_O2rot5_7, exE_O2rot7_9, ...
	  exE_O2rot9_11, exE_O2rot11_13, exE_O2rot13_15, exE_O2rot15_17, ...
	  exE_O2rot17_19, exE_O2rot19_21] * eVtoJ; 
xs_O2_rot = [xs_e_O2rot1_3, xs_e_O2rot3_5, xs_e_O2rot5_7, xs_e_O2rot7_9, ...
	     xs_e_O2rot9_11, xs_e_O2rot11_13, xs_e_O2rot13_15, xs_e_O2rot15_17, ...
	     xs_e_O2rot17_19, xs_e_O2rot19_21]; 
Q_O2_non_max = zeros(max_level_O2,1); 

% NEW
Z_O2 = zeros(max_level_O2,1);  
for v_O2 = 1:max_level_O2
  Z_O2(v_O2) = g_j_O2(v_O2) * exp(-E_j_O2(v_O2)/(k_b*Tn)); % Partition function 
end 
Z_O2 = sum(Z_O2); 
% end NEW

% Q [eV cm^3 s^-1] 
for v_O2 = 1:max_level_O2
    
  % Non-Maxwellian 
  dist_O2 = exp(-E_i_O2(v_O2)/(k_b*Tn)) * (exp(( (-E_v_O2(v_O2)/(k_b*Te_non_max*Tn)) * (Te_non_max-Tn) )-1)); 
  Q_0v_O2_SI_non_max = g_i_O2(v_O2) * dist_O2 * E_v_O2(v_O2) * 8*pi/me^2 * A * ...
		       trapz( E_bin_mid,...
			      xs_O2_rot(:,v_O2) .* E_bin_mid.*f0_exp_stubbe );% [J m^3 s^-1]
  Q_O2_non_max(v_O2) =  Q_0v_O2_SI_non_max * (1/eVtoJ) * m_to_cm^3; 
  
end

% Cooling rates
L_rot_O2    = (1/Z_O2) * nNe_cm_neg3 * nO2_cm_neg3 * sum(Q_O2_non_max); % [eV cm^-3 s^-1]
L_rot_O2_SI = L_rot_O2 * eVtoJ * cm_neg3_to_m_neg3;                     % [J m^-3 s^-1]


%% Elastic collisions
max_level_el = 3;
E_v_el = [exE_N2elastic  exE_O2elastic  exE_Oelastic]*eVtoJ; 
xs_el =  [xs_e_N2elastic xs_e_O2elastic xs_e_Oelastic]; 
Q_el_non_max = zeros(max_level_el,1); 
for v_el = 1:max_level_el
    
  % Non-Maxwellian 
  Q_0v_el_SI_non_max = E_v_el(v_el) * 8*pi/me^2 * A * ...
		       trapz( E_bin_mid, ...
			      xs_el(:,v_el) .* E_bin_mid.*f0_exp_stubbe ); % [J m^3 s^-1]
  Q_el_non_max(v_el) =  Q_0v_el_SI_non_max*(1/eVtoJ)*m_to_cm^3;            % [eV cm^3 s^-1]

end

% Cooling rates
L_el = nNe_cm_neg3 * nN2_cm_neg3 * Q_el_non_max(1) + ...
       nNe_cm_neg3 * nO2_cm_neg3 * Q_el_non_max(2) + ...
       nNe_cm_neg3 * nO_cm_neg3  * Q_el_non_max(3); % [eV cm^-3 s^-1]
L_el_SI = L_el*eVtoJ*cm_neg3_to_m_neg3;             % [J m^-3 s^-1] 


end 
