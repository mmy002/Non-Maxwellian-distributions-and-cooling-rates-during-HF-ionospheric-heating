function [f0_save,v_bins,Te_non_max,Te_grad,A_save,dI_sum,dIe,deg_Ie] = ...
	 compute_f0_stubbe_more_xs(E0_field,n_e,Tn,nO,nN2,nO2,f_HF,xs_data,max_E)
% This function computes non-Maxwellian distribution during ionospheric HF heating,
% re-implemented from Stubbe 1981. 
% INPUT: 
% E0_field: [1 x 1] Electric field amplitude [V/m]
% n_e: [1 x 1] Electron density [m^-3]
% Tn: [1 x 1] Neutral temperature [K]
% nO: [1 x 1] Density of atomic oxygen [m^-3]
% nN2: [1 x 1] Density of molecular nitrogen [m^-3]
% nO2: [1 x 1] Density of molecular oxygen [m^-3]
% f_HF: [1 x 1] Radio wave frequency [Hz] 
% x_data: cell[1 x 27] Load and process data for cross sections
% max_E: [1 x 1] Maximum energy in the energy grid [eV]
% OUTPUT:
% f0_save: [energies x 1] Distribution function [no unit]
% v_bins: [bin edges x 1] Bin energy of speed [m/s], eneries = bin egdes -1
% Te_non_max: [1 x 1] Temperature from second moment of distribution [K]
% Te_grad: [energies x 1] Gradient temperature 
% A_save: [1 x 1] Normalization constant 
% dI_sum: [energies x 1] Sum over all different cross sections for gain 
% and loss for each bin [#e/(m^3*s)]
% dIe: [energies x levels] Number of electrons degrading from higher to lower
% energies [#e/(m^3*s)]
% deg_Ie: [energies x levels] Degraded electrons at lower energies
% distributed into their respective bins [#e/(m^3*s)]

% This functions includes cross sections for: vibrational excitation of N2
% for 0 to 1,2,3,4,5,6,7,8,9,10, electronic excitation of O2 and O, fine
% structure excitation of O, and vibrational excitation of O2 for 0 to 1,2,3,4,5,6,7

%% Some constant 
me      = 9.10938291e-31;   % Electron mass [kg]
eVtoJ   = 1.602176634e-19;  % Convert from [eV] to [J] 
k_b     = 1.380649e-23;     % Boltzmann constant [J/K] 
q_e     = 1.602176634e-19;  % Elementary charge of 1 electron [C]
cm_to_m = 0.01;             % Convert from [cm] to [m]
g_to_kg = 1/1000;           % Convert from [g] to [kg]
K_to_eV = 8.6175e-5;        % Convert from K to eV

%% Geophysical and ionospheric
Te_start = 1000;            % Start temperature for the iteration scheme

%% Radio wave heating parameters
w_HF = 2*pi*f_HF; % Angular frequency [1/s]
         
%% Energy/speed bins 

% Resolution - more points = smaller resolution 
res = 8000; 

% Edges for speed bins [m/s]
v_bins = (linspace(0,sqrt((2*eVtoJ*max_E)/me),res))';

% Convert speed bins to energy bins [J]
E_bins = (me*v_bins.^2)/2;

% Resolution of speed (constant) [m/s]
delta_v = v_bins(2)-v_bins(1); 

% Resolution of energy (varying) [J]
delta_E = zeros(length(E_bins)-1,1); 
delta_E_eV = zeros(length(E_bins)-1,1); 
for i = 1:length(E_bins)-1
  delta_E(i,1) = E_bins(i+1)-E_bins(i); 
  delta_E_eV(i,1) = (E_bins(i+1)-E_bins(i))/eVtoJ; 
end 

% Velocities at the bin-edges
v_bin_higher = v_bins(2:end); 
v_bin_lower  = v_bins(1:end-1); 

% Energy/speed in the middle fo the bin
v_bin_mid = v_bins(1:end-1)/2 + v_bins(2:end)/2; % [m/s]
E_bin_mid = (me*v_bin_mid.^2)/2; % [J]

%% Maxwellian distribution: [#e/m^3 * s/m]
f_max_v = n_e*((me/(2*pi*k_b*Te_start))^(3/2)*4*pi*v_bin_mid.^2.*exp(-me*v_bin_mid.^2/(2*k_b*Te_start))); 

%% Cross sections 

E1_O2 = 2239*K_to_eV; % Energy of first vibrational level of O2 [eV] 
deltaE_O2 = 17.3*K_to_eV; % [eV]
max_level_O2 = 7; 
E_v_O2 = zeros(max_level_O2,1); 
for v_O2 = 1:max_level_O2
  E_v_O2(v_O2) = v_O2*E1_O2 - v_O2*(v_O2-1)*deltaE_O2; % [eV]
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

% Excitation energy [J] 
ex_E = [exE_N2vib0_1; exE_N2vib0_2; exE_N2vib0_3; exE_N2vib0_4; exE_N2vib0_5; ...
	exE_N2vib0_6; exE_N2vib0_7; exE_N2vib0_8; exE_N2vib0_9; exE_N2vib0_10; ...
	exE_O2a1Dg; exE_O2b1Sgp; ...
	exE_Ofine_1_0; exE_Ofine_2_0; exE_Ofine_2_1; exE_O1S; exE_O1D;...
	exE_O2vib_01; exE_O2vib_02; exE_O2vib_03; exE_O2vib_04; exE_O2vib_05; ...
	exE_O2vib_06; exE_O2vib_07]*eVtoJ; 

% Neutral densities [m^-3] in same order as the cross sections. Input for degrading electrons 
nN_i = [nN2; nN2; nN2; nN2; nN2; nN2; nN2; nN2; nN2; nN2; nO2; nO2; nO; nO; nO; nO; nO; ...
	nO2; nO2; nO2; nO2; nO2; nO2; nO2]; 

% Neutral densities [m^-3]. Input for electron-meutral collsion frequency
% and eq. 21 in Stubbe 1981
N_n  = nN2 + nO2 + nO;   % Gives same results as Monte Carlo Simulation
N_n2 = nN2 + nO2;        % Gives same results as Monte Carlo Simulation
r_N  = N_n/N_n2;         % Ratio

% Cross sections [m^2]. Energy in must be in [J]. 
E_in = E_bin_mid; 
xs_e_N2vib0_1  = sigma_N2_vib_01( E_in, xs_data{4}(:,1),  xs_data{4}(:,2)); 
xs_e_N2vib0_2  = sigma_N2_vib_02( E_in, xs_data{5}(:,1),  xs_data{5}(:,2));
xs_e_N2vib0_3  = sigma_N2_vib_03( E_in, xs_data{6}(:,1),  xs_data{6}(:,2)); 
xs_e_N2vib0_4  = sigma_N2_vib_04( E_in, xs_data{7}(:,1),  xs_data{7}(:,2));
xs_e_N2vib0_5  = sigma_N2_vib_05( E_in, xs_data{8}(:,1),  xs_data{8}(:,2)); 
xs_e_N2vib0_6  = sigma_N2_vib_06( E_in, xs_data{9}(:,1),  xs_data{9}(:,2));
xs_e_N2vib0_7  = sigma_N2_vib_07( E_in, xs_data{10}(:,1), xs_data{10}(:,2)); 
xs_e_N2vib0_8  = sigma_N2_vib_08( E_in, xs_data{11}(:,1), xs_data{11}(:,2)); 
xs_e_N2vib0_9  = sigma_N2_vib_09( E_in, xs_data{12}(:,1), xs_data{12}(:,2)); 
xs_e_N2vib0_10 = sigma_N2_vib_010(E_in, xs_data{13}(:,1), xs_data{13}(:,2)); 
xs_e_O2a1Dg    = sigma_O2_a1_delta_g(E_in, xs_data{19}(:,1), xs_data{19}(:,2)); 
xs_e_O2b1Sgp   = sigma_O2_b1_delta_g(E_in, xs_data{20}(:,1), xs_data{20}(:,2)); 
xs_e_Ofine_1_0 = sigma_O_fine_str_bell(E_in, xs_data{21}(:,1), xs_data{21}(:,2));  
xs_e_Ofine_2_0 = sigma_O_fine_str_bell(E_in, xs_data{22}(:,1), xs_data{22}(:,2)); 
xs_e_Ofine_2_1 = sigma_O_fine_str_bell(E_in, xs_data{23}(:,1), xs_data{23}(:,2));  
xs_e_O3P1D     = sigma_O_3P_1D(  E_in, xs_data{24}(:,1), xs_data{24}(:,2)); 
xs_e_O3P1S     = sigma_O_3P_1S(  E_in, xs_data{25}(:,1), xs_data{25}(:,2)); 
xs_e_O2vib_01  = sigma_O2_vib_01(E_in, xs_data{26}(:,1), xs_data{26}(:,2)); 
xs_e_O2vib_02  = sigma_O2_vib_02(E_in, xs_data{27}(:,1), xs_data{27}(:,2)); 
xs_e_O2vib_03  = sigma_O2_vib_03(E_in, xs_data{28}(:,1), xs_data{28}(:,2)); 
xs_e_O2vib_04  = sigma_O2_vib_04(E_in, xs_data{29}(:,1), xs_data{29}(:,2)); 
xs_e_O2vib_05  = sigma_O2_vib_05(E_in, xs_data{30}(:,1), xs_data{30}(:,2)); 
xs_e_O2vib_06  = sigma_O2_vib_06(E_in, xs_data{31}(:,1), xs_data{31}(:,2)); 
xs_e_O2vib_07  = sigma_O2_vib_07(E_in, xs_data{32}(:,1), xs_data{32}(:,2)); 
xs_i = [xs_e_N2vib0_1, xs_e_N2vib0_2, xs_e_N2vib0_3, xs_e_N2vib0_4, xs_e_N2vib0_5, ...
	xs_e_N2vib0_6, xs_e_N2vib0_7, xs_e_N2vib0_8, xs_e_N2vib0_9, xs_e_N2vib0_10, ...
	xs_e_O2a1Dg, xs_e_O2b1Sgp, xs_e_Ofine_1_0, xs_e_Ofine_2_0, xs_e_Ofine_2_1, xs_e_O3P1D, xs_e_O3P1S,...
	xs_e_O2vib_01, xs_e_O2vib_02, xs_e_O2vib_03, xs_e_O2vib_04, xs_e_O2vib_05, ...
	xs_e_O2vib_06 xs_e_O2vib_07];

% Set cross sections to zero below threshold for excitation 
for m = 1:length(ex_E)
    [~,ind_trh] = min(abs(E_bin_mid-ex_E(m))); % J
    xs_i(1:ind_trh,m) = 0; 
end 

%% Constants from Stubbe 1981 (which are hopefully correct)
gamma        = -(q_e*E0_field)/me;             % Q*V/(m*kg)
b            = 2.8e-23;                        % [cm*s]
sigma0_B0    = 2.49e-32;                       % [g*cm^4*s^-2]
b_SI         = b*cm_to_m;                      % [m*s]
sigma0_B0_SI = sigma0_B0*g_to_kg*(cm_to_m)^4;  % [kg*m^4*s^-2]
nu           = b_SI*N_n*v_bin_mid.^2;          % Electron-neutral collision frequency [1/s]

%% Solve eq. 21 from Stubbe 1981 by a numerical iteration technique
num_iter = 80; 
f0_iterate = zeros(length(v_bin_mid),num_iter); 
f0_iterate(:,1) = f_max_v; % Form: ne*A*4*pi*v^2*exp(x) and units: [#e/m^3 * s/m]
f0 = zeros(length(v_bin_mid),num_iter); 
G_v = zeros(length(v_bin_mid),num_iter); 

% The iteration 
for j = 1:num_iter

  % Input: Histogram [#e/m^3] per speed bin 
  histogram_in = f0_iterate(:,j)*delta_v; 
  
  % Electrons degrading in energy due to inelastic collisions. Same-ish as eq. 14 in Stubbe 1981 
  [dI_sum,dIe,deg_Ie,~] = degrading_electrons_works3(ex_E,xs_i,E_bins,v_bin_mid,...
						     nN_i,histogram_in); % dI_Sum: [#e/(m^3*s)]
  
  % Convert output from degrading electron function to something somewhat sensible
  df_dt_vib = dI_sum./(4*pi*v_bin_mid.^2*delta_v); % Form: A*exp(x) and units: [#e/(m^3*s) * s^3/m^3]
  
  % This is Eq. 18 in Stubbe 1981: The inelastic collision integral integrated from 0 to u
  G_v(:,j) = cumsum((v_bin_higher.^3/3 - v_bin_lower.^3/3).*df_dt_vib);
  
  % Compute the ugly expression, eq. 21 in Stubbe 1981 
  ledd_oppe = me^2*G_v(:,j); 
  f0_in = (f0_iterate(:,j)./(4*pi*v_bin_mid.^2)); % Form: A*exp(x) and units: [#e/m^3 * s^3/m^3]
  ledd_nede = N_n2*4*sigma0_B0_SI*v_bin_mid.*f0_in; 
  G_v_over_f0 = ledd_oppe./ledd_nede; 
  heating_term = (me^2*b_SI*r_N*gamma.^2)/(24*sigma0_B0_SI) * (v_bin_mid.^3./(w_HF^2 + (nu).^2));
  ugly_exp = (me*v_bin_mid + G_v_over_f0)./(k_b*Tn + heating_term);   
  f0(:,j) = exp(-cumtrapz(v_bin_mid,ugly_exp)); 
  
  % Normalization 
  A = n_e/trapz(v_bin_mid,n_e*4*pi*v_bin_mid.^2.*f0(:,j)); 
  
  % Break after the given number of iterations 
  if j >= num_iter
    break 
  end 
  
  % Distribution: Form: ne*A*4*pi*v^2*exp(x): [#e/m^3 * s/m]
  %f0_iterate(:,j+1) = n_e*A*4*pi*v_bin_mid.^2.*f0(:,j); 
  C = n_e*A*4*pi*v_bin_mid.^2;
  alpha_c = 0.01;
  gamma_c = 0.3; 
  f0_iterate(:,j+1) = (1-alpha_c)*(C.*f0(:,j)).^gamma_c.*f0_iterate(:,j).^(1-gamma_c) + ...
		      alpha_c*f0_iterate(:,j); 
  
end 

%% For saving data
f0_save = f0(:,end); 
A_save = A; 
 
%% Compute temperature of distribution based on Stubbe 1981 version
some_cont = (4*pi*me)/(3*k_b); 
Te_non_max = some_cont*trapz(v_bin_mid,v_bin_mid.^4*A.*f0(:,end));

%% Energy grid with constant delta E

% Bin edges
E_bins_new = (linspace(0,eVtoJ*max_E,res))'; % [J]

% Resolution of energy (constant) [J]
delta_E_new = E_bins_new(end)-E_bins_new(end-1); 

% Energy/speed in the middle fo the bin
E_bin_mid_new = E_bins_new(1:end-1)/2 + E_bins_new(2:end)/2; % [J]

%% NEW distributions for energy grid with constant delta E
f0_A_save_new = interp1(E_bin_mid,A.*f0(:,end),E_bin_mid_new); % Non-Maxwellian 

%% Temperature by gradient
Te_grad = -f0_A_save_new./gradient(f0_A_save_new,delta_E_new*(1/eVtoJ)*(1/K_to_eV)); % delta E in K  

end
