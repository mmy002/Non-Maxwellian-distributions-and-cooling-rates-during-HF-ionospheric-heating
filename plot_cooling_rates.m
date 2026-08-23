%% Analyse and plot data for cooling rates 
% Plot figures for the following paper: non-Maxwellian electron distribution in the D
% region during artificial heating (Paper II): Electron cooling rates

%% Graphical setup
ftsza  = 11.5; 
ftszl  = 14; 
ftszlg = 10.5;
ftszt  = 15; 

%% Some constant 
me                = 9.10938291e-31;   % Electron mass [kg]
eVtoJ             = 1.602176634e-19;  % Convert from [eV] to [J] 
k_b               = 1.380649e-23;     % Boltzmann constant [J/K] 
k_b_eV            = 8.617333e-5;      % Boltzmann constant [eV/K] 
q_e               = 1.602176634e-19;  % Elementary charge of 1 electron [C]
cm_to_m           = 0.01;             % Convert from [cm] to [m]
g_to_kg           = 1/1000;           % Convert from [g] to [kg]
m_neg3_to_cm_neg3 = 1e-6;             % Convert densities from m^-3 to cm^-3
m_to_cm           = 100;              % Convert [m] to [cm]

%% Sort data according to altitude and E-field strength
cd(cooling_rate_data_directory); 
dfiles = dir('*.mat');
load(dfiles(1).name)

% Extract height and electric field from files
get_z0_E0 = cell(numel(dfiles),2); 
for i1 = 1:numel(dfiles)
  
  load(dfiles(i1).name)
  my_filename = dfiles(i1).name; 
  [~,N,~] = fileparts(my_filename); 
  a = strsplit(N,'-');
  get_z0_E0{i1,1} = a{3}; 
  get_z0_E0{i1,2} = a{5}; 
  
end 

% Convert from cell to double
get_z0_E0_new = zeros(size(get_z0_E0)); 
for k = 1:size(get_z0_E0,1)
  get_z0_E0_new(k,1) = str2double(get_z0_E0{k,1}); 
  get_z0_E0_new(k,2) = str2double(get_z0_E0{k,2});
end 

% Sort files 
[sort_dfiles,ind_sort_files] = sortrows(get_z0_E0_new); 
dfiles_cell = struct2cell(dfiles)'; 
dfiles_new1 = cell(size(dfiles,1),6); 
for i1 = 1:numel(dfiles)
  index = find(i1 == ind_sort_files); 
  dfiles_new1{index,1} = dfiles_cell{i1,1}; 
  dfiles_new1{index,2} = dfiles_cell{i1,2}; 
  dfiles_new1{index,3} = dfiles_cell{i1,3}; 
  dfiles_new1{index,4} = dfiles_cell{i1,4}; 
  dfiles_new1{index,5} = dfiles_cell{i1,5}; 
  dfiles_new1{index,6} = dfiles_cell{i1,6};  
end 

% Convert from cell to struct
fields = ["name", "folder", "date", "bytes", "isdir", "datenum"];
for i1 = 1:numel(dfiles)
  dfiles_new2 = cell2struct(dfiles_new1,fields,2);
end 

%% Find height z0 and electric field strength E0 in data 
z0 = unique(sort_dfiles(:,1)); 
E0_field_height = zeros(length(z0),size(sort_dfiles,1)/length(z0)); 
for l = 1:length(z0)
  for k = 1:size(sort_dfiles,1)/length(z0)
    n = l-1; 
    m = k + size(sort_dfiles,1)/length(z0)*n; 
    E0_field_height(l,k) = sort_dfiles(m,2); 
  end 
end

%% Index at different heights and electric field: 
% row are height index
% columns are electric field index
ind = zeros(length(E0_field_height),length(z0))'; 
v = 0; 
for j = 1:length(z0)
  a = 1+v; 
  ind(j,:) = a:a+length(E0_field_height)-1; 
  v = v+length(E0_field_height); 
end 

%% rehape height vector z0 to match the number of files 
z_intm = zeros(size(ind)); 
for m = 1:length(z0)
  z_intm(m,:) = z0(m); 
end
z = reshape(z_intm',size(dfiles)); 

%% Load and process data for cross sections
% With more O2 vib. cross sections and new cross section from Bell for fs O 
% and new cross sections for O2 and N2 rot. 
xs_data = get_data_even_more_xs(); 

%% Load data (go first to the folder containing the data...)
load(dfiles_new2(1).name)

% Non-Maxwellian electron temperature 
Te_NM_2nd = zeros(1,length(z0)*length(E0_field_height));
Te_NM_L = zeros(1,length(z0)*length(E0_field_height));

% Q
Q_0v_non_max_N2      = zeros(1,length(z0)*length(E0_field_height));
Q_0v_maxw_N2_2nd     = zeros(1,length(z0)*length(E0_field_height));
Q_0v_maxw_N2_L       = zeros(1,length(z0)*length(E0_field_height));
Q_0v_non_max_O2      = zeros(1,length(z0)*length(E0_field_height));
Q_0v_maxw_O2_2nd     = zeros(1,length(z0)*length(E0_field_height));
Q_0v_maxw_O2_L       = zeros(1,length(z0)*length(E0_field_height));
Q_non_max_O_fs       = zeros(1,length(z0)*length(E0_field_height));
Q_maxw_O_fs_2nd      = zeros(1,length(z0)*length(E0_field_height));
Q_maxw_O_fs_L        = zeros(1,length(z0)*length(E0_field_height));
Q_non_max_a1delta_O2 = zeros(1,length(z0)*length(E0_field_height));
Q_maxw_a1delta_O2_L  = zeros(1,length(z0)*length(E0_field_height));
Q_non_max_b1sigma_O2 = zeros(1,length(z0)*length(E0_field_height));
Q_maxw_b1sigma_O2_L  = zeros(1,length(z0)*length(E0_field_height));
Q_non_max_rot_N2     = zeros(1,length(z0)*length(E0_field_height));
Q_maxw_rot_N2_L      = zeros(1,length(z0)*length(E0_field_height));
Q_non_max_rot_O2     = zeros(1,length(z0)*length(E0_field_height));
Q_maxw_rot_O2_L      = zeros(1,length(z0)*length(E0_field_height));
Q_non_max_el         = zeros(1,length(z0)*length(E0_field_height));
Q_maxw_el_L          = zeros(1,length(z0)*length(E0_field_height));

% Maxwellian cooling rates from Pavlov tables 
L_oxygen_10_pavlov_tables_2nd  = zeros(1,length(z0)*length(E0_field_height));
L_oxygen_20_pavlov_tables_2nd  = zeros(1,length(z0)*length(E0_field_height));
L_oxygen_21_pavlov_tables_2nd  = zeros(1,length(z0)*length(E0_field_height));
L_O_fs_pavlov_tables_2nd       = zeros(1,length(z0)*length(E0_field_height));
L_a1delta_O2_pavlov_tables_2nd = zeros(1,length(z0)*length(E0_field_height));
L_b1sigma_O2_pavlov_tables_2nd = zeros(1,length(z0)*length(E0_field_height));
L_rot_N2_pavlov_tables_2nd     = zeros(1,length(z0)*length(E0_field_height));
L_rot_O2_pavlov_tables_2nd     = zeros(1,length(z0)*length(E0_field_height));
L_elas_pavlov_tables_2nd       = zeros(1,length(z0)*length(E0_field_height));
L_oxygen_10_pavlov_tables_L    = zeros(1,length(z0)*length(E0_field_height));
L_oxygen_20_pavlov_tables_L    = zeros(1,length(z0)*length(E0_field_height));
L_oxygen_21_pavlov_tables_L    = zeros(1,length(z0)*length(E0_field_height));
L_O_fs_pavlov_tables_L         = zeros(1,length(z0)*length(E0_field_height));
L_a1delta_O2_pavlov_tables_L   = zeros(1,length(z0)*length(E0_field_height));
L_b1sigma_O2_pavlov_tables_L   = zeros(1,length(z0)*length(E0_field_height));
L_vib_N2_pavlov_tables_L       = zeros(1,length(z0)*length(E0_field_height));
L_vib_O2_pavlov_tables_L       = zeros(1,length(z0)*length(E0_field_height));
L_rot_N2_pavlov_tables_L       = zeros(1,length(z0)*length(E0_field_height));
L_rot_O2_pavlov_tables_L       = zeros(1,length(z0)*length(E0_field_height));
L_elas_pavlov_tables_L         = zeros(1,length(z0)*length(E0_field_height));

% Cooling rates for N2 vib.
L_vib_N2_nonmax   = zeros(1,length(z0)*length(E0_field_height));
L_vib_N2_maxw_2nd = zeros(1,length(z0)*length(E0_field_height));
L_vib_N2_maxw_L   = zeros(1,length(z0)*length(E0_field_height));

% Cooling rates for O2 vib.
L_vib_O2_nonmax   = zeros(1,length(z0)*length(E0_field_height));
L_vib_O2_maxw_2nd = zeros(1,length(z0)*length(E0_field_height));
L_vib_O2_maxw_L   = zeros(1,length(z0)*length(E0_field_height));

% Cooling rates from inelastic collision integral from non-maxwellia distribution 
L_N2_vib_dIe    = zeros(1,length(z0)*length(E0_field));
L_Ofine_1_0_dIe = zeros(1,length(z0)*length(E0_field_height));
L_Ofine_2_0_dIe = zeros(1,length(z0)*length(E0_field_height));
L_Ofine_2_1_dIe = zeros(1,length(z0)*length(E0_field_height));
L_non_max_O_fs  = zeros(1,length(z0)*length(E0_field_height));
L_O2a1Dg_dIe    = zeros(1,length(z0)*length(E0_field_height));
L_O2b1Sgp_dIe   = zeros(1,length(z0)*length(E0_field_height));
L_O3P1D_dIe     = zeros(1,length(z0)*length(E0_field_height));
L_O3P1S_dIe     = zeros(1,length(z0)*length(E0_field_height));
L_O2_vib_dIe    = zeros(1,length(z0)*length(E0_field));

% Calculated non-Maxwellian cooling rates 
L_N2_rot = zeros(1,length(z0)*length(E0_field_height));
L_O2_rot = zeros(1,length(z0)*length(E0_field_height));
L_el     = zeros(1,length(z0)*length(E0_field_height));

% Absorped power = cooling rates 
L = zeros(1,length(z0)*length(E0_field_height));

for i1 = 1:numel(dfiles)
  
  load(dfiles_new2(i1).name)
  
  % Non-Maxwellian temperature [K]
  Te_NM_2nd(i1) = Te_nm_2nd;
  Te_NM_L(i1) = Te_nm_L; 
  
  % Call MATLAB MSIS [1/m^3]
  [rho,~] = neutral_output_MSIS(z(i1)*1e3); 
  
  % Neutral densities MATLAB MSIS
  nO = rho(:,2)';   % [1/m^3] 
  nN2 = rho(:,3)';    % [1/m^3] 
  nO2 = rho(:,4)';    % [1/m^3]
  
  % Electron density as a function of height [1/m^3]
  N_e = n_e; 
  
  % Neutral temperature [K]
  T_n = Tn; 
  
  % Q_0v [ev*cm^3/s]
  Q_0v_maxw_N2_2nd(i1) = sum(Q_maxw_2nd{1,1});   % OBS: newer cross sections than Pavlov N2 vib. !!!! 
  Q_0v_maxw_N2_L(i1) = sum(Q_maxw_L{1,1});       % OBS: newer cross sections than Pavlov N2 vib. !!!! 
  Q_0v_maxw_O2_2nd(i1) = sum(Q_maxw_2nd{1,2}); 
  Q_0v_maxw_O2_L(i1) = sum(Q_maxw_L{1,2});  
  
  % Absorped power=sum of cooling rates
  L(i1) = L_loop; % [J m^-3 s^-1]
  
  % Pavlov tables cooling rates [eV/(cm^3*s)]
  % At temperature T_{2nd}
  [~,~,L_oxygen_10_pavlov_tables_2nd(i1),L_oxygen_20_pavlov_tables_2nd(i1),L_oxygen_21_pavlov_tables_2nd(i1),...
   L_rot_N2_pavlov_tables_2nd(i1),L_rot_O2_pavlov_tables_2nd(i1),L_elas_pavlov_tables_2nd(i1),...
   L_a1delta_O2_pavlov_tables_2nd(i1),L_b1sigma_O2_pavlov_tables_2nd(i1)] ...
       = cooling_rates_tables_pavlov2(Te_NM_2nd(i1),Tn,N_e,nN2,nO2,nO);
  L_O_fs_pavlov_tables_2nd(i1) = L_oxygen_10_pavlov_tables_2nd(i1) + ...
				 L_oxygen_20_pavlov_tables_2nd(i1) + ...
				 L_oxygen_21_pavlov_tables_2nd(i1);
  % At temperature T_L
  [L_vib_N2_pavlov_tables_L(i1),L_vib_O2_pavlov_tables_L(i1),L_oxygen_10_pavlov_tables_L(i1),...
   L_oxygen_20_pavlov_tables_L(i1),L_oxygen_21_pavlov_tables_L(i1),...
   L_rot_N2_pavlov_tables_L(i1),L_rot_O2_pavlov_tables_L(i1),L_elas_pavlov_tables_L(i1),...
   L_a1delta_O2_pavlov_tables_L(i1),L_b1sigma_O2_pavlov_tables_L(i1)] ...
    = cooling_rates_tables_pavlov2(Te_NM_L(i1),Tn,N_e,nN2,nO2,nO);
  L_O_fs_pavlov_tables_L(i1) = L_oxygen_10_pavlov_tables_L(i1) + ...
			       L_oxygen_20_pavlov_tables_L(i1) + ...
			       L_oxygen_21_pavlov_tables_L(i1);

  % Cooling rates [ev*cm^-3/s]
  L_vib_N2_nonmax(i1)   = L_non_max(1,1);
  L_vib_N2_maxw_2nd(i1) = L_maxw_2nd(1,1); % OBS: newer cross sections than Pavlov N2 vib. !!!!
  L_vib_N2_maxw_L(i1)   = L_maxw_L(1,1);   % OBS: newer cross sections than Pavlov N2 vib. !!!!
  L_vib_O2_nonmax(i1)   = L_non_max(2,1);
  L_vib_O2_maxw_2nd(i1) = L_maxw_2nd(2,1);
  L_vib_O2_maxw_L(i1)   = L_maxw_L(2,1);

  % Degrading electrons due to collisions: Sum of loss and gain [#e/(m^3*s)]
  sum_dIe_N2_vib    = sum(dIe(:,1:10),2) + sum(deg_Ie(:,1:10),2);  
  sum_dIe_Ofine_1_0 = sum(dIe(:,13),2) + sum(deg_Ie(:,13),2); 
  sum_dIe_Ofine_2_0 = sum(dIe(:,14),2) + sum(deg_Ie(:,14),2); 
  sum_dIe_Ofine_2_1 = sum(dIe(:,15),2) + sum(deg_Ie(:,15),2); 
  sum_dIe_O2a1Dg    = sum(dIe(:,11),2) + sum(deg_Ie(:,11),2); 
  sum_dIe_O2b1Sgp   = sum(dIe(:,12),2) + sum(deg_Ie(:,12),2); 
  sum_dIe_O3P1D     = sum(dIe(:,16),2) + sum(deg_Ie(:,16),2); 
  sum_dIe_O3P1S     = sum(dIe(:,17),2) + sum(deg_Ie(:,17),2); 
  sum_dIe_O2_vib    = sum(dIe(:,18:end),2) + sum(deg_Ie(:,18:end),2);
  
  % Energy/speed in the middle of the bin
  v_bin_mid = v_bins_stubbe(1:end-1)/2 + v_bins_stubbe(2:end)/2; % [m/s]
  E_bin_mid = (me*v_bin_mid.^2)/2; % [J]
  E_bin_mid_eV = (me*v_bin_mid.^2)/(2*eVtoJ); % [eV]
  
  % Cooling rates from inelastic collision integral E(eV)*histogram [ev/(cm^3*s)]
  L_N2_vib_dIe(i1)    = sum(-E_bin_mid_eV.*sum_dIe_N2_vib)*m_neg3_to_cm_neg3; 
  L_Ofine_1_0_dIe(i1) = sum(-E_bin_mid_eV.*sum_dIe_Ofine_1_0)*m_neg3_to_cm_neg3; 
  L_Ofine_2_0_dIe(i1) = sum(-E_bin_mid_eV.*sum_dIe_Ofine_2_0)*m_neg3_to_cm_neg3; 
  L_Ofine_2_1_dIe(i1) = sum(-E_bin_mid_eV.*sum_dIe_Ofine_2_1)*m_neg3_to_cm_neg3; 
  L_non_max_O_fs(i1)  = L_Ofine_1_0_dIe(i1)+L_Ofine_2_0_dIe(i1)+L_Ofine_2_1_dIe(i1);
  L_O2a1Dg_dIe(i1)    = sum(-E_bin_mid_eV.*sum_dIe_O2a1Dg)*m_neg3_to_cm_neg3; 
  L_O2b1Sgp_dIe(i1)   = sum(-E_bin_mid_eV.*sum_dIe_O2b1Sgp)*m_neg3_to_cm_neg3; 
  L_O3P1D_dIe(i1)     = sum(-E_bin_mid_eV.*sum_dIe_O3P1D)*m_neg3_to_cm_neg3; 
  L_O3P1S_dIe(i1)     = sum(-E_bin_mid_eV.*sum_dIe_O3P1S)*m_neg3_to_cm_neg3;
  L_O2_vib_dIe(i1)    = sum(-E_bin_mid_eV.*sum_dIe_O2_vib)*m_neg3_to_cm_neg3; 
  
  % Cooling rates for rotational and elastic [ev/(cm^3*s)] 
  [L_N2_rot(i1),L_O2_rot(i1),L_el(i1)] = L_rot_el_new(v_bins_stubbe,...
						      f0_exp_stubbe,...
						      A,...
						      xs_data,...
						      N_e,...
						      nN2,...
						      nO2,...
						      nO,...
						      Te_NM_2nd(i1),...
						      Tn); 
  L_N2_rot(i1) = L_N2_rot(i1)*m_neg3_to_cm_neg3*(1/eVtoJ); 
  L_O2_rot(i1) = L_O2_rot(i1)*m_neg3_to_cm_neg3*(1/eVtoJ); 
  L_el(i1)     = L_el(i1)*m_neg3_to_cm_neg3*(1/eVtoJ);  
  
  % Q [ev*cm^3/s]
  Q_0v_non_max_N2(i1)      = L_N2_vib_dIe(i1)/(N_e*m_neg3_to_cm_neg3*nN2*m_neg3_to_cm_neg3); 
  Q_0v_non_max_O2(i1)      = L_O2_vib_dIe(i1)/(N_e*m_neg3_to_cm_neg3*nO2*m_neg3_to_cm_neg3); 
  Q_non_max_O_fs(i1)       = (L_Ofine_1_0_dIe(i1)+L_Ofine_2_0_dIe(i1)+L_Ofine_2_1_dIe(i1))/(N_e*m_neg3_to_cm_neg3*nO*m_neg3_to_cm_neg3); 
  Q_maxw_O_fs_2nd(i1)      = L_O_fs_pavlov_tables_2nd(i1)/(N_e*m_neg3_to_cm_neg3*nO*m_neg3_to_cm_neg3); 
  Q_maxw_O_fs_L(i1)        = L_O_fs_pavlov_tables_L(i1)/(N_e*m_neg3_to_cm_neg3*nO*m_neg3_to_cm_neg3); 
  Q_non_max_a1delta_O2(i1) = L_O2a1Dg_dIe(i1)/(N_e*m_neg3_to_cm_neg3*nO2*m_neg3_to_cm_neg3); 
  Q_maxw_a1delta_O2_L(i1)  = L_a1delta_O2_pavlov_tables_L(i1)/(N_e*m_neg3_to_cm_neg3*nO2*m_neg3_to_cm_neg3); 
  Q_non_max_b1sigma_O2(i1) = L_O2b1Sgp_dIe(i1)/(N_e*m_neg3_to_cm_neg3*nO2*m_neg3_to_cm_neg3); 
  Q_maxw_b1sigma_O2_L(i1)  = L_b1sigma_O2_pavlov_tables_L(i1)/(N_e*m_neg3_to_cm_neg3*nO2*m_neg3_to_cm_neg3); 
  Q_non_max_rot_N2(i1)     = L_N2_rot(i1)/(N_e*m_neg3_to_cm_neg3*nN2*m_neg3_to_cm_neg3); 
  Q_maxw_rot_N2_L(i1)      = L_rot_N2_pavlov_tables_L(i1)/(N_e*m_neg3_to_cm_neg3*nN2*m_neg3_to_cm_neg3);
  Q_non_max_rot_O2(i1)     = L_O2_rot(i1)/(N_e*m_neg3_to_cm_neg3*nO2*m_neg3_to_cm_neg3); 
  Q_maxw_rot_O2_L(i1)      = L_rot_O2_pavlov_tables_L(i1)/(N_e*m_neg3_to_cm_neg3*nO2*m_neg3_to_cm_neg3);
  Q_non_max_el(i1)         = L_el(i1)/(N_e*m_neg3_to_cm_neg3*(nN2+nO2+nO)*m_neg3_to_cm_neg3); 
  Q_maxw_el_L(i1)          = L_elas_pavlov_tables_L(i1)/(N_e*m_neg3_to_cm_neg3*(nN2+nO2+nO)*m_neg3_to_cm_neg3);
  
end 

%% Choose electric fields to plot for the different heights
pick_E0_70_km  = 1:size(E0_field_height,2)-1; 
pick_E0_80_km  = 1:size(E0_field_height,2)-1; 
pick_E0_90_km  = 2:size(E0_field_height,2); 
pick_E0_100_km = 2:size(E0_field_height,2); 

%% N2: M->2nd

figure

subplot(2,2,1)
yyaxis left
plot(E0_field_height(1,pick_E0_70_km),log10(Q_0v_non_max_N2(ind(1,pick_E0_70_km))),'b','LineWidth',2) 
hold on;
plot(E0_field_height(1,pick_E0_70_km),log10(Q_0v_maxw_N2_2nd(ind(1,pick_E0_70_km))),'b-.','LineWidth',2)
hold off; 
ax = gca; 
ax.YColor = 'b'; 
ylabel('log Q [eV cm^3 s^{-1}]','FontSize',ftszl)
ylim([-15 -8])
yticks([-15 -14 -13 -12 -11 -10 -9 -8])
title('70 km','FontSize',ftszt) 
yyaxis right
plot(E0_field_height(1,pick_E0_70_km),...
     Te_NM_2nd(ind(1,pick_E0_70_km)),...
     'r--','LineWidth',2) 
ylabel('Te [K]')
ylim([0 6000])
xlim([0 4])
legend('Q_{NM}','Q_{M}','T_{2nd}',...
       'Location','northwest',...
       'FontSize',ftszlg)
ax = gca; 
ax.YColor = 'r';
set(gca,...
    'fontsize',ftsza,...
    'box','off',...
    'tickdir','both',...
    'TickLength',[0.02 0.01])

subplot(2,2,2)
yyaxis left
plot(E0_field_height(2,pick_E0_80_km),log10(Q_0v_non_max_N2(ind(2,pick_E0_80_km))),'b','LineWidth',2) 
hold on;
plot(E0_field_height(2,pick_E0_80_km),log10(Q_0v_maxw_N2_2nd(ind(2,pick_E0_80_km))),'b-.','LineWidth',2)
hold off; 
ax = gca; 
ax.YColor = 'b'; 
ylabel('log Q [eV cm^3 s^{-1}]','FontSize',ftszl)
ylim([-15 -8])
yticks([-15 -14 -13 -12 -11 -10 -9 -8])
title('80 km','FontSize',ftszt) 
yyaxis right
plot(E0_field_height(2,pick_E0_80_km),...
     Te_NM_2nd(ind(2,pick_E0_80_km)),...
     'r--','LineWidth',2) 
legend()
ylabel('Te [K]')
ylim([0 6000])
xlim([0 4])
legend('Q_{NM}','Q_{M}','T_{2nd}',...
       'Location','northwest',...
       'FontSize',ftszlg)
ax = gca; 
ax.YColor = 'r';
set(gca,...
    'fontsize',ftsza,...
    'box','off',...
    'tickdir','both',...
    'TickLength',[0.02 0.01])

subplot(2,2,3)
yyaxis left
plot(E0_field_height(3,pick_E0_90_km),log10(Q_0v_non_max_N2(ind(3,pick_E0_90_km))),'b','LineWidth',2) 
hold on;
plot(E0_field_height(3,pick_E0_90_km),log10(Q_0v_maxw_N2_2nd(ind(3,pick_E0_90_km))),'b-.','LineWidth',2)
hold off; 
ax = gca; 
ax.YColor = 'b'; 
ylabel('log Q [eV cm^3 s^{-1}]','FontSize',ftszl)
xlabel('E_0 [V/m]','FontSize',ftszl)
ylim([-15 -8])
yticks([-15 -14 -13 -12 -11 -10 -9 -8])
title('90 km','FontSize',ftszt) 
yyaxis right
plot(E0_field_height(3,pick_E0_90_km),...
     Te_NM_2nd(ind(3,pick_E0_90_km)),...
     'r--','LineWidth',2) 
legend()
ylabel('Te [K]')
ylim([0 6000])
xlim([0 4])
legend('Q_{NM}','Q_{M}','T_{2nd}',...
       'Location','northwest',...
       'FontSize',ftszlg)
ax = gca; 
ax.YColor = 'r';
set(gca,...
    'fontsize',ftsza,...
    'box','off',...
    'tickdir','both',...
    'TickLength',[0.02 0.01])

subplot(2,2,4)
yyaxis left
plot(E0_field_height(4,pick_E0_100_km),log10(Q_0v_non_max_N2(ind(4,pick_E0_100_km))),'b','LineWidth',2) 
hold on;
plot(E0_field_height(4,pick_E0_100_km),log10(Q_0v_maxw_N2_2nd(ind(4,pick_E0_100_km))),'b-.','LineWidth',2)
hold off; 
ax = gca; 
ax.YColor = 'b'; 
ylabel('log Q [eV cm^3 s^{-1}]','FontSize',ftszl)
xlabel('E_0 [V/m]','FontSize',ftszl)
ylim([-15 -8])
yticks([-15 -14 -13 -12 -11 -10 -9 -8])
title('100 km','FontSize',ftszt) 
yyaxis right
plot(E0_field_height(4,pick_E0_100_km),Te_NM_2nd(ind(4,pick_E0_100_km)),'r--','LineWidth',2) 
legend()
ylabel('Te [K]')
ylim([0 6000])
xlim([0 4])
legend('Q_{NM}','Q_{M}','T_{2nd}',...
       'Location','northwest',...
       'FontSize',ftszlg)
ax = gca; 
ax.YColor = 'r';
set(gca,...
    'fontsize',ftsza,...
    'box','off',...
    'tickdir','both',...
    'TickLength',[0.02 0.01]) 

sgtitle('Vibrational excitation of molecular nitrogen')

%% O2: M->2nd 

figure

subplot(2,2,1)
yyaxis left
plot(E0_field_height(1,pick_E0_70_km),log10(Q_0v_non_max_O2(ind(1,pick_E0_70_km))),'b','LineWidth',2) 
hold on;
plot(E0_field_height(1,pick_E0_70_km),log10(Q_0v_maxw_O2_2nd(ind(1,pick_E0_70_km))),'b-.','LineWidth',2)
hold off; 
ax = gca; 
ax.YColor = 'b'; 
ylabel('log Q [eV cm^3 s^{-1}]','FontSize',ftszl)
ylim([-13 -8])
yticks([-13 -12 -11 -10 -9 -8])
title('70 km','FontSize',ftszt) 
yyaxis right
plot(E0_field_height(1,pick_E0_70_km),Te_NM_2nd(ind(1,pick_E0_70_km)),'r--','LineWidth',2) 
legend()
ylabel('Te [K]')
ylim([0 6000])
xlim([0 4])
legend('Q_{NM}','Q_{M}','T_{2nd}',...
       'Location','northwest',...
       'FontSize',ftszlg)
ax = gca; 
ax.YColor = 'r';
set(gca,...
    'fontsize',ftsza,...
    'box','off',...
    'tickdir','both',...
    'TickLength',[0.02 0.01])

subplot(2,2,2)
yyaxis left
plot(E0_field_height(2,pick_E0_80_km),log10(Q_0v_non_max_O2(ind(2,pick_E0_80_km))),'b','LineWidth',2) 
hold on;
plot(E0_field_height(2,pick_E0_80_km),log10(Q_0v_maxw_O2_2nd(ind(2,pick_E0_80_km))),'b-.','LineWidth',2)
hold off; 
ax = gca; 
ax.YColor = 'b'; 
ylabel('log Q [eV cm^3 s^{-1}]','FontSize',ftszl)
ylim([-13 -8])
yticks([-13 -12 -11 -10 -9 -8])
title('80 km','FontSize',ftszt) 
yyaxis right
plot(E0_field_height(2,pick_E0_80_km),Te_NM_2nd(ind(2,pick_E0_80_km)),'r--','LineWidth',2) 
legend()
ylabel('Te [K]')
ylim([0 6000])
xlim([0 4])
legend('Q_{NM}','Q_{M}','T_{2nd}',...
       'Location','northwest',...
       'FontSize',ftszlg)
ax = gca; 
ax.YColor = 'r';
set(gca,...
    'fontsize',ftsza,...
    'box','off',...
    'tickdir','both',...
    'TickLength',[0.02 0.01])

subplot(2,2,3)
yyaxis left
plot(E0_field_height(3,pick_E0_90_km),log10(Q_0v_non_max_O2(ind(3,pick_E0_90_km))),'b','LineWidth',2) 
hold on;
plot(E0_field_height(3,pick_E0_90_km),log10(Q_0v_maxw_O2_2nd(ind(3,pick_E0_90_km))),'b-.','LineWidth',2)
hold off; 
ax = gca; 
ax.YColor = 'b'; 
ylabel('log Q_ [eV cm^3 s^{-1}]','FontSize',ftszl)
xlabel('E_0 [V/m]','FontSize',ftszl)
ylim([-13 -8])
yticks([-13 -12 -11 -10 -9 -8])
title('90 km','FontSize',ftszt) 
yyaxis right
plot(E0_field_height(3,pick_E0_90_km),...
     Te_NM_2nd(ind(3,pick_E0_90_km)),...
     'r--','LineWidth',2) 
legend()
ylabel('Te [K]')
ylim([0 6000])
xlim([0 4])
legend('Q_{NM}','Q_{M}','T_{2nd}',...
       'Location','northwest',...
       'FontSize',ftszlg)
ax = gca; 
ax.YColor = 'r';
set(gca,...
    'fontsize',ftsza,...
    'box','off',...
    'tickdir','both',...
    'TickLength',[0.02 0.01])

subplot(2,2,4)
yyaxis left
plot(E0_field_height(4,pick_E0_100_km),log10(Q_0v_non_max_O2(ind(4,pick_E0_100_km))),'b','LineWidth',2) 
hold on;
plot(E0_field_height(4,pick_E0_100_km),log10(Q_0v_maxw_O2_2nd(ind(4,pick_E0_100_km))),'b-.','LineWidth',2)
hold off; 
ax = gca; 
ax.YColor = 'b'; 
ylabel('log Q [eV cm^3 s^{-1}]','FontSize',ftszl)
xlabel('E_0 [V/m]','FontSize',ftszl)
ylim([-13 -8])
yticks([-13 -12 -11 -10 -9 -8])
title('100 km','FontSize',ftszt) 
yyaxis right
plot(E0_field_height(4,pick_E0_100_km),Te_NM_2nd(ind(4,pick_E0_100_km)),'r--','LineWidth',2) 
legend()
ylabel('Te [K]')
ylim([0 6000])
xlim([0 4])
legend('Q_{NM}','Q_{M}','T_{2nd}',...
       'Location','northwest',...
       'FontSize',ftszlg)
ax = gca; 
ax.YColor = 'r';
set(gca,...
    'fontsize',ftsza,...
    'box','off',...
    'tickdir','both',...
    'TickLength',[0.02 0.01]) 

sgtitle('Vibrational excitation of molecular oxygen')

%% Fs O: M->2nd 

figure 

subplot(2,2,1)
yyaxis left
plot(E0_field_height(1,pick_E0_70_km),log10(Q_non_max_O_fs(ind(1,pick_E0_70_km))),'b','LineWidth',2) 
hold on;
plot(E0_field_height(1,pick_E0_70_km),log10(Q_maxw_O_fs_2nd(ind(1,pick_E0_70_km))),'b-.','LineWidth',2.0)
hold off; 
ax = gca; 
ax.YColor = 'b'; 
ylabel('log Q [eV cm^3 s^{-1}]','FontSize',ftszl)
ylim([-11.6 -10.8])
title('70 km','FontSize',ftszt) 
yyaxis right
plot(E0_field_height(1,pick_E0_70_km),...
     Te_NM_2nd(ind(1,pick_E0_70_km)),...
     'r--','LineWidth',2) 
legend()
ylabel('Te [K]')
ylim([0 6000])
xlim([0 4])
legend('Q_{NM}','Q_{M}','T_{2nd}',...
       'Location','northwest',...
       'FontSize',ftszlg)
ax = gca; 
ax.YColor = 'r';
set(gca,...
    'fontsize',ftsza,...
    'box','off',...
    'tickdir','both',...
    'TickLength',[0.02 0.01]) 

subplot(2,2,2)
yyaxis left
plot(E0_field_height(2,pick_E0_80_km),log10(Q_non_max_O_fs(ind(2,pick_E0_80_km))),'b','LineWidth',2) 
hold on;
plot(E0_field_height(2,pick_E0_80_km),log10(Q_maxw_O_fs_2nd(ind(2,pick_E0_80_km))),'b-.','LineWidth',2.0)
hold off; 
ax = gca; 
ax.YColor = 'b'; 
ylabel('log Q [eV cm^3 s^{-1}]','FontSize',ftszl)
ylim([-11.6 -10.8])
title('80 km','FontSize',ftszt) 
yyaxis right
plot(E0_field_height(2,pick_E0_80_km),Te_NM_2nd(ind(2,pick_E0_80_km)),'r--','LineWidth',2) 
legend()
ylabel('Te [K]')
ylim([0 6000])
xlim([0 4])
legend('Q_{NM}','Q_{M}','T_{2nd}',...
       'Location','northwest',...
       'FontSize',ftszlg)
ax = gca; 
ax.YColor = 'r';
set(gca,...
    'fontsize',ftsza,...
    'box','off',...
    'tickdir','both',...
    'TickLength',[0.02 0.01]) 

subplot(2,2,3)
yyaxis left
plot(E0_field_height(3,pick_E0_90_km),log10(Q_non_max_O_fs(ind(3,pick_E0_90_km))),'b','LineWidth',2) 
hold on;
plot(E0_field_height(3,pick_E0_90_km),log10(Q_maxw_O_fs_2nd(ind(3,pick_E0_90_km))),'b-.','LineWidth',2.0)
hold off; 
ax = gca; 
ax.YColor = 'b'; 
ylabel('log Q [eV cm^3 s^{-1}]','FontSize',ftszl)
xlabel('E_0 [V/m]','FontSize',ftszl)
ylim([-11.6 -10.8])
title('90 km','FontSize',ftszt) 
yyaxis right
plot(E0_field_height(3,pick_E0_90_km),Te_NM_2nd(ind(3,pick_E0_90_km)),'r--','LineWidth',2) 
ylabel('Te [K]')
ylim([0 6000])
xlim([0 4])
legend('Q_{NM}','Q_{M}','T_{2nd}',...
       'Location','northwest',...
       'FontSize',ftszlg)
ax = gca; 
ax.YColor = 'r';
set(gca,...
    'fontsize',ftsza,...
    'box','off',...
    'tickdir','both',...
    'TickLength',[0.02 0.01]) 

subplot(2,2,4)
yyaxis left
plot(E0_field_height(4,pick_E0_100_km),log10(Q_non_max_O_fs(ind(4,pick_E0_100_km))),'b','LineWidth',2) 
hold on;
plot(E0_field_height(4,pick_E0_100_km),log10(Q_maxw_O_fs_2nd(ind(4,pick_E0_100_km))),'b-.','LineWidth',2.0)
hold off; 
ax = gca; 
ax.YColor = 'b'; 
ylabel('log Q [eV cm^3 s^{-1}]','FontSize',ftszl)
xlabel('E_0 [V/m]','FontSize',ftszl)
ylim([-11.6 -10.8])
title('100 km','FontSize',ftszt) 
yyaxis right
plot(E0_field_height(4,pick_E0_100_km),Te_NM_2nd(ind(4,pick_E0_100_km)),'r--','LineWidth',2) 
legend()
ylabel('Te [K]')
ylim([0 6000])
xlim([0 4])
legend('Q_{NM}','Q_{M}','T_{2nd}',...
       'Location','northwest',...
       'FontSize',ftszlg)
ax = gca; 
ax.YColor = 'r';
set(gca,...
    'fontsize',ftsza,...
    'box','off',...
    'tickdir','both',...
    'TickLength',[0.02 0.01]) 

sgtitle('Excitation of fine structures in atomic oxygen')

%% N2: M->L

figure

subplot(2,2,1)
yyaxis left
plot(E0_field_height(1,pick_E0_70_km),...
     log10(Q_0v_non_max_N2(ind(1,pick_E0_70_km))),...
     'b','LineWidth',2) 
hold on;
plot(E0_field_height(1,pick_E0_70_km),log10(Q_0v_maxw_N2_L(ind(1,pick_E0_70_km))),'b-.','LineWidth',2)
hold off; 
ax = gca; 
ax.YColor = 'b'; 
ylabel('log Q [eV cm^3 s^{-1}]','FontSize',ftszl)
ylim([-15 -9])
yticks([-15 -14 -13 -12 -11 -10 -9])
title('70 km','FontSize',ftszt) 
yyaxis right
plot(E0_field_height(1,pick_E0_70_km),...
     Te_NM_L(ind(1,pick_E0_70_km)),...
     'r--','LineWidth',2) 
ylabel('Te [K]')
ylim([0 4000])
xlim([0 4])
legend('Q_{NM}','Q_{M}','T_L',...
       'Location','northwest',...
       'FontSize',ftszlg)
ax = gca; 
ax.YColor = 'r';
set(gca,...
    'fontsize',ftsza,...
    'box','off',...
    'tickdir','both',...
    'TickLength',[0.02 0.01])

subplot(2,2,2)
yyaxis left
plot(E0_field_height(2,pick_E0_80_km),log10(Q_0v_non_max_N2(ind(2,pick_E0_80_km))),'b','LineWidth',2) 
hold on;
plot(E0_field_height(2,pick_E0_80_km),log10(Q_0v_maxw_N2_L(ind(2,pick_E0_80_km))),'b-.','LineWidth',2)
hold off; 
ax = gca; 
ax.YColor = 'b'; 
ylabel('log Q [eV cm^3 s^{-1}]','FontSize',ftszl)
ylim([-15 -9])
yticks([-15 -14 -13 -12 -11 -10 -9])
title('80 km','FontSize',ftszt) 
yyaxis right
plot(E0_field_height(2,pick_E0_80_km),...
     Te_NM_L(ind(2,pick_E0_80_km)),...
     'r--','LineWidth',2) 
legend()
ylabel('Te [K]')
ylim([0 4000])
xlim([0 4])
legend('Q_{NM}','Q_{M}','T_L',...
       'Location','northwest',...
       'FontSize',ftszlg)
ax = gca; 
ax.YColor = 'r';
set(gca,...
    'fontsize',ftsza,...
    'box','off',...
    'tickdir','both',...
    'TickLength',[0.02 0.01])

subplot(2,2,3)
yyaxis left
plot(E0_field_height(3,pick_E0_90_km),log10(Q_0v_non_max_N2(ind(3,pick_E0_90_km))),'b','LineWidth',2) 
hold on;
plot(E0_field_height(3,pick_E0_90_km),log10(Q_0v_maxw_N2_L(ind(3,pick_E0_90_km))),'b-.','LineWidth',2)
hold off; 
ax = gca; 
ax.YColor = 'b'; 
ylabel('log Q [eV cm^3 s^{-1}]','FontSize',ftszl)
xlabel('E_0 [V/m]','FontSize',ftszl)
ylim([-15 -9])
yticks([-15 -14 -13 -12 -11 -10 -9])
title('90 km','FontSize',ftszt) 
yyaxis right
plot(E0_field_height(3,pick_E0_90_km),Te_NM_L(ind(3,pick_E0_90_km)),'r--','LineWidth',2) 
legend()
ylabel('Te [K]')
ylim([0 4000])
xlim([0 4])
legend('Q_{NM}','Q_{M}','T_L',...
       'Location','northwest',...
       'FontSize',ftszlg)
ax = gca; 
ax.YColor = 'r';
set(gca,...
    'fontsize',ftsza,...
    'box','off',...
    'tickdir','both',...
    'TickLength',[0.02 0.01])

subplot(2,2,4)
yyaxis left
plot(E0_field_height(4,pick_E0_100_km),log10(Q_0v_non_max_N2(ind(4,pick_E0_100_km))),'b','LineWidth',2) 
hold on;
plot(E0_field_height(4,pick_E0_100_km),log10(Q_0v_maxw_N2_L(ind(4,pick_E0_100_km))),'b-.','LineWidth',2)
hold off; 
ax = gca; 
ax.YColor = 'b'; 
ylabel('log Q [eV cm^3 s^{-1}]','FontSize',ftszl)
xlabel('E_0 [V/m]','FontSize',ftszl)
ylim([-15 -9])
yticks([-15 -14 -13 -12 -11 -10 -9])
title('100 km','FontSize',ftszt) 
yyaxis right
plot(E0_field_height(4,pick_E0_100_km),Te_NM_L(ind(4,pick_E0_100_km)),'r--','LineWidth',2) 
legend()
ylabel('Te [K]')
ylim([0 4000])
xlim([0 4])
legend('Q_{NM}','Q_{M}','T_L',...
       'Location','northwest',...
       'FontSize',ftszlg)
ax = gca; 
ax.YColor = 'r';
set(gca,...
    'fontsize',ftsza,...
    'box','off',...
    'tickdir','both',...
    'TickLength',[0.02 0.01]) 

sgtitle('Vibrational excitation of molecular nitrogen')

%% O2: M->L

figure

subplot(2,2,1)
yyaxis left
plot(E0_field_height(1,pick_E0_70_km),log10(Q_0v_non_max_O2(ind(1,pick_E0_70_km))),'b','LineWidth',2) 
hold on;
plot(E0_field_height(1,pick_E0_70_km),log10(Q_0v_maxw_O2_L(ind(1,pick_E0_70_km))),'b-.','LineWidth',2)
hold off; 
ax = gca; 
ax.YColor = 'b'; 
ylabel('log Q [eV cm^3 s^{-1}]','FontSize',ftszl)
ylim([-13 -8])
yticks([-13 -12 -11 -10 -9 -8])
title('70 km','FontSize',ftszt) 
yyaxis right
plot(E0_field_height(1,pick_E0_70_km),Te_NM_L(ind(1,pick_E0_70_km)),'r--','LineWidth',2) 
legend()
ylabel('Te [K]')
ylim([0 4000])
xlim([0 4])
legend('Q_{NM}','Q_{M}','T_L',...
       'Location','northwest',...
       'FontSize',ftszlg)
ax = gca; 
ax.YColor = 'r';
set(gca,...
    'fontsize',ftsza,...
    'box','off',...
    'tickdir','both',...
    'TickLength',[0.02 0.01])

subplot(2,2,2)
yyaxis left
plot(E0_field_height(2,pick_E0_80_km),...
     log10(Q_0v_non_max_O2(ind(2,pick_E0_80_km))),...
     'b','LineWidth',2) 
hold on;
plot(E0_field_height(2,pick_E0_80_km),...
     log10(Q_0v_maxw_O2_L(ind(2,pick_E0_80_km))),...
     'b-.','LineWidth',2)
hold off; 
ax = gca; 
ax.YColor = 'b'; 
ylabel('log Q [eV cm^3 s^{-1}]','FontSize',ftszl)
ylim([-13 -8])
yticks([-13 -12 -11 -10 -9 -8])
title('80 km','FontSize',ftszt) 
yyaxis right
plot(E0_field_height(2,pick_E0_80_km),...
     Te_NM_L(ind(2,pick_E0_80_km)),...
     'r--','LineWidth',2) 
legend()
ylabel('Te [K]')
ylim([0 4000])
xlim([0 4])
legend('Q_{NM}','Q_{M}','T_L',...
       'Location','northwest',...
       'FontSize',ftszlg)
ax = gca; 
ax.YColor = 'r';
set(gca,...
    'fontsize',ftsza,...
    'box','off',...
    'tickdir',...
    'both',...
    'TickLength',[0.02 0.01])

subplot(2,2,3)
yyaxis left
plot(E0_field_height(3,pick_E0_90_km),...
     log10(Q_0v_non_max_O2(ind(3,pick_E0_90_km))),...
     'b','LineWidth',2) 
hold on;
plot(E0_field_height(3,pick_E0_90_km),...
     log10(Q_0v_maxw_O2_L(ind(3,pick_E0_90_km))),...
     'b-.','LineWidth',2)
hold off; 
ax = gca; 
ax.YColor = 'b'; 
ylabel('log Q_ [eV cm^3 s^{-1}]','FontSize',ftszl)
xlabel('E_0 [V/m]','FontSize',ftszl)
ylim([-13 -8])
yticks([-13 -12 -11 -10 -9 -8])
title('90 km','FontSize',ftszt) 
yyaxis right
plot(E0_field_height(3,pick_E0_90_km),...
     Te_NM_L(ind(3,pick_E0_90_km)),...
     'r--','LineWidth',2) 
legend()
ylabel('Te [K]')
ylim([0 4000])
xlim([0 4])
legend('Q_{NM}','Q_{M}','T_L',...
       'Location','northwest',...
       'FontSize',ftszlg)
ax = gca; 
ax.YColor = 'r';
set(gca,...
    'fontsize',ftsza,...
    'box','off',...
    'tickdir','both',...
    'TickLength',[0.02 0.01])

subplot(2,2,4)
yyaxis left
plot(E0_field_height(4,pick_E0_100_km),log10(Q_0v_non_max_O2(ind(4,pick_E0_100_km))),'b','LineWidth',2) 
hold on;
plot(E0_field_height(4,pick_E0_100_km),log10(Q_0v_maxw_O2_L(ind(4,pick_E0_100_km))),'b-.','LineWidth',2)
hold off; 
ax = gca; 
ax.YColor = 'b'; 
ylabel('log Q [eV cm^3 s^{-1}]','FontSize',ftszl)
xlabel('E_0 [V/m]','FontSize',ftszl)
ylim([-13 -8])
yticks([-13 -12 -11 -10 -9 -8])
title('100 km','FontSize',ftszt) 
yyaxis right
plot(E0_field_height(4,pick_E0_100_km),Te_NM_L(ind(4,pick_E0_100_km)),'r--','LineWidth',2) 
legend()
ylabel('Te [K]')
ylim([0 4000])
xlim([0 4])
legend('Q_{NM}','Q_{M}','T_L',...
       'Location','northwest',...
       'FontSize',ftszlg)
ax = gca; 
ax.YColor = 'r';
set(gca,...
    'fontsize',ftsza,...
    'box','off',...
    'tickdir','both',...
    'TickLength',[0.02 0.01]) 

sgtitle('Vibrational excitation of molecular oxygen')

%% Fs O: M->L

figure 

subplot(2,2,1)
yyaxis left
plot(E0_field_height(1,pick_E0_70_km),...
     log10(Q_non_max_O_fs(ind(1,pick_E0_70_km))),...
     'b','LineWidth',2) 
hold on;
plot(E0_field_height(1,pick_E0_70_km),...
     log10(Q_maxw_O_fs_L(ind(1,pick_E0_70_km))),...
     'b-.','LineWidth',2.0)
hold off; 
ax = gca; 
ax.YColor = 'b'; 
ylabel('log Q [eV cm^3 s^{-1}]','FontSize',ftszl)
ylim([-11.6 -10.8])
title('70 km','FontSize',ftszt) 
yyaxis right
plot(E0_field_height(1,pick_E0_70_km),...
     Te_NM_L(ind(1,pick_E0_70_km)),...
     'r--','LineWidth',2) 
legend()
ylabel('Te [K]')
ylim([0 4000])
xlim([0 4])
legend('Q_{NM}','Q_{M}','T_L',...
       'Location','northwest',...
       'FontSize',ftszlg)
ax = gca; 
ax.YColor = 'r';
set(gca,...
    'fontsize',ftsza,...
    'box','off',...
    'tickdir','both',...
    'TickLength',[0.02 0.01]) 

subplot(2,2,2)
yyaxis left
plot(E0_field_height(2,pick_E0_80_km),...
     log10(Q_non_max_O_fs(ind(2,pick_E0_80_km))),...
     'b','LineWidth',2) 
hold on;
plot(E0_field_height(2,pick_E0_80_km),...
     log10(Q_maxw_O_fs_L(ind(2,pick_E0_80_km))),...
     'b-.','LineWidth',2.0)
hold off; 
ax = gca; 
ax.YColor = 'b'; 
ylabel('log Q [eV cm^3 s^{-1}]','FontSize',ftszl)
ylim([-11.6 -10.8])
title('80 km','FontSize',ftszt) 
yyaxis right
plot(E0_field_height(2,pick_E0_80_km),...
     Te_NM_L(ind(2,pick_E0_80_km)),...
     'r--','LineWidth',2) 
legend()
ylabel('Te [K]')
ylim([0 4000])
xlim([0 4])
legend('Q_{NM}','Q_{M}','T_L',...
       'Location','northwest',...
       'FontSize',ftszlg)
ax = gca; 
ax.YColor = 'r';
set(gca,...
    'fontsize',ftsza,...
    'box','off',...
    'tickdir','both',...
    'TickLength',[0.02 0.01]) 

subplot(2,2,3)
yyaxis left
plot(E0_field_height(3,pick_E0_90_km),...
     log10(Q_non_max_O_fs(ind(3,pick_E0_90_km))),...
     'b','LineWidth',2) 
hold on;
plot(E0_field_height(3,pick_E0_90_km),...
     log10(Q_maxw_O_fs_L(ind(3,pick_E0_90_km))),...
     'b-.','LineWidth',2.0)
hold off; 
ax = gca; 
ax.YColor = 'b'; 
ylabel('log Q [eV cm^3 s^{-1}]','FontSize',ftszl)
xlabel('E_0 [V/m]','FontSize',ftszl)
ylim([-11.6 -10.8])
title('90 km','FontSize',ftszt) 
yyaxis right
plot(E0_field_height(3,pick_E0_90_km),...
     Te_NM_L(ind(3,pick_E0_90_km)),...
     'r--','LineWidth',2) 
legend()
ylabel('Te [K]')
ylim([0 4000])
xlim([0 4])
legend('Q_{NM}','Q_{M}','T_L',...
       'Location','northwest',...
       'FontSize',ftszlg)
ax = gca; 
ax.YColor = 'r';
set(gca,...
    'fontsize',ftsza,...
    'box','off',...
    'tickdir','both',...
    'TickLength',[0.02 0.01]) 

subplot(2,2,4)
yyaxis left
plot(E0_field_height(4,pick_E0_100_km),...
     log10(Q_non_max_O_fs(ind(4,pick_E0_100_km))),...
     'b','LineWidth',2) 
hold on;
plot(E0_field_height(4,pick_E0_100_km),...
     log10(Q_maxw_O_fs_L(ind(4,pick_E0_100_km))),...
     'b-.','LineWidth',2.0)
hold off; 
ax = gca; 
ax.YColor = 'b'; 
ylabel('log Q [eV cm^3 s^{-1}]','FontSize',ftszl)
xlabel('E_0 [V/m]','FontSize',ftszl)
ylim([-11.6 -10.8])
title('100 km','FontSize',ftszt) 
yyaxis right
plot(E0_field_height(4,pick_E0_100_km),...
     Te_NM_L(ind(4,pick_E0_100_km)),...
     'r--','LineWidth',2) 
legend()
ylabel('Te [K]')
ylim([0 4000])
xlim([0 4])
legend('Q_{NM}','Q_{M}','T_L',...
       'Location','northwest',...
       'FontSize',ftszlg)
ax = gca; 
ax.YColor = 'r';
set(gca,...
    'fontsize',ftsza,...
    'box','off',...
    'tickdir','both',...
    'TickLength',[0.02 0.01]) 

sgtitle('Excitation of fine structures in atomic oxygen')


%% Create data for bar plot 

ch_z0 = [3 3]; % 100 km 
ch_E0 = [3 12]; % E_0 2.5 and 2.5  V/m 

data_NM = zeros(length(ch_z0),10); 
data_M_2nd = zeros(length(ch_z0),10); 
data_M_L = zeros(length(ch_z0),10); 
for i = 1:length(ch_z0)

    a1 = ch_z0(i);
    b1 = ch_E0(i); 

    % Non-Maxwellian
    O_fs_NM = L_Ofine_1_0_dIe(ind(a1,b1)) + ...
              L_Ofine_2_0_dIe(ind(a1,b1)) + ...
              L_Ofine_2_1_dIe(ind(a1,b1)); 
    data_NM(i,:) = [L_vib_N2_nonmax(ind(a1,b1)) L_vib_O2_nonmax(ind(a1,b1))  L_N2_rot(ind(a1,b1)) ...
        L_O2_rot(ind(a1,b1)) L_el(ind(a1,b1)) O_fs_NM  L_O2a1Dg_dIe(ind(a1,b1))  ...
        L_O2b1Sgp_dIe(ind(a1,b1)) L_O3P1S_dIe(ind(a1,b1)) L_O3P1D_dIe(ind(a1,b1))]; 
  
    % Maxwellian at T_2nd
    O_fs_M_2nd = L_oxygen_10_pavlov_tables_2nd(ind(a1,b1))+L_oxygen_20_pavlov_tables_2nd(ind(a1,b1))+...
        L_oxygen_21_pavlov_tables_2nd(ind(a1,b1)); 
    data_M_2nd(i,:) = [L_vib_N2_maxw_2nd(ind(a1,b1)) L_vib_O2_maxw_2nd(ind(a1,b1))  ...
        L_rot_N2_pavlov_tables_2nd(ind(a1,b1)) L_rot_O2_pavlov_tables_2nd(ind(a1,b1)) ...
        L_elas_pavlov_tables_2nd(ind(a1,b1)) O_fs_M_2nd L_a1delta_O2_pavlov_tables_2nd(ind(a1,b1))  ...
        L_b1sigma_O2_pavlov_tables_2nd(ind(a1,b1)) 0 0]; 

    % Maxwellian at T_L
    O_fs_M_L = L_oxygen_10_pavlov_tables_L(ind(a1,b1))+L_oxygen_20_pavlov_tables_L(ind(a1,b1))+...
        L_oxygen_21_pavlov_tables_L(ind(a1,b1));
    data_M_L(i,:) = [L_vib_N2_pavlov_tables_L(ind(a1,b1)) L_vib_O2_pavlov_tables_L(ind(a1,b1))  ...
        L_rot_N2_pavlov_tables_L(ind(a1,b1)) L_rot_O2_pavlov_tables_L(ind(a1,b1)) ...
        L_elas_pavlov_tables_L(ind(a1,b1)) O_fs_M_L L_a1delta_O2_pavlov_tables_L(ind(a1,b1))  ...
        L_b1sigma_O2_pavlov_tables_L(ind(a1,b1)) 0 0]; 

end 

legend_text = {'N_2 vib.','O_2 vib.',...
               'N_2 rot.','O_2 rot.',...
               'ela.','O fs','O2a1','O2b1'}; 
X = categorical({'M:T_e=2nd','Non-M.','M:T_e=L'}); 
X = reordercats(X,{'M:T_e=2nd','Non-M.','M:T_e=L'});


%% Scatter and bar plot

% Graphical setup
ftsza_2  = 11.5; 
ftszl_2  = 14.5; 
ftszlg_2 =  9.5;
ftszt_2  = 16; 

% Total Maxwellian cooling rates at Te = T_2nd
L_M_tot = (L_vib_N2_maxw_2nd + ...
           L_vib_O2_maxw_2nd + ...
           L_oxygen_10_pavlov_tables_2nd +...
           L_oxygen_20_pavlov_tables_2nd + ...
           L_oxygen_21_pavlov_tables_2nd +...
           L_rot_N2_pavlov_tables_2nd + ...
           L_rot_O2_pavlov_tables_2nd + ...
           L_elas_pavlov_tables_2nd + ...
           L_a1delta_O2_pavlov_tables_2nd + ...
           L_b1sigma_O2_pavlov_tables_2nd ) * eVtoJ*(1/m_neg3_to_cm_neg3); 


% For plotting 
my_values = [8 9 10 11 12 13 14 15 ...
    16 17 18 19 20 21 22 23 24 25 26 27 28 ...
    31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 ...
    46 47 48 49 50 51 52 53 54 55 56 57 58 59 60]; 
ratio_L = zeros(1,length(my_values)); 
ratio_Te = zeros(1,length(my_values)); 
L_plot = zeros(1,length(my_values)); 
z_plot = zeros(1,length(my_values)); 
for h = 1:length(my_values)
    ratio_L(h) = L_M_tot(my_values(h))./L(my_values(h)); 
    ratio_Te(h) = Te_NM_2nd(my_values(h))./Te_NM_L(my_values(h)); 
    L_plot(h) = L(my_values(h)); 
    z_plot(h) = z(my_values(h)); 
end 

figure
subplot(2,2,1)
scatter(L_plot,z_plot',37,ratio_L,'filled')
set(gca,'xscale','log')
set(gca,...
    'fontsize',ftsza_2,...
    'box','on',...
    'tickdir','both',...
    'TickLength',[0.015 0.008]) 
xlabel('Q_{abs} [J m^{-3} s^{-1}]','FontSize',ftszl_2)
ylabel('Height [km]','FontSize',ftszl_2)
colorbar
colormap(parula)
xlim([1e-11 1e-5])
xticks([1e-11 1e-9 1e-7 1e-5])
title('(a)','FontSize',ftszt_2)
ylim([60 110])
yticks([60 70 80 90 100 110])

subplot(2,2,2)
scatter(L_plot,z_plot',37,ratio_Te,'filled')
set(gca,'xscale','log')
set(gca,...
    'fontsize',ftsza_2,...
    'box','on',...
    'tickdir','both',....
    'TickLength',[0.015 0.008]) 
xlabel('Q_{abs} [J m^{-3} s^{-1}]','FontSize',ftszl_2)
colorbar
colormap("turbo")
xlim([1e-11 1e-5])
xticks([1e-11 1e-9 1e-7 1e-5])
title('(b)','FontSize',ftszt_2)
ylim([60 110])
yticks([60 70 80 90 100 110])

subplot(2,2,3)
bar(X, ...
    [data_M_2nd(1,:)*eVtoJ*(1/m_neg3_to_cm_neg3); ...
     data_NM(1,:)*eVtoJ*(1/m_neg3_to_cm_neg3); ...
     data_M_L(1,:)*eVtoJ*(1/m_neg3_to_cm_neg3)],...
    "stacked")
set(gca,...
    'fontsize',ftsza_2,...
    'box','off',...
    'tickdir','both',...
    'TickLength',[0.02 0.01]) 
title('(c)','FontSize',ftszt_2)
ylabel('L_{tot} [J m^{-3} s^{-1}]','FontSize',ftszl_2)
subplot(2,2,4) 
bar(X,...
    [data_M_2nd(2,:)*eVtoJ*(1/m_neg3_to_cm_neg3); ...
     data_NM(2,:)*eVtoJ*(1/m_neg3_to_cm_neg3); ...
     data_M_L(2,:)*eVtoJ*(1/m_neg3_to_cm_neg3)],...
    "stacked")
set(gca,...
    'fontsize',ftsza_2,...
    'box','off',...
    'tickdir','both',...
    'TickLength',[0.02 0.01]) 
title('(d)','FontSize',ftszt_2)
legend(legend_text,...
       'Location','northeast',...
       'FontSize',ftszlg_2,...
       'NumColumns',2)
colororder('glow12')

%% Plot ratio non-Mawellian vs Maxwellian 

% Graphical setup
ftszl_3 = 16.0; 
ftszlg_3 = 10.1;
ftszt_3 = 17; 

pick_E0_70_km_2 = 7:size(E0_field_height,2); 
pick_E0_80_km_2 = 5:size(E0_field_height,2); 
pick_E0_90_km_2 = 5:size(E0_field_height,2); 
pick_E0_100_km_2 = 5:size(E0_field_height,2); 

col1 = [0,204,204]/255; 
col2 = [0,0,204]/255; 
col3 = [127,0,255]/255; 
col4 = [204,0,102]/255; 

figure 

subplot(2,2,1)
semilogx(L(ind(1,pick_E0_70_km_2)),Q_0v_non_max_N2(ind(1,pick_E0_70_km_2))./....
    Q_0v_maxw_N2_L(ind(1,pick_E0_70_km_2)),'--','LineWidth',2,'color',col1)
hold on;
semilogx(L(ind(2,pick_E0_80_km_2)),Q_0v_non_max_N2(ind(2,pick_E0_80_km_2))./...
    Q_0v_maxw_N2_L(ind(2,pick_E0_80_km_2)),':','LineWidth',2,'color',col2)
semilogx(L(ind(3,pick_E0_90_km_2)),Q_0v_non_max_N2(ind(3,pick_E0_90_km_2))./...
    Q_0v_maxw_N2_L(ind(3,pick_E0_90_km_2)),'-.','LineWidth',2,'color',col3)
semilogx(L(ind(4,pick_E0_100_km_2)),Q_0v_non_max_N2(ind(4,pick_E0_100_km_2))./...
    Q_0v_maxw_N2_L(ind(4,pick_E0_100_km_2)),'LineWidth',2,'color',col4)
hold off; 
set(gca,...
    'fontsize',ftsza,...
    'box','on',...
    'tickdir','both',...
    'TickLength',[0.017 0.009]) 
ylabel('Ratio','FontSize',ftszl_3)
xlim([1e-11 1e-5])
xticks([1e-11 1e-10 1e-9 1e-8 1e-7 1e-6 1e-5])
title('N2 vib.','FontSize',ftszt_3) 
legend('70 km','80 km','90 km','100 km',...
       'Location','northwest',...
       'FontSize',ftszlg_3)

subplot(2,2,2)
semilogx(L(ind(1,pick_E0_70_km_2)),Q_0v_non_max_O2(ind(1,pick_E0_70_km_2))./.....
    Q_0v_maxw_O2_L(ind(1,pick_E0_70_km_2)),'--','LineWidth',2,'color',col1)
hold on;
semilogx(L(ind(2,pick_E0_80_km_2)),Q_0v_non_max_O2(ind(2,pick_E0_80_km_2))./...
    Q_0v_maxw_O2_L(ind(2,pick_E0_80_km_2)),':','LineWidth',2,'color',col2)
semilogx(L(ind(3,pick_E0_90_km_2)),Q_0v_non_max_O2(ind(3,pick_E0_90_km_2))./...
    Q_0v_maxw_O2_L(ind(3,pick_E0_90_km_2)),'-.','LineWidth',2,'color',col3)
semilogx(L(ind(4,pick_E0_100_km_2)),Q_0v_non_max_O2(ind(4,pick_E0_100_km_2))./...
    Q_0v_maxw_O2_L(ind(4,pick_E0_100_km_2)),'LineWidth',2,'color',col4)
hold off; 
set(gca,...
    'fontsize',ftsza,...
    'box','on',...
    'tickdir','both',...
    'TickLength',[0.017 0.009]) 
xlim([1e-11 1e-5])
xticks([1e-11 1e-10 1e-9 1e-8 1e-7 1e-6 1e-5])
title('O2 vib.','FontSize',ftszt_3) 
legend('70 km','80 km','90 km','100 km',...
       'Location','northwest',...
       'FontSize',ftszlg_3)

subplot(2,2,3)
semilogx(L(ind(1,:)),Q_non_max_rot_N2(ind(1,:))./Q_maxw_rot_N2_L(ind(1,:)),'--','LineWidth',2,'color',col1)
hold on;
semilogx(L(ind(2,:)),Q_non_max_rot_N2(ind(2,:))./Q_maxw_rot_N2_L(ind(2,:)),':','LineWidth',2,'color',col2)
semilogx(L(ind(3,:)),Q_non_max_rot_N2(ind(3,:))./Q_maxw_rot_N2_L(ind(3,:)),'-.','LineWidth',2,'color',col3)
semilogx(L(ind(4,:)),Q_non_max_rot_N2(ind(4,:))./Q_maxw_rot_N2_L(ind(4,:)),'LineWidth',2,'color',col4)
hold off; 
set(gca,...
    'fontsize',ftsza,...
    'box','on',...
    'tickdir','both',...
    'TickLength',[0.017 0.009]) 
xlabel('Q_{abs} [J m^{-3} s^{-1}]','FontSize',ftszl_2)
ylabel('Ratio','FontSize',ftszl_3)
xlim([1e-11 1e-5])
xticks([1e-11 1e-10 1e-9 1e-8 1e-7 1e-6 1e-5])
title('N_2 rot.','FontSize',ftszt_3) 
legend('70 km','80 km','90 km','100 km',...
       'Location','northwest',...
       'FontSize',ftszlg_3)

subplot(2,2,4)
semilogx(L(ind(1,:)),Q_non_max_O_fs(ind(1,:))./Q_maxw_O_fs_L(ind(1,:)),'--','LineWidth',2,'color',col1)
hold on;
semilogx(L(ind(2,:)),Q_non_max_O_fs(ind(2,:))./Q_maxw_O_fs_L(ind(2,:)),':','LineWidth',2,'color',col2)
semilogx(L(ind(3,:)),Q_non_max_O_fs(ind(3,:))./Q_maxw_O_fs_L(ind(3,:)),'-.','LineWidth',2,'color',col3)
semilogx(L(ind(4,:)),Q_non_max_O_fs(ind(4,:))./Q_maxw_O_fs_L(ind(4,:)),'LineWidth',2,'color',col4)
hold off; 
set(gca,...
    'fontsize',ftsza,...
    'box','on',...
    'tickdir','both',...
    'TickLength',[0.017 0.009]) 
xlabel('Q_{abs} [J m^{-3} s^{-1}]','FontSize',ftszl_2)
xlim([1e-11 1e-5])
xticks([1e-11 1e-10 1e-9 1e-8 1e-7 1e-6 1e-5])
title('O fs','FontSize',ftszt_3) 
legend('70 km','80 km','90 km','100 km',...
       'Location','northwest',...
       'FontSize',ftszlg_3)











