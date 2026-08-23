%% Plot figures for non-Maxwellian distributions
% Plot figures for the following paper: non-Maxwellian electron distribution in the D
% region during artificial heating (Paper I): Model development and
% electron temperature

%% Graphics setup
font_a  = 12.0; % Font for tics
font_lb = 14; % Font size for label
font_lg = 11.5; % Font size for legends
font_t  = 15; % Font size for title

%% Some constant 
me                = 9.10938291e-31;  % Electron mass [kg]
eVtoJ             = 1.602176634e-19; % Convert from [eV] to [J] 
k_b               = 1.380649e-23;    % Boltzmann constant [J/K] 
k_b_eV            = 8.617333e-5;     % Boltzmann constant [eV/K] 
q_e               = 1.602176634e-19; % Elementary charge of 1 electron [C]
cm_to_m           = 0.01;            % Convert from [cm] to [m]
g_to_kg           = 1/1000;          % Convert from [g] to [kg]
m_neg3_to_cm_neg3 = 1e-6;            % Convert densities from [m^-3] to [cm^-3]

%% Geophysical and ionospheric
f_HF = 4.6e6;                 % Radio wave frequency [Hz] 
w_HF = 2*pi*f_HF;             % Angular frequency of radio wave [rad/s] 
c    = 2.99792458e8;          % Speed of light [m/s]
n_e  = 1e10;                  % Electron density [m^-3]
Tn   = 200;                   % Neutral temperature [K]

%% Load and process data for cross sections
xs_data = get_data_less_xs(); 

%% Iteration for 70, 100 and 120 km 
z0_1         = [70e3 85e3 100e3];    % Altitude [m]
E0_field_1   = [0.8 1.5 2.5];   % Electric field [V/m]
f_E_stubbe_1 = cell(length(z0_1),length(E0_field_1)); 
f0_stubbe_1  = cell(length(z0_1),length(E0_field_1)); 
A_stubbe_1   = zeros(length(z0_1),length(E0_field_1)); 
f_max_E_1    = cell(length(z0_1),length(E0_field_1)); 
Te_1         = zeros(length(z0_1),length(E0_field_1)); 

wbh = waitbar(0,'Working, working, step by step');
for j = 1:length(z0_1)
  for i = 1:length(E0_field_1)
    try
      waitbar((i+length(E0_field_1)*(j-1))/(length(z0_1)*length(E0_field_1)));
    catch
      % whatever
    end

    % Call MATLAB MSIS
    [rho,~] = neutral_output_MSIS(z0_1(j)); 
    
    % Neutral densities MATLAB MSIS
    nO = rho(:,2)';   % [1/m^3] 
    nN2 = rho(:,3)';  % [1/m^3] 
    nO2 = rho(:,4)';  % [1/m^3]
    
    % Max electron energy [eV]
    max_E = 6.0; 
    
    % Compute non-Maxwellian distribution 
    [f0_exp_stubbe,v_bins_stubbe,Te_nm_2nd,~,A,~,~,~] = ...
    compute_f0_stubbe_less_xs(E0_field_1(i),n_e,Tn,nO,nN2,nO2,f_HF,xs_data,max_E);
    
    % Energy/speed in the middle fo the bin
    v_bin_mid = v_bins_stubbe(1:end-1)/2 + v_bins_stubbe(2:end)/2; % [m/s]
    E_bin_mid_eV = (me*v_bin_mid.^2)/(2*eVtoJ); % [eV]
    
    % Speed distribution [#e/m^3 * s/m]
    f0_speed_distri = n_e*A*4*pi*v_bin_mid.^2.*f0_exp_stubbe;
    
    % Convert speed distribution to energy distribution [#e/m^3 * 1/eV]
    f_E_stubbe_1{j,i} = convert_fv2fE(v_bins_stubbe,f0_speed_distri);
    
    % Distribution witn no unit
    f0_stubbe_1{j,i} = f0_exp_stubbe; 
    
    % Normalization constant
    A_stubbe_1(j,i) = A; 
    
    % Electron temperature [K] (second moment of distribution function)
    Te_1(j,i) = Te_nm_2nd;
    
    % Maxwellian speed distribution [#e/m^3 * s/m]
    A_max = (me/(2*pi*k_b*Te_1(j,i)))^(3/2);
    f_max_v = n_e*A_max*4*pi*v_bin_mid.^2.*exp(-me*v_bin_mid.^2/(2*k_b*Te_1(j,i))); 
    
    % Maxwellian energy distribution[#e/m^3 * 1/eV] 
    f_max_E_1{j,i} = convert_fv2fE(v_bins_stubbe,f_max_v);  
    
  end 
end 
close(wbh)

%% Plot figure 2

% Some colors
col_gray = [170, 170, 170]/255; 
col1     = [ 76, 153,   0]/255; 
col3     = [ 82, 213,   0]/255; 
col4     = [  0, 145, 147]/255; 
col6     = [  0, 213, 214]/255; 
col7     = [ 20, 120, 199]/255; 
col9     = [ 45, 140, 243]/255; 

% Limit axis
lim_x_min = 1e-2; 
lim_x_max = 5;
lim_y_min = 1e2; 
lim_y_max = 1e12; 

figure

subplot(1,3,1)
loglog(E_bin_mid_eV,f_max_E_1{1,1},'LineWidth',2.0,'color',col_gray)
hold on;
loglog(E_bin_mid_eV,f_max_E_1{1,3},'-.','LineWidth',2.0,'color',col_gray)
loglog(E_bin_mid_eV,f_E_stubbe_1{1,1},'LineWidth',2.0,'color',col1)
loglog(E_bin_mid_eV,f_E_stubbe_1{1,3},'-.','LineWidth',2.0,'color',col3)
hold off; 
set(gca,...
    'fontsize',font_a,...
    'box','on',...
    'tickdir','both',...
    'TickLength',[0.009 0.008])
legend(sprintf('%0.0f %s',Te_1(1,1),'K'),...
       sprintf('%0.0f %s',Te_1(1,3),'K'),...
       '0.8 V/m',...
       '2.5 V/m',...
       'Location','southwest',...
       'FontSize',font_lg)
title('70 km','FontSize',font_t)
xlabel('Energy [eV]','FontSize',font_lb)
ylabel('f_0(E) [m^{-3} eV^{ -1}]','FontSize',font_lb)
xlim([lim_x_min lim_x_max])
xticks([1e-2 1e-1 1e0 5])
ylim([lim_y_min lim_y_max])

subplot(1,3,2)
loglog(E_bin_mid_eV,f_max_E_1{2,1},'LineWidth',2.0,'color',col_gray)
hold on;
loglog(E_bin_mid_eV,f_max_E_1{2,3},'-.','LineWidth',2.0,'color',col_gray)
loglog(E_bin_mid_eV,f_E_stubbe_1{2,1},'LineWidth',2.0,'color',col4)
loglog(E_bin_mid_eV,f_E_stubbe_1{2,3},'-.','LineWidth',2.0,'color',col6)
hold off; 
set(gca,...
    'fontsize',font_a,...
    'box','on',...
    'tickdir','both',...
    'TickLength',[0.009 0.008])
legend(sprintf('%0.0f %s',Te_1(2,1),'K'),...
       sprintf('%0.0f %s',Te_1(2,3),'K'),...
       '0.8 V/m',...
       '2.5 V/m',...
       'Location','southwest',...
       'FontSize',font_lg)
title('85 km','FontSize',font_t)
xlabel('Energy [eV]','FontSize',font_lb)
xlim([lim_x_min lim_x_max])
xticks([1e-2 1e-1 1e0 5])
ylim([lim_y_min lim_y_max])

subplot(1,3,3)
loglog(E_bin_mid_eV,f_max_E_1{3,1},'LineWidth',2.0,'color',col_gray)
hold on;
loglog(E_bin_mid_eV,f_max_E_1{3,3},'-.','LineWidth',2.0,'color',col_gray)
loglog(E_bin_mid_eV,f_E_stubbe_1{3,1},'LineWidth',2.0,'color',col7)
loglog(E_bin_mid_eV,f_E_stubbe_1{3,3},'-.','LineWidth',2.0,'color',col9)
hold off; 
set(gca,...
    'fontsize',font_a,...
    'box','on',...
    'tickdir','both',...
    'TickLength',[0.009 0.008])
legend(sprintf('%0.0f %s',Te_1(3,1),'K'),...
       sprintf('%0.0f %s',Te_1(3,3),'K'),...
       '0.8 V/m',...
       '2.5 V/m',...
       'Location','southwest',...
       'FontSize',font_lg)
title('100 km','FontSize',font_t)
xlabel('Energy [eV]','FontSize',font_lb)
xlim([1e-2 lim_x_max])
xticks([lim_x_min 1e-1 1e0 5])
ylim([lim_y_min lim_y_max])

%% Distribution at 90 km and E0 = 2 V/m

z0_2 = 90e3; 
E0_field_2 = 2.0;

% Call MATLAB MSIS
[rho,temperature] = neutral_output_MSIS(z0_2); 

% Neutral densities MATLAB MSIS
nO = rho(:,2)';   % [1/m^3] 
nN2 = rho(:,3)';    % [1/m^3] 
nO2 = rho(:,4)';    % [1/m^3]

% Max electron energy [eV]
max_E = 6.0; 

% Compute non-Maxwellian distribution 
[f0_exp_stubbe,v_bins_stubbe,Te_nm_2nd,~,A,~,~,~] = ...
    compute_f0_stubbe_less_xs(E0_field_2,n_e,Tn,nO,nN2,nO2,f_HF,xs_data,max_E);

% Energy/speed in the middle fo the bin
v_bin_mid = v_bins_stubbe(1:end-1)/2 + v_bins_stubbe(2:end)/2; % [m/s]
E_bin_mid = (me*v_bin_mid.^2)/2; % [J]
E_bin_mid_eV = (me*v_bin_mid.^2)/(2*eVtoJ); % [eV]

% Distribution [s^3/m^6]
f0 = A*f0_exp_stubbe;

% Electron temperature (second moment of distribution function)
Te_2 = Te_nm_2nd;

% Maxwellian distribution [s^3/m^6]
A_max = (me/(2*pi*k_b*Te_2))^(3/2);
f_max = A_max*exp(-me*v_bin_mid.^2/(2*k_b*Te_2));  

%% Call funnction for interpolation of cross section data 
xs_e_N2elastic = sigma_N2_el(E_bin_mid,xs_data{1}(:,1),xs_data{1}(:,2));
xs_e_O2elastic = sigma_O2_el(E_bin_mid,xs_data{2}(:,1),xs_data{2}(:,2));
xs_e_Oelastic  = sigma_O_el_2(E_bin_mid,xs_data{3}(:,1),xs_data{3}(:,2));
xs_e_N2vib0_1  = sigma_N2_vib_01(E_bin_mid,xs_data{4}(:,1),xs_data{4}(:,2)); 
xs_e_N2vib0_2  = sigma_N2_vib_02(E_bin_mid,xs_data{5}(:,1),xs_data{5}(:,2));
xs_e_N2vib0_3  = sigma_N2_vib_03(E_bin_mid,xs_data{6}(:,1),xs_data{6}(:,2)); 
xs_e_N2vib0_4  = sigma_N2_vib_04(E_bin_mid,xs_data{7}(:,1),xs_data{7}(:,2));
xs_e_N2vib0_5  = sigma_N2_vib_05(E_bin_mid,xs_data{8}(:,1),xs_data{8}(:,2)); 
xs_e_N2vib0_6  = sigma_N2_vib_06(E_bin_mid,xs_data{9}(:,1),xs_data{9}(:,2));
xs_e_N2vib0_7  = sigma_N2_vib_07(E_bin_mid,xs_data{10}(:,1),xs_data{10}(:,2)); 
xs_e_N2vib0_8  = sigma_N2_vib_08(E_bin_mid,xs_data{11}(:,1),xs_data{11}(:,2)); 
xs_e_N2vib0_9  = sigma_N2_vib_09(E_bin_mid,xs_data{12}(:,1),xs_data{12}(:,2)); 
xs_e_N2vib0_10 = sigma_N2_vib_010(E_bin_mid,xs_data{13}(:,1),xs_data{13}(:,2)); 
xs_e_N2rot0_2  = sigma_N2_rot_02(E_bin_mid,xs_data{15}(:,1),xs_data{15}(:,2)); 
xs_e_N2rot0_4  = sigma_N2_rot_04(E_bin_mid,xs_data{16}(:,1),xs_data{16}(:,2));
xs_e_N2rot0_6  = sigma_N2_rot_06(E_bin_mid,xs_data{17}(:,1),xs_data{17}(:,2)); 
xs_e_N2rot0_8  = sigma_N2_rot_08(E_bin_mid,xs_data{18}(:,1),xs_data{18}(:,2));
xs_e_O2a1Dg    = sigma_O2_a1_delta_g(E_bin_mid,xs_data{19}(:,1),xs_data{19}(:,2)); 
xs_e_O2b1Sgp   = sigma_O2_b1_delta_g(E_bin_mid,xs_data{20}(:,1),xs_data{20}(:,2)); 
xs_e_Ofine_1_0 = sigma_O_fine_str_10(E_bin_mid,xs_data{21}(:,1),xs_data{21}(:,2));  
xs_e_Ofine_2_0 = sigma_O_fine_str_20(E_bin_mid,xs_data{22}(:,1),xs_data{22}(:,2)); 
xs_e_Ofine_2_1 = sigma_O_fine_str_21(E_bin_mid,xs_data{23}(:,1),xs_data{23}(:,2));  
xs_e_O3P1D     = sigma_O_3P_1D(E_bin_mid,xs_data{24}(:,1),xs_data{24}(:,2)); 
xs_e_O3P1S     = sigma_O_3P_1S(E_bin_mid,xs_data{25}(:,1),xs_data{25}(:,2)); 

% Cross sections from Itikawa 1989
[E_old,xs_old] = p_data_O2_vib(); 

%% Plot figure 3 for 90 km and E0=2.0 V/m with cross sections 

% Some color
col60 = [0,0,120]/255; 

figure

subplot(2,1,1)
semilogy(E_bin_mid_eV,f0,'b','LineWidth',2.0)
hold on; 
semilogy(E_bin_mid_eV,f_max,'m-.','LineWidth',2.0)
hold off; 
xlim([0.005 2.6])
xticks([0.1 0.5 1 1.5 2])
ylim([1e-28 1e-16])
yticks([1e-28 1e-24 1e-20 1e-16])
set(gca,...
    'fontsize',font_a,...
    'box','on',...
    'tickdir','both',...
    'TickLength',[0.0035 0.0015])
ylabel('f_0 [s^3 m^{- 6}]','FontSize',font_lb)
title('(a)','FontSize',font_t)
xline(0.1)
xline(1)
xline(1.72)
text(0.03,1e-24,'1','FontSize',15)
text(0.5,1e-24,'2','FontSize',15)
text(1.3,1e-24,'3','FontSize',15)
text(2,1e-24,'4','FontSize',15)

subplot(2,1,2)
ph1 = semilogy(repmat(E_bin_mid_eV,1,8),...
               [xs_e_N2vib0_1 xs_e_N2vib0_2 xs_e_N2vib0_3 xs_e_N2vib0_4...
                xs_e_N2vib0_5 xs_e_N2vib0_6 xs_e_N2vib0_7 xs_e_N2vib0_8],...
                'LineWidth',2.5); 
cmp1 = colormap('winter'); 
cmlines(ph1,cmp1)
hold on;
semilogy(E_old*(1/eVtoJ),xs_old,'-.','LineWidth',2.5,'color',col60); 
ph3 = semilogy(repmat(E_bin_mid_eV,1,3),...
               [xs_e_Ofine_2_1 xs_e_Ofine_2_0 xs_e_Ofine_1_0], ...
               'LineWidth',2.5);
cmp3 = colormap('cool'); 
cmlines(ph3,cmp3)
hold on; 
ph4 = semilogy(repmat(E_bin_mid_eV,1,4), ...
               [xs_e_O2a1Dg xs_e_O2b1Sgp xs_e_O3P1D xs_e_O3P1S],...
               ':',...
               'LineWidth',2.5);
cmp4 = colormap('autumn'); 
cmlines(ph4,cmp4)
hold off; 
set(gca,...
    'fontsize',font_a,...
    'box','on',...
    'tickdir','both',...
    'TickLength',[0.005 0.002])
xlim([0.005 2.6])
xticks([0.1 0.5 1 1.5 2])
ylim([1e-24 1e-19])
yticks([1e-24 1e-23 1e-22 1e-21 1e-20 1e-19])
xlabel('Energy [eV]','FontSize',font_lb)
ylabel('Cross section [m^2]','FontSize',font_lb)
title('(b)','FontSize',font_t)
xline(0.1)
xline(1)
xline(1.72)


%% Load MC data

% Histogram [#e/m^3] for all the time steps
N_E_80_km  = load('N_E_80_km.mat'); 
N_E_80_km  = N_E_80_km.N_E; 
N_E_120_km = load('N_E_120_km.mat'); 
N_E_120_km = N_E_120_km.N_E; 
N_E_100_km = load('N_E_100_km.mat'); 
N_E_100_km = N_E_100_km.N_E; 
N_v_80_km  = load('N_v_80_km.mat'); 
N_v_80_km  = N_v_80_km.N_v; 
N_v_100_km = load('N_v_100_km.mat'); 
N_v_100_km = N_v_100_km.N_v; 
N_v_120_km = load('N_v_120_km.mat'); 
N_v_120_km = N_v_120_km.N_v; 
avg_t_80_km  = load('avg_t_80_km.mat'); 
avg_t_80_km  = avg_t_80_km.avg_t; 
avg_t_100_km = load('avg_t_100_km.mat'); 
avg_t_100_km = avg_t_100_km.avg_t; 
avg_t_120_km = load('avg_t_120_km.mat'); 
avg_t_120_km = avg_t_120_km.avg_t; 

% Index at equilibrium
ind_eq_80_km = 31; % at 0.003 s
ind_eq_100_km = 41; % at 0.3 s
ind_eq_120_km = 41; % 10.0 s

% Histogram [#e/m^3] as a function of v_bins_mid_MC: average over all time steps
N_v_avg_80_km = sum(N_v_80_km(:,ind_eq_80_km:end),2)/(length(avg_t_80_km)-ind_eq_80_km+1);
N_v_avg_100_km = sum(N_v_100_km(:,ind_eq_100_km:end),2)/(length(avg_t_100_km)-ind_eq_100_km+1);
N_v_avg_120_km = sum(N_v_120_km(:,ind_eq_120_km:end),2)/(length(avg_t_120_km)-ind_eq_120_km+1);

% Energy bins [eV]
E_bins_MC = 0:0.025:4;

% Bin-edges 
v_bins_MC = linspace(0,sqrt((2*eVtoJ*E_bins_MC(end))/me),length(E_bins_MC)); % [m/s]

% Convert speed bins to energy bins 
E_bins_conv_MC = (me*v_bins_MC.^2)/(2*eVtoJ); % eV

% Energies/speed in the middle of the bin
v_bin_mid_MC = v_bins_MC(1:end-1)/2 + v_bins_MC(2:end)/2; % [m/s]
E_bin_mid_MC = E_bins_MC(1:end-1)/2 + E_bins_MC(2:end)/2; 
E_bin_mid_conv_MC = ((me*v_bin_mid_MC.^2)/(2*eVtoJ))'; % [eV]

% Resolution of speed 
delta_v_MC = zeros(length(E_bins_conv_MC)-1,1); 
for i = 1:length(E_bins_conv_MC)-1
    delta_v_MC(i,1) = v_bins_MC(i+1)-v_bins_MC(i); % [eV] 
end 

% Resolution of converted energy
delta_E_conv_MC = zeros(length(E_bins_conv_MC)-1,1); 
for i = 1:length(E_bins_conv_MC)-1
    delta_E_conv_MC(i,1) = E_bins_conv_MC(i+1)-E_bins_conv_MC(i); % [eV] 
end 

% Convert to energy distribution [#e/m^3 1/eV]
f_E_MC_avg_80_km = (N_v_avg_80_km*10000)./delta_E_conv_MC;
f_E_MC_avg_100_km = (N_v_avg_100_km*10000)./delta_E_conv_MC;
f_E_MC_avg_120_km = (N_v_avg_120_km*10000)./delta_E_conv_MC;
f_E_MC_80_km = (N_v_80_km*10000)./delta_E_conv_MC;
f_E_MC_100_km = (N_v_100_km*10000)./delta_E_conv_MC;
f_E_MC_120_km = (N_v_120_km*10000)./delta_E_conv_MC;

% Convert to speed distribution [#e/m^3 s/m]
f_v_MC_avg_80_km = (N_v_avg_80_km*10000)./delta_v_MC;
f_v_MC_avg_100_km = (N_v_avg_100_km*10000)./delta_v_MC;
f_v_MC_avg_120_km = (N_v_avg_120_km*10000)./delta_v_MC;

%% Plot figure 4: MC energy distribution and energy
% as a function of time for 80, 100 and 120 km

aftsz = 12.0; % font-size for axe
font_lb = 13; % Font size for label
font_lg = 10.5; % Font size for legends
font_t = 14; % Font size for title

% Pick some nice indices to plot
pick_ind_legend_80_km = [1 2 8 31 length(avg_t_80_km)]; 
pick_ind_legend_100_km = [1 2 10 41 length(avg_t_100_km)]; 
pick_ind_legend_120_km = [1 2 4 29 length(avg_t_120_km)];

%fig = figure; 
figure
subplot(2,3,1)
phNE = semilogy(E_bin_mid_MC,N_E_80_km(:,pick_ind_legend_80_km));
set(phNE,'linewidth',2)
set(gca,'fontsize',aftsz)
xlabel('Energy [eV]','fontsize',font_lb)
ylabel('[#e]','fontsize',font_lb)
title('80 km','FontSize',font_t)
legend(strcat(num2str(1000*avg_t_80_km(pick_ind_legend_80_km)','%2.1f'),' ms'),'FontSize',font_lg)
xlim([0 2.7])
xticks([0 0.5 1 1.5 2 2.5])
subplot(2,3,2)
phNE = semilogy(E_bin_mid_MC,N_E_100_km(:,pick_ind_legend_100_km));
set(phNE,'linewidth',2)
set(gca,'fontsize',aftsz)
xlabel('Energy [eV]','fontsize',font_lb)
title('100 km','FontSize',font_t)
legend(strcat(num2str(avg_t_100_km(pick_ind_legend_100_km)','%2.2f'),' s'),'FontSize',font_lg)
set(gca,'fontsize',aftsz)
xlim([0 2.7])
xticks([0 0.5 1 1.5 2 2.5])
subplot(2,3,3)
phNE = semilogy(E_bin_mid_MC,N_E_120_km(:,pick_ind_legend_120_km));
set(phNE,'linewidth',2)
set(gca,'fontsize',aftsz)
xlabel('Energy [eV]','fontsize',font_lb)
title('120 km','FontSize',font_t)
legend(strcat(num2str(avg_t_120_km(pick_ind_legend_120_km)','%2.1f'),' s'),'FontSize',font_lg)
xlim([0 2.7])
xticks([0 0.5 1 1.5 2 2.5])
subplot(2,3,4)
imagesc(avg_t_80_km*1000,E_bin_mid_MC,log10(N_E_80_km))
set(gca,'fontsize',aftsz)
xlabel('Time [ms]','fontsize',font_lb)
ylabel('Energy [eV]','fontsize',font_lb)
ylim([0 2])
colormap('parula')
clim([0 6])

subplot(2,3,5)
imagesc(avg_t_100_km,E_bin_mid_MC,log10(N_E_100_km))
set(gca,'fontsize',aftsz)
xlabel('Time [s]','fontsize',font_lb)
ylim([0 2])
colormap('turbo')
clim([0 6])

subplot(2,3,6)
imagesc(avg_t_120_km,E_bin_mid_MC,log10(N_E_120_km))
set(gca,'fontsize',aftsz)
xlabel('Time [s]','fontsize',font_lb)
ylim([0 2])
colormap('hot')
label = '[#e]'; 
lg = 'log'; 
clim([0 6])
colorbar_labeled(label,lg);

% uncomment to save figure
% fig.PaperPositionMode = 'manual';
% orient(fig,'landscape')
% print(fig,'fig04.pdf','-dpdf')

%% Compute average and standard deviation of the time-dependent Monte Carlo energy distribution

std_80_km = zeros(size(f_E_MC_avg_80_km,1),1); 
for i = 1:size(f_E_MC_avg_80_km,1)
    std_80_km(i) = sqrt(sum((f_E_MC_80_km(i,ind_eq_80_km:end) - ...
        f_E_MC_avg_80_km(i)).^2)/(length(avg_t_80_km(ind_eq_80_km:end))-1)); 
end 

std_100_km = zeros(size(f_E_MC_avg_100_km,1),1); 
for i = 1:size(f_E_MC_avg_100_km,1)
    std_100_km(i) = sqrt(sum((f_E_MC_100_km(i,ind_eq_100_km:end) - ...
        f_E_MC_avg_100_km(i)).^2)/(length(avg_t_100_km(ind_eq_100_km:end))-1)); 
end


std_120_km = zeros(size(f_E_MC_avg_120_km,1),1); 
for i = 1:size(f_E_MC_avg_120_km,1)
    std_120_km(i) = sqrt(sum((f_E_MC_120_km(i,ind_eq_120_km:end) - ...
        f_E_MC_avg_120_km(i)).^2)/(length(avg_t_120_km(ind_eq_120_km:end))-1)); 
end

%% Compute error bar

% Lower and upper errorbars
L_80_km  = f_E_MC_avg_80_km  - std_80_km;
L_100_km = f_E_MC_avg_100_km - std_100_km;
L_120_km = f_E_MC_avg_120_km - std_120_km;
for i = 1:length(f_E_MC_avg_80_km)
    if L_80_km(i) < 0
        L_80_km(i) = max(1e-1,f_E_MC_avg_80_km(i) - std_80_km(i));
    end 
    if L_100_km(i) < 0
        L_100_km(i) = max(1e-1,f_E_MC_avg_100_km(i) - std_100_km(i));
    end 

    if L_120_km(i) < 0
        L_120_km(i) = max(1e-1,f_E_MC_avg_120_km(i) - std_120_km(i));
    end  
end 
U_80_km  = f_E_MC_avg_80_km  + std_80_km;
U_100_km = f_E_MC_avg_100_km + std_100_km;
U_120_km = f_E_MC_avg_120_km + std_120_km;

%% Solution of Boltzmann equation and 80, 100 and 120 km with ERP 200 MW
z0_3         = [80e3 100e3 120e3]; 
E0_field_3   = [1.4 1.12 0.93];
f_E_stubbe_3 = zeros(7999,length(z0_3));
for j = 1:length(z0_3)

        % Call MATLAB MSIS
        [rho,temperature] = neutral_output_MSIS(z0_3(j)); 

        % Neutral densities MATLAB MSIS
        nO = rho(:,2)';   % [1/m^3] 
        nN2 = rho(:,3)';    % [1/m^3] 
        nO2 = rho(:,4)';    % [1/m^3]

        % Max electron energy [eV]
        max_E = 6.0; 

        % Compute non-Maxwellian distribution 
        [f0_exp_stubbe,v_bins_stubbe,Te_nm_2nd,~,A,~,~,~] = ...
            compute_f0_stubbe_less_xs(E0_field_3(j),...
                                      n_e,Tn,nO,nN2,nO2,...
                                      f_HF,xs_data,max_E);

        % Energy/speed in the middle fo the bin
        v_bin_mid = v_bins_stubbe(1:end-1)/2 + v_bins_stubbe(2:end)/2; % [m/s]
        E_bin_mid_eV = (me*v_bin_mid.^2)/(2*eVtoJ); % [eV]

        % Speed distribution [#e/m^3 * s/m]
        f0_speed_distri = n_e*A*4*pi*v_bin_mid.^2.*f0_exp_stubbe;

        % Convert speed distribution to energy distribution [#e/m^3 * 1/eV]
        f_E_stubbe_3(:,j) = convert_fv2fE(v_bins_stubbe,f0_speed_distri);

end 


%% Second moment temperature of MC 
Te_non_max_MC_120_km = (4*pi*me)/(3*k_b) * ...
                       trapz(v_bin_mid_MC',...
                             (v_bin_mid_MC'.^2.*f_v_MC_avg_120_km)/(n_e*4*pi));
Te_non_max_MC_100_km = (4*pi*me)/(3*k_b) * ...
                       trapz(v_bin_mid_MC',...
                             (v_bin_mid_MC'.^2.*f_v_MC_avg_100_km)/(n_e*4*pi));
Te_non_max_MC_80_km = (4*pi*me)/(3*k_b) * ...
                      trapz(v_bin_mid_MC',...
                            (v_bin_mid_MC'.^2.*f_v_MC_avg_80_km)/(n_e*4*pi));

%% Plot figure 5: Solution of Boltzmann equation and MC run 
% at 80, 100 and 120 km with error bars (Jackife)

% Some colors
col_oransje = [255 128   0]./255; 
col_bla =     [  0   0 255]./255; 
col_gron =    [  0 204   0]./255; 

figure
subplot(1,3,1)
phS_80 = loglog(E_bin_mid_eV,f_E_stubbe_3(:,1),...
                '-.',...
                'LineWidth',2.0,...
                'Color',col_bla);
hold on; 

[Pa_80_km,Li_80_km] = JackKnife(E_bin_mid_conv_MC(f_E_MC_avg_80_km(:)>0)',....
                    f_E_MC_avg_80_km(f_E_MC_avg_80_km(:)>0)',...
                    L_80_km(f_E_MC_avg_80_km(:)>0)',...
                    U_80_km(f_E_MC_avg_80_km(:)>0)',...
                    col_gron,col_oransje);
hold off;
set(gca,...
    'fontsize',font_a,...
    'box','on',...
    'tickdir','both',...
    'TickLength',[0.008 0.008])
legend([phS_80,Li_80_km,Pa_80_km],...
       'Stubbe','MC-run','\pm \sigma',...
       'Location','southwest',...
       'FontSize',font_lg)
title('80 km','FontSize',font_t)
xlabel('Energy [eV]','FontSize',font_lb)
ylabel('f_0(E) [m^{-3} eV^{ -1}]','FontSize',font_lb)
xlim([0.05 3])
xticks([0.1 0.5 1 2])
ylim([5e2 1e12])
set(gca,'xscale','log')
set(gca,'yscale','log')

subplot(1,3,2)
phS_100 = loglog(E_bin_mid_eV,f_E_stubbe_3(:,2),...
                 '-.',...
                 'LineWidth',2.0,...
                 'Color',col_bla);
hold on; 

[Pa_100_km,Li_100_km] = JackKnife(E_bin_mid_conv_MC(f_E_MC_avg_100_km(:)>0)',....
                    f_E_MC_avg_100_km(f_E_MC_avg_100_km(:)>0)',...
                    L_100_km(f_E_MC_avg_100_km(:)>0)',...
                    U_100_km(f_E_MC_avg_100_km(:)>0)',...
                    col_gron,col_oransje);
hold off;
set(gca,...
    'fontsize',font_a,...
    'box','on',...
    'tickdir','both',...
    'TickLength',[0.008 0.008])
legend([phS_100,Li_100_km,Pa_100_km],...
       'Stubbe','MC-run','\pm \sigma',...
       'Location','southwest',...
       'FontSize',font_lg)
title('100 km','FontSize',font_t)
xlabel('Energy [eV]','FontSize',font_lb)
xlim([0.05 3])
xticks([0.1 0.5 1 2])
ylim([5e2 1e12])
set(gca,'xscale','log')
set(gca,'yscale','log')

subplot(1,3,3)
phS_120 = loglog(E_bin_mid_eV,f_E_stubbe_3(:,3),...
                 '-.','LineWidth',2.0,'Color',col_bla);
hold on; 

[Pa_120_km,Li_120_km] = JackKnife(E_bin_mid_conv_MC(f_E_MC_avg_120_km(:)>0)',....
                    f_E_MC_avg_120_km(f_E_MC_avg_120_km(:)>0)',...
                    L_120_km(f_E_MC_avg_120_km(:)>0)',...
                    U_120_km(f_E_MC_avg_120_km(:)>0)',...
                    col_gron,col_oransje);
hold off;
set(gca,...
    'fontsize',font_a,...
    'box','on',...
    'tickdir','both',...
    'TickLength',[0.008 0.008])
legend([phS_120,Li_120_km,Pa_120_km],...
       'Stubbe','MC-run','\pm \sigma',...
       'Location','southwest',...
       'FontSize',font_lg)
title('120 km','FontSize',font_t)
xlabel('Energy [eV]','FontSize',font_lb)
xlim([0.05 3])
xticks([0.1 0.5 1 2])
ylim([5e2 1e12])
set(gca,'xscale','log')
set(gca,'yscale','log')

%% Plot cross sections in figure A1

figure
subplot(1,3,1)
ph1 = semilogy(repmat(E_bin_mid_eV,1,8),[xs_e_N2vib0_1 xs_e_N2vib0_2 ...
    xs_e_N2vib0_3 xs_e_N2vib0_4 xs_e_N2vib0_5 xs_e_N2vib0_6 ...
    xs_e_N2vib0_7 xs_e_N2vib0_8],'LineWidth',2.5); 
cmp1 = colormap('winter'); 
cmlines(ph1,cmp1)
hold on;
semilogy(E_old*(1/eVtoJ),xs_old,'-.','LineWidth',2.5,'color',col60); 
hold off; 
legend('1','2','3','4','5','6','7','8','9','Location','northwest')
xlim([0 3])
xticks([0 1 2 3])
ylim([2e-24 1e-18])
yticks([1e-24 1e-23 1e-22 1e-21 1e-20 1e-19 1e-18])
xlabel('Energy [eV]','FontSize',font_lb)
ylabel('Cross sections [m^2]','FontSize',font_lb)
title('a)','FontSize',font_t)
set(gca,...
    'fontsize',font_a,...
    'box','on',...
    'tickdir','both','TickLength',[0.007 0.004])

col1 = [0,0,204]/255; 
col2 = [0,150,220]/255; 
col3 = [116,201,0]/255; 
col4 = [128,245,0]/255; 
subplot(1,3,2)
ph5 = semilogy(repmat(E_bin_mid_eV,1,3),...
               [xs_e_N2elastic xs_e_O2elastic xs_e_Oelastic],...
               'LineWidth',2.5); 
cmp5 = colormap('copper'); 
cmlines(ph5,cmp5)
hold on;
semilogy(E_bin_mid_eV,xs_e_N2rot0_2,'-.','LineWidth',2.5,'color',col1); 
semilogy(E_bin_mid_eV,xs_e_N2rot0_4,'-.','LineWidth',2.5,'color',col2); 
semilogy(E_bin_mid_eV,xs_e_N2rot0_6,'-.','LineWidth',2.5,'color',col3); 
semilogy(E_bin_mid_eV,xs_e_N2rot0_8,'-.','LineWidth',2.5,'color',col4); 
hold off; 
legend('1','2','3','4','5','6','7','Location','southeast')
xlim([0 3])
xticks([0 1 2 3])
ylim([2e-24 1e-18])
yticks([1e-24 1e-23 1e-22 1e-21 1e-20 1e-19 1e-18])
xlabel('Energy [eV]','FontSize',font_lb)
title('b)','FontSize',font_t)
set(gca,...
    'fontsize',font_a,...
    'box','on',...
    'tickdir','both',...
    'TickLength',[0.007 0.004])
subplot(1,3,3)
ph3 = semilogy(repmat(E_bin_mid_eV,1,3),...
               [xs_e_Ofine_2_1 xs_e_Ofine_2_0 xs_e_Ofine_1_0],...
               'LineWidth',2.5);
cmp3 = colormap('cool'); 
cmlines(ph3,cmp3)
hold on; 
ph4 = semilogy(repmat(E_bin_mid_eV,1,4),...
               [xs_e_O2a1Dg xs_e_O2b1Sgp xs_e_O3P1D xs_e_O3P1S],...
               ':','LineWidth',2.5);
cmp4 = colormap('autumn'); 
cmlines(ph4,cmp4)
hold off; 
legend('1','2','3','4','5','6','Location','northeast')
xlim([0 3])
xticks([0 1 2 3])
ylim([2e-24 1e-18])
yticks([1e-24 1e-23 1e-22 1e-21 1e-20 1e-19 1e-18])
xlabel('Energy [eV]','FontSize',font_lb)
title('c)','FontSize',font_t)
set(gca,...
    'fontsize',font_a,...
    'box','on',...
    'tickdir','both',...
    'TickLength',[0.007 0.004])



