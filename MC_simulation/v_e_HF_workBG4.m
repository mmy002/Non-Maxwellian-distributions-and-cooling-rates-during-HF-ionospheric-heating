%% Script for MC-trajectory calculations
% Includes only collisions between electrons and neutrals 
% Runs for 1 milliions electrons (OBS: takes a lot of time...)

%% Physical constants
c0 = 2.99792458e8;             % Speed of light [m/s]
h = 6.62607015e-34;            % Plank's constant [Js]
kB = 1.380649e-23;             % Boltzmann constant [J/K]   
m_e = 9.1093835611e-31;        % electron rest mass [kg]
q_e = 1.602176620898e-19;      % elementary charge [C]

%% Geophysical and ionospheric
B = 5e-05;         % Magnetic field-strength [T]
n_e = 1e10;        % electron density [/m^3]
nu_e = 1e6;        % electron collision-frequency [/s]
w_e = w_e_gyro(B); % electron gyro-frequency [rad/s]
T_e = 300;         % electron temperature [K]

%% Heating
w_HF = 2*pi*4.6e6; % HF-frequency of Heating [rad/s]
ERP = 200e6;       % Effective radiative power [W]
z0 = 80e3;         % Altitude [m]

%% Radio-wave-propagation etc
beta_HF = 1;                                       % Reduction-factor due to absorption up to z0, 1 <-> no loss
epsilon1 = dielectric_T_hf(w_HF,n_e,nu_e,B);       % HF Dielectric tensor at hight
sigma1 = conductivity_T_hf(w_HF,n_e,nu_e,B);       % HF conductivity tensor
n_refr1 = n_refractive_indx(w_HF,sigma1,epsilon1); % Refractive index
[Eo1_0,Eo1_0b] = E0_HF_pump(n_refr1{1},1,ERP,z0);  % E-fields at height [V/m]
theta2B = 12/180*pi;                               % Angle of k-vector to B (rad)

%% Load data for cross sections
xs_data = get_data_withN2vib(); 

%% Call MATLAB MSIS for neutral densities
[rho,~] = neutral_output_MSIS(z0); 

% Neutral densities MATLAB MSIS
nO = rho(:,2)';   % [1/m^3] 
nN2 = rho(:,3)';    % [1/m^3] 
nO2 = rho(:,4)';    % [1/m^3]

Nn = [nN2; nO2; nO;
      nN2; nN2; nN2; nN2; nN2; nN2; nN2; nN2;
      nN2; nN2; nN2; nN2; ...
      nO2; nO2;
      nO; nO; nO; nO; nO;
      nO2]; 

%% Modeling with collisions, using parfor to parallelize the workings
% This speeds things up with a factor of 40 or thereabouts
% downside is that we get the results in 100 separate files with 1e5
% velocity evolutions in each file
n_p = 10;%10000; % Number of electron for each loop
n_colls = 300; 
filename = 'BG3-run1';
delete(gcp('nocreate'))
pool = parpool(48)

% Excitation energy for inelastic collisions [eV] from Itakawa 1986, 1989, 1990
exE_N2elastic = 0; 
exE_O2elastic = 0; 
exE_Oelastic = 0; 
exE_N2vib0_1 = 0.2888; 
exE_N2vib0_2 = 0.5742; 
exE_N2vib0_3 = 0.8559; 
exE_N2vib0_4 = 1.1342; 
exE_N2vib0_5 = 1.4088; 
exE_N2vib0_6 = 1.6800; 
exE_N2vib0_7 = 1.9475; 
exE_N2vib0_8 = 2.2115; 
exE_N2rot0_2 = 2*2.4668e-4;
exE_N2rot0_4 = 4*2.4668e-4; 
exE_N2rot0_6 = 6*2.4668e-4; 
exE_N2rot0_8 = 8*2.4668e-4; 
exE_O2a1Dg = 0.9770; 
exE_O2b1Sgp = 1.6270; 
exE_Ofine_1_0 = 0.0085; 
exE_Ofine_2_0 = 0.0281; 
exE_Ofine_2_1 = 0.0196; 
exE_O1S = 4.1900; 
exE_O1D = 1.9670; 
exE_O2vib = 0.193; 
dE = [exE_N2elastic; exE_O2elastic; exE_Oelastic;
      exE_N2vib0_1; exE_N2vib0_2; exE_N2vib0_3; exE_N2vib0_4; exE_N2vib0_5; exE_N2vib0_6; exE_N2vib0_7; exE_N2vib0_8;
      exE_N2rot0_2; exE_N2rot0_4; exE_N2rot0_6; exE_N2rot0_8;
      exE_O2a1Dg; exE_O2b1Sgp;
      exE_Ofine_1_0; exE_Ofine_2_0; exE_Ofine_2_1; exE_O1S; exE_O1D;
      exE_O2vib]; 

% Folder stuff
t_start = clock;
dirname = sprintf('MC-Run-2-%s',sprintf('UT-%d%02d%02d-%02d%02d%02d',floor(clock))); 
mkdir(dirname) % Create new folder
addpath(dirname) % Add created folder to matlab path 
cd(dirname) % Go to folder with data for the MC run
pwd 

% Max number of iterations 
max_iter = 1;%100; 

parameters = cell(2,3); 
parameters{1,1} = '#electrons'; 
parameters{1,2} = '#coll'; 
parameters{1,3} = 'height_km'; 
parameters{2,1} = n_p*max_iter; 
parameters{2,2} = n_colls; 
parameters{2,3} = z0/1000; 

tic
for i_iter = 1:max_iter

    if rem(i_iter,1)==0
    fprintf('Starting with loop %d at %s, %2.2f done\n',...
          i_iter,...
          datestr(now,'HH:MM:SS'),...
          (i_iter-1)/100)
    end

    t_all = zeros(n_p,n_colls+1);
    E_all = zeros(n_p,n_colls+1);
    v_pars = zeros(n_p,n_colls+1);
    v_perps = zeros(n_p,n_colls+1);
    v_x = zeros(n_p,n_colls+1);
    v_y = zeros(n_p,n_colls+1);
    dnCol = 250;

    %for iP = 1:n_p
    parfor (iP = 1:n_p, 48 )
        [t_a,E_a,v_par,v_perp,v_X,v_Y] = f_e_workerBG3(n_colls*dnCol,B,Eo1_0,w_HF,theta2B,T_e,Nn,dE,xs_data);
        t_all(iP,:)   = t_a(1:dnCol:end);
        E_all(iP,:)   = E_a(1:dnCol:end);
        v_pars(iP,:)  = v_par(1:dnCol:end);
        v_perps(iP,:) = v_perp(1:dnCol:end);
        v_x(iP,:)     = v_X(1:dnCol:end);
        v_y(iP,:)     = v_Y(1:dnCol:end);
    end 

    avg_t = mean(t_all);
    currfilename = sprintf('%s-%04d.mat',filename,i_iter);
    save(currfilename,...
    't_all',...
    'E_all',...
    'v_pars',...
    'v_perps',...
    'v_x',...
    'v_y',...
    'dE',...
    'avg_t',...
    'parameters');

end
toc
