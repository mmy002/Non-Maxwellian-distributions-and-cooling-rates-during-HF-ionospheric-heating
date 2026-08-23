%% Calculate cooling rates for different heights and amplitude of the electric field: 
% Add new function to compute more correct cooling rates for N2 and O2 rot. 
% Add new cross sections for N2 and O2 rot. 

%% Geophysical and ionospheric
z0   = [70e3 80e3 90e3 100e3];  % Altitude [m]
f_HF = 4.6e6;                   % Radio wave frequency [Hz] 
w_HF = 2*pi*f_HF;               % Angular frequency of radio wave [rad/s] 
c    = 2.99792458e8;            % Speed of light [m/s]   

% Electric field amplitude [V/m]
E0_field_height = [0.5 0.75 0.8 0.9 1 1.25 1.5 1.75 2.0 2.2 2.5 3.0 3.5 4.0 7.0;...
		   0.5 0.75 0.8 0.9 1 1.25 1.5 1.75 2.0 2.2 2.5 3.0 3.5 4.0 5.0;...
		   0.3 0.5 0.75 0.8 0.9 1 1.25 1.5 1.75 2.0 2.2 2.5 3.0 3.5 4.0;...
		   0.3 0.5 0.75 0.8 0.9 1 1.25 1.5 1.75 2.0 2.2 2.5 3.0 3.5 4.0];  

% Load electron density data based on model from Baumann et al. 2013 
% at https://zenodo.org/records/5668082,
% which is called: night_data_autumn.txt
height_non_MSP = (60:1:120)*1e3; % [m]
Ne_dust_no_dust = readmatrix('night_data_autumn.txt'); 
Ne_no_MSP      = Ne_dust_no_dust(:,2); % No MSP

% Electron density at height region z0
Ne_no_MSP_z0 = zeros(1,length(z0)); 
for kk = 1:length(z0)
  ind_in = find(height_non_MSP==z0(kk)); 
  Ne_no_MSP_z0(kk) = Ne_no_MSP(ind_in);  
end

%% Load and process data for cross sections
% With more O2 vib. cross sections and new cross section from Bell
% for fs O and new cross sections for O2 and N2 rot. 
xs_data = get_data_even_more_xs(); 

%% Compute cooling rates
% in a folder with file name for one height and one electric field
filename     = 'cooling-rate';
height_unit  = 'km'; 
E_field_unit = 'V_m'; 

% Folder stuff
t_start = clock;
dirname = sprintf('Cooling-rate-run-%s',sprintf('UT-%d%02d%02d-%02d%02d%02d',floor(clock))); 
mkdir(dirname) % Create new folder
addpath(dirname) % Add created folder to matlab path 
cd(dirname) % Go to folder with data for the MC run
pwd 


tic
for j = 1:length(z0)
    for i = 1:size(E0_field_height,2)

      % Call MATLAB MSIS
      [rho,temperature] = neutral_output_MSIS(z0(j)); 
      
      % Neutral densities MATLAB MSIS
      nO = rho(:,2)';   % [1/m^3] 
      nN2 = rho(:,3)';  % [1/m^3] 
      nO2 = rho(:,4)';  % [1/m^3]

      % Electron density [m^-3]
      n_e = Ne_no_MSP_z0(j);
      
      % Neutral temperature [K]
      Tn = temperature(2); 
      
      % Max electron energy [eV]
      max_E = 6.0; 
      
      % Electric field at height j
      E0_field = E0_field_height(j,i); 
      
      % Compute non-Maxwellian distributions 
      [f0_exp_stubbe,v_bins_stubbe,Te_nm_2nd,Te_grad,A,dI_sum,dIe,deg_Ie] = ...
      compute_f0_stubbe_more_xs(E0_field,n_e,Tn,nO,nN2,nO2,f_HF,xs_data,max_E);
      
      % Cooling rates [J/(m^3*s)]
      L_stubbe = L_stubbe_tot(v_bins_stubbe,dIe,deg_Ie);     
      [L_rot_N2,L_rot_O2,L_el] = L_rot_el_new(v_bins_stubbe,f0_exp_stubbe,A,xs_data,...
					      n_e,nN2,nO2,nO,Te_nm_2nd,Tn); 
      L_loop = L_stubbe + L_rot_N2 + L_rot_O2 + L_el;
      
      % Compute electron temperature from Maxwellian cooling rates with the same total energy as non-Maxwellian
      if ~isnan(L_loop) 
        c_Te = @(Te) find_Te_LNM_LM(L_loop,Te,Tn,n_e,nO,nN2,nO2); % This is Q-L 
        try 
          Te_nm_L = fzero(c_Te,[100 6000]); % MATLAB fzero finds Te for Q(Te)-L(Te)=0
        catch
          try
            Te_nm_L = fzero(c_Te,[100 30000]); % MATLAB fzero finds Te for Q(Te)-L(Te)=0
          catch
            T_nm_L = 50; 
          end 
        end 
      elseif isnan(L_loop)
        Te_nm_L = 0; 
      end 
      
      % Call function to compute L and Q 
      [L_non_max,L_maxw_2nd,L_maxw_L,Q_non_max,Q_maxw_2nd,Q_maxw_L] = analytical_cooling_rates...
									(v_bins_stubbe,...
									 f0_exp_stubbe,...
									 Te_nm_2nd,...
									 Te_nm_L,...
									 A,...
									 n_e,...
									 nN2,...
									 nO2,...
									 Tn,...
									 xs_data);  
      
      % Save stuff
      currfilename = sprintf('%s-%d-%s-%.2f-%s.mat',filename,z0(j)/1000,height_unit,E0_field,E_field_unit);
      save(currfilename,...
           'Tn',...
           'n_e',...
           'f0_exp_stubbe',...
           'A',...
           'v_bins_stubbe',...
           'dIe',...
           'deg_Ie',...
           'Te_nm_2nd',...
           'Te_nm_L',...
           'Te_grad',...
           'L_loop',...
           'L_non_max',...
           'L_maxw_2nd',...
           'L_maxw_L',...
           'Q_non_max',...
           'Q_maxw_2nd',...
           'Q_maxw_L',...
           'max_E',...
           'E0_field');
      
    end 
end 
toc


