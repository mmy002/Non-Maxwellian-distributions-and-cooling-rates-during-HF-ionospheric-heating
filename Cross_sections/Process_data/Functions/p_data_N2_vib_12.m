function [E_old,xs_old] = p_data_N2_vib_12()

% Convert from cm to m
cm_to_m = 0.01; 

% Convert from [eV] to [J] 
eVtoJ = 1.602176634e-19;    

% Figure 1 in Cambell et al. 2004 with recommended (solid line) vib. excitations. Units: 10^-17 cm^2 per energy eV
N2_vib_v_1_2_cambell_et_al = load('N2_vib_cross_sections_v_1_2_cambell_et_al.mat'); 
E_old = (N2_vib_v_1_2_cambell_et_al.N2_vib_cross_sections_v_1_2_cambell_et_al(:,1))*eVtoJ ;
xs_old = (N2_vib_v_1_2_cambell_et_al.N2_vib_cross_sections_v_1_2_cambell_et_al(:,2))*1e-17*(cm_to_m)^2; 

% Remove negative data 
E_old(1) = []; 
xs_old(1) = []; 

E_old(1) = []; 
xs_old(1) = [];

% %%
% %%%
% 
% me = 9.10938291e-31; % Electron mass [kg]
% res_u_e = 1000; % Resolution of electron speed
% u_e_low = 0.2e5; % Lower limit of speed [m/s]
% u_e_high = 5.9e6; % Higher limit of speed [m/s]
% u_e = (linspace(u_e_low,u_e_high,res_u_e))'; % Electron speed [m/s]
% E = 0.5*u_e.^2*me; % Electron energy [J]
% E_eV = E*(1/eVtoJ); % Electron energy [eV]
% 
% % Interpolation 
% cross_sect = exp(interp1(E_old,log(xs_old),E,'pchip')); 
% 
% % Set to zero below excitation energy and above sample data 
% for i = 1:length(E)
%     if E(i)*(1/eVtoJ) < E_old*(1/eVtoJ)
%         cross_sect(i) = 0; 
%     end 
%     if E(i)*(1/eVtoJ) > E_old*(1/eVtoJ)
%         cross_sect(i) = 0; 
%     end 
% end 
% 
% figure
% plot(E_old*(1/eVtoJ),xs_old/(1e-17*(cm_to_m)^2),'r*')
% hold on; 
% plot(E_eV,cross_sect/(1e-17*(cm_to_m)^2),'b')
% hold off; 
% xlim([0 5])

end 