function L_sum = L_stubbe_tot(v_bins_stubbe,dIe,deg_Ie)
% Compute non-Maxwellian cooling rates from inelatic collision integral

%% Constant
me = 9.10938291e-31; % Electron mass [kg]

%% Non-Maxwellian cooling rates 

% Degrading electrons due to collisions: Sum of loss and gain [#e/(m^3*s)]  
sum_dIe_N2_vib  = sum(dIe(:,1:10),2) + sum(deg_Ie(:,1:10),2); 
sum_dIe_O2a1Dg  = sum(dIe(:,11),2)   + sum(deg_Ie(:,11),2);
sum_dIe_O2b1Sgp = sum(dIe(:,12),2)   + sum(deg_Ie(:,12),2);   
sum_dIe_O_fs_10 = sum(dIe(:,13),2)   + sum(deg_Ie(:,13),2); 
sum_dIe_O_fs_20 = sum(dIe(:,14),2)   + sum(deg_Ie(:,14),2); 
sum_dIe_O_fs_21 = sum(dIe(:,15),2)   + sum(deg_Ie(:,15),2); 
sum_dIe_O3P1D   = sum(dIe(:,16),2)   + sum(deg_Ie(:,16),2); 
sum_dIe_O3P1S   = sum(dIe(:,17),2)   + sum(deg_Ie(:,17),2);
sum_dIe_O2_vib  = sum(dIe(:,18:end),2) + sum(deg_Ie(:,18:end),2); 

% Energy/speed in the middle fo the bin
v_bin_mid = v_bins_stubbe(1:end-1)/2 + v_bins_stubbe(2:end)/2; % [m/s]
E_bin_mid = (me*v_bin_mid.^2)/2; % [J]

% J/(m^3*s)
L_N2_vib  = sum(-E_bin_mid.*sum_dIe_N2_vib);
L_O2a1Dg  = sum(-E_bin_mid.*sum_dIe_O2a1Dg); 
L_O2b1Sgp = sum(-E_bin_mid.*sum_dIe_O2b1Sgp);
L_O_fs_10 = sum(-E_bin_mid.*sum_dIe_O_fs_10); 
L_O_fs_20 = sum(-E_bin_mid.*sum_dIe_O_fs_20); 
L_O_fs_21 = sum(-E_bin_mid.*sum_dIe_O_fs_21);  
L_O3P1D   = sum(-E_bin_mid.*sum_dIe_O3P1D); 
L_O3P1S   = sum(-E_bin_mid.*sum_dIe_O3P1S);
L_O2_vib  = sum(-E_bin_mid.*sum_dIe_O2_vib); 
L_sum     = L_N2_vib + ...
	    L_O2_vib + ...
	    L_O2a1Dg + L_O2b1Sgp + ...
	    L_O_fs_10 + L_O_fs_20 + L_O_fs_21 + ...
	    L_O3P1D + L_O3P1S;  


end 
