function f0_E_stubbe = convert_fv2fE(v_bins_stubbe,f0_v_stubbe)
% This function converts from a speed distribution to a energy distribution

me = 9.10938291e-31;        % Electron mass [kg]
eVtoJ = 1.602176634e-19;    % Convert from [eV] to [J] 

% Convert speed bins to energy bins 
E_bins_stubbe = (me*v_bins_stubbe.^2)/2; % [J]

% Resolution of speed(constant) [m/s]
delta_v_stubbe = v_bins_stubbe(2)-v_bins_stubbe(1); % Resolution of speed(constant) [m/s]

% Resolution of energy (varying)
delta_E_eV = zeros(length(E_bins_stubbe)-1,1); 
for i = 1:length(E_bins_stubbe)-1
  delta_E_eV(i,1) = (E_bins_stubbe(i+1)-E_bins_stubbe(i))/eVtoJ; % [eV]
end 


% Convert speed distribution to energy distribution [#e/m^3 * 1/eV]
f0_E_stubbe = (f0_v_stubbe*delta_v_stubbe)./delta_E_eV; 

end 
