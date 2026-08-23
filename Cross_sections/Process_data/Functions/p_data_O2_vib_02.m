function [E_old,xs_old] = p_data_O2_vib_02()

% Convert from cm to m
cm_to_m = 0.01; 

% Convert from [eV] to [J] 
eVtoJ = 1.602176634e-19;  

delta_E_th_itiawka = 50e-3; % [J]

% Table 1 in Allan 1995
energy_allan = [0.214 0.338 0.460 0.579 0.696 0.812 0.925 1.036 1.144 1.251 1.355 1.458 1.558 1.656 ...
    1.752 1.846 1.937 2.026 2.113 2.197]'*eVtoJ; % J
O2_vib_allan = ([0 0 1.6 40 88 95 76 47 25 12 5.6 2.4 0.94 0.6 0 0 0 0 0 0]/delta_E_th_itiawka)'...
    *1e-20*cm_to_m^2;

% Add zeros between data
num_E = 1; 
dE = zeros(length(O2_vib_allan),1); 
kk = zeros(20,1); 
hh = zeros(20,1); 
k = -1; 
h = 0;
E_allan_zero_pad = zeros(length(O2_vib_allan)+(length(O2_vib_allan)),1); 
O2_vib_allan_zero_pad = zeros(length(O2_vib_allan)+(length(O2_vib_allan)),1); 

% Find resolution of the energy, which varies
for j =1:length(O2_vib_allan)-1
    dE(j) = energy_allan(j+1) - energy_allan(j); % Delta E   
end 
dE(end) = dE(end-1); 

%Insert zeros between data points 
for i = 1:length(O2_vib_allan)
        k = k+2; % Index of oddetall, where the energy from the data goes 
        kk(i) = k; 
        h = h+2; % index of partall, where new energies is put between then original energies -> zero for xs
        hh(i) = h; 

    % Orginal data 
    E_allan_zero_pad(k) = energy_allan(i);  
    O2_vib_allan_zero_pad(k) = O2_vib_allan(i);

    % -> zeros inserted between cross sectio data 
    E_allan_zero_pad(h) = (linspace(energy_allan(i),energy_allan(i)+dE(i)/2,num_E))'; 
    O2_vib_allan_zero_pad(h) = 0; 

end


E_old = E_allan_zero_pad;
xs_old = O2_vib_allan_zero_pad; 

end