function [L_vib_N2,L_vib_O2,L_oxygen_10,L_oxygen_20,L_oxygen_21,L_rot_N2,L_rot_O2,...
	  L_elas,L_delta_O2,L_sigma_O2] = cooling_rates_tables_pavlov2(Te,Tn,nNe,nN2,nO2,nO)
% This functions models cooling rates of electron colliding 
% with the neutral species molecular nitrogen, molecular oxygen 
% and atomic oxygen, both elastic and inelastic collisions. 
% NEW: add cooling rates for electronic excitation of O_2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: Electron and neutral temperature in [K] and electron 
% and neutral densities in [m^-3 ]. Row/column vectors or scalars.  
% .........................................................................
% Output: Cooling rates between electrons and different neutrals [eV/(cm^3*s)] 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert densities from m^-3 to cm^-3
m_neg3_to_cm_neg3 = 1e-6; 
nNe_cm_neg3       = nNe * m_neg3_to_cm_neg3;  
nN2_cm_neg3       = nN2 * m_neg3_to_cm_neg3; 
nO2_cm_neg3       = nO2 * m_neg3_to_cm_neg3; 
nO_cm_neg3        = nO  * m_neg3_to_cm_neg3; 

%% Vibrational and rotational excitation of N2: 
% From Pavlov (1998a): " New electron energy transfer rates for
% vibrational excitation of N2". Ann. Geophysicea 16, 176-182. 

E1       = 3353; % Energy of first vibrational level of N2, [K]
T_vib_N2 = Tn;   % Assume that the vibrational temperature is equal to neutral temperature.

% Implement eq. 19-> log(Q_0v) and 20-> log(Q_1v) from Pavlov (1998a)
log_Q_0v = zeros(length(1:10),1); 
log_Q_1v = zeros(length(2:9),1); 

% Coefficents from table 3 of Pavlov (1998a)
coef_matrix_Q_1v = [2  -3.413 7.326e-3 -2.200e-6 3.128e-10 -1.702e-14;
		    3  -4.160 7.803e-3 -2.352e-6 3.352e-10 -1.828e-14; 
		    4  -5.193 8.360e-3 -2.526e-6 3.606e-10 -1.968e-14;
		    5  -5.939 8.807e-3 -2.669e-6 3.806e-10 -2.073e-14;
		    6  -8.261 1.010e-2 -3.039e-6 4.318e-10 -2.347e-14;
		    7  -8.185 1.010e-2 -3.039e-6 4.318e-10 -2.347e-14;
		    8 -10.823 1.199e-2 -3.620e-6 5.159e-10 -2.810e-14;
		    9 -11.273 1.283e-2 -3.879e-6 5.534e-10 -3.016e-14];

for v1 = 1:8
  log_Q_1v(v1) = coef_matrix_Q_1v(v1,2) + ...
		 coef_matrix_Q_1v(v1,3) * Te + ...
		 coef_matrix_Q_1v(v1,4) * (Te)^2 + ...
		 coef_matrix_Q_1v(v1,5) * (Te)^3 + ...
		 coef_matrix_Q_1v(v1,6) * (Te)^4 - ...
		 16;  % Q_1v have unit [eVcm^3s^-1]
end

Q_1v = (10.^(log_Q_1v)); %[eVcm^3s^-1]
%Q_1v_out = Q_1v'; 

BB = zeros(8,1); 

for j = 1:8
  BB(j) = Q_1v(j)*(1 - exp(j*E1*((1/Te) - (1/T_vib_N2)))); 
end 

sum_BB = sum(BB); % sum from v = 2 to 9 for Q_1v

if Te <= 1500
  % Coefficents from table 2 of Pavlov (1998a)
  coef_matrix_Q_0v_2 = [1 -6.426 3.151e-2 -4.075e-5 2.439e-8 -5.497e-12]; 
  
  log_Q_0v(1) = coef_matrix_Q_0v_2(2) + ...
		coef_matrix_Q_0v_2(3) * Te + ...
		coef_matrix_Q_0v_2(4) * (Te)^2 + ...
		coef_matrix_Q_0v_2(5) * (Te)^3 + ...
		coef_matrix_Q_0v_2(6) * (Te)^4 - ...
		16;  % Q_0v have unit [eVcm^3s^-1]
  
  Q_0v(1) = (10.^(log_Q_0v(1))); %[eVcm^3s^-1]
  Q_0v(2:10) = zeros(9,1); 
  %Q_0v_out = Q_0v'; 
  
  AA = zeros(1,1); 
  for i = 1
    AA(i) = Q_0v(i)*(1 - exp(i*E1*((1/Te) - (1/T_vib_N2)))); 
  end 
  sum_AA = sum(AA); % sum from v = 1 for Q_0v_2
  
  % Electron cooling rate for Boltzmann distribution. [eV/cm^3s^1]: Equation
  % 11 from paper: 
  L_vib_N2 = nNe_cm_neg3 * nN2_cm_neg3 * (1-exp(-E1/T_vib_N2)) * sum_AA + ...
	     nNe_cm_neg3 * nN2_cm_neg3 * (1-exp(-E1/T_vib_N2)) * exp(-E1/T_vib_N2) * sum_BB; 
  
else 
  
  % Coefficents from table 1 of Pavlov (1998a)
  coef_matrix_Q_0v = [1    2.025 8.782e-4  2.954e-7 -9.562e-11  7.252e-15;
		      2   -7.066 1.001e-2 -3.066e-6  4.436e-10 -2.449e-14;
		      3   -8.211 1.092e-2 -3.369e-6  4.891e-10 -2.706e-14;
		      4   -9.713 1.204e-2 -3.732e-6  5.431e-10 -3.008e-14;
		      5  -10.353 1.243e-2 -3.850e-6  5.600e-10 -3.100e-14;
		      6  -10.819 1.244e-2 -3.771e-6  5.385e-10 -2.936e-14;
		      7  -10.183 1.185e-2 -3.570e-6  5.086e-10 -2.769e-14;
		      8  -12.698 1.309e-2 -3.952e-6  5.636e-10 -3.071e-14;
		      9  -14.710 1.409e-2 -4.249e-6  6.058e-10 -3.300e-14;
		      10 -17.538 1.600e-2 -4.916e-6  7.128e-10 -3.941e-14]; 
  for v0 = 1:10
    log_Q_0v(v0) = coef_matrix_Q_0v(v0,2) + ...
		   coef_matrix_Q_0v(v0,3) * Te + ...
		   coef_matrix_Q_0v(v0,4) * (Te)^2 + ...
		   coef_matrix_Q_0v(v0,5) * (Te)^3 + ...
		   coef_matrix_Q_0v(v0,6) * (Te)^4 - ...
		   16;  % Q_0v have unit [eVcm^3s^-1]
  end 
  
  Q_0v = (10.^(log_Q_0v)); %[eVcm^3s^-1]
  %Q_0v_out = Q_0v'; 
  
  AA = zeros(10,1); 
  for i = 1:10
    AA(i) = Q_0v(i)*(1 - exp(i*E1*((1/Te) - (1/T_vib_N2)))); 
  end 
  sum_AA = sum(AA); % sum from v = 1 to 10 for Q_0v
  
  % Electron cooling rate for Boltzmann distribution. [eV/cm^3s^1]: Equation
  % 11 from Pavlov (1998a)
  L_vib_N2 = nNe_cm_neg3 * nN2_cm_neg3 * (1-exp(-E1/T_vib_N2)) * sum_AA + ...
	     nNe_cm_neg3 * nN2_cm_neg3 * (1-exp(-E1/T_vib_N2)) * exp(-E1/T_vib_N2) * sum_BB; 
  
end 

%% Fine structures exication of atomic oxygen: 
% From Pavlov and Berringston (1999): "Cooling rates of thermal electrons by 
% electron impact excitation of fine structure levels of atomic oxygen". 
% Ann. Geophysicea 17, 919-924. 

D   = 5 + exp(-326.6*(1./Tn)) + 3*exp(-227.7*(1./Tn));
S10 = (8.249*1e-16*Te.^0.6.*exp(-227.7*(1./Tn))); % [eVcm^3s^-1]
S20 = (1.191*1e-11);                              % [eVcm^3s^-1]
S21 = (1.863*1e-11);                              % [eVcm^3s^-1]

L_oxygen_10 = nNe_cm_neg3 .* nO_cm_neg3.*(1./D) .* (S10.*(1-exp(98.9*((1./Te)-(1./Tn))))); % [eVcm^-3s^-1]
L_oxygen_20 = nNe_cm_neg3 .* nO_cm_neg3.*(1./D) .* (S20*(1-exp(326.6*((1./Te)-(1./Tn))))); % [eVcm^-3s^-1]
L_oxygen_21 = nNe_cm_neg3 .* nO_cm_neg3.*(1./D) .* (S21*(1-exp(227.7*((1./Te)-(1./Tn))))); % [eVcm^-3s^-1]

%L_oxygen_15 = [L_oxygen_20; L_oxygen_21; L_oxygen_10]; 

%% Vibrational and rotational excitation of O2: 
% From Pavlov (1998b): " New electron energy transfer and cooling rates for 
% vibrational excitation of O2", Ann. Geophysicea 16, 1007-1013. 

% Energy of the first vibrational level of O2
E1_O2 = 2239; % [K] 

% Table 1 from Pavlov (1998b)
coef_matrix_Q0v_02_Pavlov = [1, 8.56e-15, 303, 10.00, -0.20, 1.05e-3, 3150; 
			     2, 1.15e-18, 305, 18.78, -0.25, 9.24e-4, 3450; 
			     3, 1.77e-23, 301, 29.00, -0.20, 6.16e-4, 3150;
			     4, 7.05e-25, 301, 30.07,  0.31, 1.00e-3, 1800;
			     5, 2.14e-28, 299, 37.52, -0.42, 5.28e-4, 3000; 
			     6, 2.94e-31, 300, 43.27, -0.39, 7.85e-4, 3000;
			     7, 1.08e-35, 299, 53.10, -0.90, 5.70e-4, 2900]; 

% Eq. 11 from Pavlov (1998b)
Q_0v_O2_pavlov = [0 0 0 0 0 0 0];
for vv = 1:7
  Q_0v_O2_pavlov(vv) = coef_matrix_Q0v_02_Pavlov(vv,2) * ...
		       exp((1 - coef_matrix_Q0v_02_Pavlov(vv,3)*(1/Te)) * ...
			   (coef_matrix_Q0v_02_Pavlov(vv,4) + ...
			    coef_matrix_Q0v_02_Pavlov(vv,5) * sin(coef_matrix_Q0v_02_Pavlov(vv,6) * ...
								  (Te - coef_matrix_Q0v_02_Pavlov(vv,7)))));
end 

% With Q_0v from Pavlov (1998b)
A_O2_pavlov = [0 0 0 0 0 0 0]; 
for g = 1:7
  A_O2_pavlov(g) = Q_0v_O2_pavlov(g)*(1-exp(g*E1_O2*((1/Te) - (1/Tn)))); 
end 
sum_A_O2_pavlov = sum(A_O2_pavlov); % sum from v = 1 to 4 for Q_0v

% Eq. 8. from Pavlov (1998b)
L_vib_O2 = nNe_cm_neg3*nO2_cm_neg3*sum_A_O2_pavlov; % [eVcm^-3s^-1] 

%% Rotational excitation of O2 with eq. 16 from Pavlov (1998b)
C_O2 = 5.2e-15; % [eVcm^3s^-1K^-0.5]
L_rot_O2 = C_O2.*nO2_cm_neg3.*nNe_cm_neg3.*(Te - Tn).*Te.^(-0.5); % [eVcm^-3s^-1] 

%% Rotational excitation of N2: We use eq. A2 from Pavlov (1998a)
C = 3.51e-14; % [eVcm^3s^-1K^-0.5]
L_rot_N2 = C.*nN2_cm_neg3.*nNe_cm_neg3.*(Te-Tn).*Te.^(-0.5); % [eVcm^-3s^-1] 

%% Electronic excitation of O_2: Units are in [eVcm^-3s^-1]. 

f1           = (13200 + 1410*sin(0.000241*(Te-500)))*(1 + exp((Te-14011)/1048)); 
G_Te         = 1.143e-14*exp(f1*(1/1500 - 1/Te)); 
E_1_sigma_O2 = 0.977; 
g1_g0        = 5/3; 
L_delta_O2   = nNe_cm_neg3.*nO2_cm_neg3*G_Te*(1 - exp(-11400/Tn)*g1_g0*exp(E_1_sigma_O2/Te)); 

f2           = (19225 + 560*sin(0.000383*(Te-1000)))*(1 + exp((Te-16382)/1760)); 
H_Te         = 1.616e-16*exp(f2*(1/1500 - 1/Te)); 
E_2_sigma_O2 = 1.627; 
g2_g0        = 1/3; 
L_sigma_O2   = nNe_cm_neg3.*nO2_cm_neg3*H_Te*(1 - exp(-11400/Tn)*g2_g0*exp(E_2_sigma_O2/Te)); 

%% Elastic collisions: 
% From R. W. Schunk and A. F. Nagy (1978): "Electron temperature in the f-region of the ionoshpere: 
% Theory and observations". Rev. Geophys. 16., 355-399. doi:10.1029/RG016i003p00355

% Units are in [eVcm^-3s^-1]. Equation 43a, 43b, 43c from Schunk and Nagy (1978)
LN2elas = nNe_cm_neg3 .* nN2_cm_neg3 * 1.77e-19 .*Te       .*(Te-Tn) .* (1-1.21e-4.*Te);     % Elastic collisions of N2 with neutral particles
LO2elas = nNe_cm_neg3 .* nO2_cm_neg3 * 1.21e-18 .*sqrt(Te) .*(Te-Tn) .* (1+3.6e-2*sqrt(Te)); % Elastic collisions of O2 with neutral particles
LOelas  = nNe_cm_neg3 .* nO_cm_neg3  * 7.9e-19  .*sqrt(Te) .*(Te-Tn) .* (1+5.7e-4*Te);       % Elastic collisions of O with neutral particles
L_elas  = LN2elas + LO2elas + LOelas;                                                        % Sum of cooling rates for elastic collisions 


end 
