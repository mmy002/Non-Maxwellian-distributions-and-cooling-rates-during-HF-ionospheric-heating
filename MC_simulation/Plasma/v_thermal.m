function V = v_thermal(T,m)
% V_THERMAL - thermal velocity
%   
% Calling:
%  V = v_thermal(T,m)
% Input:
%  T - Temperature (K)  [n x m]
%  m - mass of particle species (kg)
% Output
%  V - average velocity (m/s)

kB  = 1.380649e-23; % Boltzmann constant [J/K]


V = (kB*T/m).^(1/2);
