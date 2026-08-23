function f_p = f_plasma(ne)
% F_PLASMA - Plasma frequency
%   
% Callling:
%   f_p = f_plasma(ne)
% Input:
%  ne - electron density profiles [N x M] (/m^3)
% Output:
%  w - angular plasma frequencies [N x M] (radians/s)

% Old values for physical constants
% Ep_0    = 1e-9/36/pi;               % Permittivity [As/Vm]
% q_e   = 1.6021773e-19;            % elementary charge [C]
% m_e   = 9.10939e-31;              % electron rest mass [kg]
% New values for physical constants
m_e   = 9.1093835611e-31;      % electron rest mass [kg]
q_e   = 1.602176620898e-19;    % elementary charge [C]
c0    = 2.99792458e8;          % Speed of light [m/s]
my_0  = 4*pi*1e-7;             % Permeability [Vs/Am]
Ep_0  = 1/(my_0*c0^2);         % Permittivity [As/Vm]



f_p = (max(ne,0)*q_e^2/Ep_0/m_e).^.5/2/pi;
