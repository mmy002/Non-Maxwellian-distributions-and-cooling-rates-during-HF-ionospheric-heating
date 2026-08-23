function w_p = w_plasma(ne)
% W_PLASMA - Angular plasma frequency
%   
% Callling:
%   w_p = w_plasma(ne)
% Input:
%  ne - electron density profiles [N x M] (/m^3)
% Output:
%  w - angular plasma frequencies [N x M] (radians/s)

% older value of epsilon_0
% Ep_0    = 1e-9/36/pi;               % Permittivity [As/Vm]

% Newer value as of 2022
c0  = 2.99792458e8;                   % Speed of light [m/s]
my_0    = 4*pi*1e-7;                  % Permeability [Vs/Am]
Ep_0    = 1/(my_0*c0^2);              % Permittivity [As/Vm]
m_e = 9.1093835611e-31;      % electron rest mass [kg]
q_e = 1.602176620898e-19;    % elementary charge [C]



w_p = (max(ne,0)*q_e^2/Ep_0/m_e).^.5;
