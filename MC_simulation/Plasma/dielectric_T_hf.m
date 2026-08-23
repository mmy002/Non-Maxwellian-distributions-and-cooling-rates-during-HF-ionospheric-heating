function epsilon = dielectric_T_hf(w,ne,ny_e,B)
% DIELECTRIC_T_HF - plasma dielectric tensor, high frequency case
%   
% Calling:
%  epsilon = dielectric_T_hf(w,ne,ny_e,B)
% Input:
%  w - angular frequency (radian/s)
%  ne - electron density (m^-3)
%  ny_e - electron collision frequncy (/s)
%  B - Magnetic field
% Output:
%  epsilon - dielectric tensor for EM-waves in magnetized plasma
%          B-parallel dielectric-profile in epsilon(3,3), and the
%          left and right-hand dielectric-profile in the (1,1) and (2,2)
%          elements.
% 
% From Shoucri 1984

% (c) 2008 B Gustavsson, all rights reserved
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.

% Old values for physical constants
% q_e     = 1.6021773e-19;            % elementary charge [C]
% m_e     = 9.10939e-31;              % electron rest mass [kg]
% Ep_0    = 1e-9/36/pi;               % Permittivity [As/Vm]
% New values for physical constants
m_e   = 9.1093835611e-31;      % electron rest mass [kg]
q_e   = 1.602176620898e-19;    % elementary charge [C]
c0    = 2.99792458e8;          % Speed of light [m/s]
my_0  = 4*pi*1e-7;             % Permeability [Vs/Am]
Ep_0  = 1/(my_0*c0^2);         % Permittivity [As/Vm]


w_p = (ne*q_e^2/Ep_0/m_e).^.5; % plasma frequency (angular)
w_c = q_e*B/m_e;               % electron gyro frequency (angular)
w_c = repmat(w_c,1,size(ne,2));

epsilon{3,3} = 1 - w_p.^2./(w.^2+ny_e.^2);
epsilon{2,2} = 1 - (w-w_c).*w_p.^2./(w.*((w-w_c).^2+ny_e.^2));
epsilon{1,1} = 1 - (w+w_c).*w_p.^2./(w.*((w+w_c).^2+ny_e.^2));
