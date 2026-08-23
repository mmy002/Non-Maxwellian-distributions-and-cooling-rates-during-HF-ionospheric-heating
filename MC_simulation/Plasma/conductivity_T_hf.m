function sigma = conductivity_T_hf(w,ne,ny_e,B)
% CONDUCTIVITY_T_HF - High frequency electron conductivity tensor
% 
% Calling:
%  ma = conductivity_T_hf(w,ne,ny_e,B)
% Input:
%  w - angular frequency (radian/s)
%  ne - electron density (m^-3)
%  ny_e - electron collision frequncy (/s)
%  B - Magnetic field (T)
% Output:
%  sigma - conductivity-tensor, returned as a [3 x 3] cell-array with the
%          conductivities || conductivity-profile in sigma(3,3), and the
%          left and right-hand conductivities in the (1,1) and (2,2)
%          elements.
% 
% After Shoucri 1984

% (copyright) 2008 B Gustavsson, all rights reserved
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.

% q_e     = 1.6021773e-19;            % elementary charge [C]
% m_e     = 9.10939e-31;              % electron rest mass [kg]
% Ep_0    = 1e-9/36/pi;               % Permittivity [As/Vm]

m_e   = 9.1093835611e-31;      % electron rest mass [kg]
q_e   = 1.602176620898e-19;    % elementary charge [C]
c0    = 2.99792458e8;          % Speed of light [m/s]
my_0  = 4*pi*1e-7;             % Permeability [Vs/Am]
Ep_0  = 1/(my_0*c0^2);         % Permittivity [As/Vm]

w_p = (ne*q_e^2/Ep_0/m_e).^.5; % plasma frequency (angular)
w_c = q_e*B/m_e;               % electron gyro frequency (angular)
w_c = repmat(w_c,1,size(ne,2));


sigma{3,3} = Ep_0*ny_e.*w_p.^2./(w.^2+ny_e.^2);
sigma{2,2} = Ep_0*ny_e.*w_p.^2./((w-w_c).^2+ny_e.^2);
sigma{1,1} = Ep_0*ny_e.*w_p.^2./((w+w_c).^2+ny_e.^2);
