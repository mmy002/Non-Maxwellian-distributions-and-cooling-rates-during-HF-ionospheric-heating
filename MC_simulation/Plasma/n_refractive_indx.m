function n = n_refractive_indx(w,sigma,epsilon)
% N_REFRACTIVE_INDX - Complex refractive index for HF-waves in electron plasma
%   
% Calling: 
%  n = n_refractive_indx(w,sigma,epsilon)
% Input:
%  w       - angular frequency (rad/s)
%  sigma   - conductivity tensor
%  epsilon - dielectric tensor

% (c) 2008 B Gustavsson, all rights reserved
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.

% Ep_0    = 1e-9/36/pi;               % Permittivity [As/Vm]
c0  = 2.99792458e8;                   % Speed of light [m/s]
my_0    = 4*pi*1e-7;                  % Permeability [Vs/Am]
Ep_0    = 1/(my_0*c0^2);              % Permittivity [As/Vm]

if all(size(sigma) == 3)
  n{1} = sqrt(epsilon{1,1}+1i*(1./(w*Ep_0)).*sigma{1,1});
  n{2} = sqrt(epsilon{2,2}+1i*(1./(w*Ep_0)).*sigma{2,2});
  n{3} = sqrt(epsilon{3,3}+1i*(1./(w*Ep_0)).*sigma{3,3});
else
  n = sqrt(epsilon+1i*(4*pi./w).*sigma);
end
