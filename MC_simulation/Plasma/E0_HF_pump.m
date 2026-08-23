function [E0,E0b] = E0_HF_pump(epsilon,beta,ERP,h0)
% E0_HF_PUMP - E-field at lowest ionospheric altitude
%   Radial free-space propagation (E \propto 1/h0, W \propto 1/h^2)
%   
% Calling:
%  E0_HF_pump(epsilon,beta,ERP,h0)
% Input:
%  epsilon - dielectric coeficient at lowest altitude
%  beta    - low altitude absorption [0-1]
%  ERP     - effective radiated power (effect*antenna gain) (W)
%  h0      - altitude (m)
% Output:
%  E0 - Pump-wave electrical field at lowest altitude (V/m)

% (c) 2008 B Gustavsson, all rights reserved
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.


E0=0.25*(ERP*beta/1e3)^(1/2)/h0*1e3;
if nargout > 1
  E0b = abs((2/(1+epsilon(1)^(1/2)))*(60*beta*ERP)^(1/2)/h0);
end
