function w_e = w_e_gyro(B)
% W_E_GYRO - angular electron gyro frequency (non-relativistic)
%   
% Calling: 
%  w_e = w_e_gyro(B)
% Input:
%  B magnetic field [T]

m_e = 9.1093835611e-31;      % electron rest mass [kg]
q_e = 1.602176620898e-19;    % elementary charge [C]


w_e = q_e*B/m_e;
