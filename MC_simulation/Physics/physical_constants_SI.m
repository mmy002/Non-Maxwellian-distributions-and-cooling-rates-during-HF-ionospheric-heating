%---------------------------------------------------------------
% physical_constants
%
% This scriptfile sets the fysical constants needed for calculations.
% The constants are set as global variables.
%
% Global variabels:
% c0	= 2.99792458e8;		% Speed of light [m/s]
% h	= 6.62618e-34;		% Plank's constant [Js]
% kB	= 1.380662e-23;		% Boltzmann constant [J/K]
% Na	= 6.022169e26;		% Avogadros number [molecules/kmol]
% R 	= kB * Na;		% Molar gas constant [J/kmol K]
% Re	= 6.378e6;		% Radius of earth [m]
% cbgr	= 2.735;		% Cosmic background radiation [K]
%---------------------------------------------------------------

%-------------------------------------
% Copy these declarations into the file where the function where the 
% constants are needed.

c0       = 2.99792458e8;              % Speed of light [m/s]
h        = 6.62607015e-34;            % Plank's constant [Js]
kB       = 1.380649e-23;              % Boltzmann constant [J/K]
Na       = 6.022169e26;               % Avogadros number [molecules/kmol]
R        = kB * Na;                   % Molar gas constant [J/kmol K]
Re       = 6.378e6;                   % Radius of earth [m]
cbgr     = 2.735;                     % Cosmic background radiation [K]
my_0     = 4*pi*1e-7;                 % Permeability [Vs/Am]
Ep_0     = 1/(my_0*c0^2);             % Permittivity [As/Vm]
G        = 6.6726e-11;                % Gravitational constant [Nm^2/kg^2] 
% m_e    = 9.10938291e-31;            % electron rest mass [kg] pre 20180530
m_e      = 9.1093835611e-31;          % electron rest mass [kg]
m_p      = 1.672621778e-27;           % proton rest mass [kg]
m_n      = 1.674927352e-27;           % neutron rest mass [kg]
% q_e    = 1.602176565e-19;           % elementary charge [C] pre 20180530
q_e      = 1.602176620898e-19;        % elementary charge [C]
a_0      = h^2*Ep_0/(pi*m_e*q_e);     % Bohr radius [m]
my_B     = q_e*h/(2*pi*2*m_e);        % Bohr magneton [Am^2]
sigma_SB = 5.67037441918442945397e-8; % W/m^2/K^4

