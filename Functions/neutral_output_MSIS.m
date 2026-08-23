function [rho,temperature] = neutral_output_MSIS(height)
% This function process data from MATLAB MSIS
% INPUT:
% height: Altitude in [m]
% OUTPUT: 
% rho: Densities of neutral in [m^{-3}]
% temperature: Neutral temperature in [K]
  
scale_10p7flux = 1;       % Skalerer F10,7-indeksen (Mengden stråling med bølgelengde 10,7 cm)
scale_ap = 1;             % Skalerer ap-indeksen (Forstyrrelser i den horisontale delen av Jordas magnetfelt)

% Same date and time as Carsten, 8 Sep. 2010 23:55 LT or 12:15 LT, Andenes 69 degree north, 16
% degree east. 
altitude  = height;                         % Heigth over sea level [m]
latitude  = 69.0;                           % Self explanatory 
longitude = 16.0;                           % Self explanatory
year      = 2010;                           % Self explanatory
dayOfYear = 31+28+31+30+31+30+31+31+8;      % 31 days in January + 28 in February osv. for 8. sep
time = 23+55/60;                            % Local time 23:55 
UTC = 1;                                    % Time zone: UTC +1
UTseconds = (time-UTC)*3600;                % Time in universal time [s]

% Geophysical indices taken from MSIS-files 
f107Average = 126.1*scale_10p7flux;                 % F10.7 3-months average
f107Daily   = 111.8*scale_10p7flux;                 % F10.7_daily
magneticIndex = [22.4000000000000, 15, 39, 15, 39, ...   % Magnetic indices 
		 25.9000000000000, 46.6000000000000]*scale_ap;

% Getting data from Matlab MSIS atmosphere model
[temperature,rho] = atmosnrlmsise00(altitude, ...
				    repmat(latitude,     length(altitude),1), ...
				    repmat(longitude,    length(altitude),1), ...
				    repmat(year,         length(altitude),1), ...
				    repmat(dayOfYear,    length(altitude),1), ...
				    repmat(UTseconds,    length(altitude),1), ...
				    repmat(f107Average,  length(altitude),1), ...
				    repmat(f107Daily,    length(altitude),1), ...
				    repmat(magneticIndex,length(altitude),1));



