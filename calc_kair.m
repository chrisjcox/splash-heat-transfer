
function kair = calc_kair(T)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% AUTHORS:
%
%   Christopher Cox (NOAA) christopher.j.cox@noaa.gov
%
% REFERENCES:
%
%   Kannuluik & Carman (1951) The Temperature Dependence of the Thermal 
%       Conductivity of Air, Australian J. Sci. Res. 4(3), 305 - 314, 
%       https://doi.org/10.1071/CH9510305
%
% PURPOSE:
% 
%   Calculate thermal conductivity of air f(T)
%
% INPUT:
%
%   T = temperature in degC; may be a vector
%
% OUTPUT:
%
%   kair = thermal conductivity in W m^-1 C^-1
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Kannuluik & Carman polynomial
kair = 5.75e-5.*(1+0.00317.*T-0.0000021.*T.^2) ; 

% Convert from cal^-1 cm^-1 s^-1 C^-1 to W m^-1 C^-1
kair = kair * 4.1868 * 100                     ; 