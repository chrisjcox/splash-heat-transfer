
function Pr = calc_prandtl(cp,mu,kair)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% AUTHORS:
%
%   Christopher Cox (NOAA) christopher.j.cox@noaa.gov
%
% REFERENCES:
%
%   Definition of Prandtl number. See wikipedia.
%
% PURPOSE:
% 
%   Calculate Pradntl number.   
%
% INPUT:
%
%   (may be vectors)
%   cp = specfic heat of air [m^2 s^-2 K^-1]
%   mu = dynamic viscosity of air [kg m^-1 s^-1]
%   kair = thermal condictivity of air [W m^-1 C^-1]
%
% OUTPUT:
%
%   Pr = Prandtl number [dimensionless]
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Pr = (cp .* mu) ./ kair;