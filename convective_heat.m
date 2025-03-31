
function [Q,h] = convective_heat(T,Tsfc,L,Nu,kair,A)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% AUTHORS:
%
%   Christopher Cox (NOAA) christopher.j.cox@noaa.gov
%
% REFERENCES:
%
%   n/a
%
% PURPOSE:
% 
%   Calculates the convective heat transfer coefficient using the equation
%   for the Nusselt number and then applies it to the basic heat transfer
%   equation. Unconventional, but it permits h to be calculated
%   independently of Q for comparison to direct observations.
%
% INPUT:
%
%   (may be vectors)
%   T = bulk fluid temperature [K]
%   Tsfc = surface (wall) temperature [K]
%   L = characteristic length [m]
%   Nu = Nusselt number
%   kair = thermal conductivity of fluid [W m^-1 C^-1] 
%   A = area [m^2], usually 1 if Q is W m^-2
%
% OUTPUT:
%
%   h = convective heat transfer coefficient [W m^-2 K^-1]
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Back out the convective heat transfer coeficient, h, from Nu 
h = (Nu .* kair) ./ L;
% Convective heat transfer equation, Q in W/m2
Q = h .* A .* (Tsfc-T);