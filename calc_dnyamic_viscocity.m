
function mu = calc_dnyamic_viscocity(T)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% AUTHORS:
%
%   Christopher Cox (NOAA) christopher.j.cox@noaa.gov
%
% REFERENCES:
%
%   Standard application of Sutherland's formula.
%
% PURPOSE:
% 
%   Calculate change in dynamic viscosity of air f(T)
%
% INPUT:
%
%   T = temperature in K; may be a vector
%
% OUTPUT:
%
%   mu = dynamic viscosity air in kg m^-1 s^-1
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sutherland's formula for dynanmic viscosity of air

T0  = 273.15   ; % reference temperature
S   = 110.4    ; % Sutherland's temperature
mu0 = 1.716e-5 ; % reference viscocity (corresponds to T0)

% Calculate
mu = mu0.*(T./T0).^(3./2).*((T0+S)./(T+S));

