
function Re = calc_reynolds_number(rho,u,L,mu)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% AUTHORS:
%
%   Christopher Cox (NOAA) christopher.j.cox@noaa.gov
%
% REFERENCES:
%
%   Definition of Reynold's number. See wikipedia.
%
% PURPOSE:
% 
%   Calculate Reynold's number.   
%
% INPUT:
%
%   (may be vectors)
%   rho = air density [kg m^-3]
%   u = wind velocity [m s^-1]
%   L = characteristic length [m]
%   mu = dynamic viscosity of air [kg m^-1 s^-1]
%
% OUTPUT:
%
%   Re = Reynolds number [dimensionless]
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Re = (rho.*u.*L) ./ mu;