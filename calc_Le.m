
function [Le, Lei] = calc_Le(T)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% AUTHORS:
%
%   Christopher Cox (NOAA) christopher.j.cox@noaa.gov
%
% REFERENCES:
%
%   Note: The equation is widely used in the atmospheric community and
%   almost never referenced, except occasionally to other non-primary
%   sources. It's worth detailing the history here. The most easily 
%   obtainable origin is an aside in Harrison (1965) that first notes the 
%   equation of Ferrel (1886):
%
%       Le = 606.5 - 0.695*T cal/g
%          = (2.5393 - 0.0029*T)*1e6 J/kg
%
%   Harrison goes on to state that Goff and Gratch have updated the
%   formula, referencing both List (1958: here 1968) and ASHRAE (1963). The
%   update equation is given by Harrison as
%
%       Le = 597.31*(1 - 0.5637/597.31*T) ITcal/g
%          = 2.501*(1 - 0.00236/2.501*T)*1e6 J/kg
%
%   and the modern formulation is
%
%          = 2.501 - 0.00237*T)*1e6 J/kg
%
%   Errors in the conversion are trivial, ~0.01% and Ferrel vs Harrison is
%   about ~ -0.21%/C and +1.4% at 0 C.
%
%   ASHRAE (1963) ASHRAE Guide and Data Book 1963: Fundamentals and 
%       Equipment. Ch 3, "Psychometry and Air-Conditioning Theory".
%
%   Ferrel, W. (1886) Report on psychorometic tables for use in the signal
%       service. Annual Report of the Chief Signal Officer, Appendix 24,
%       pp. 233-259.
%
%   Harrison, L.P. (1965) Some fundamental considerations regarding 
%       psychrometry. In Humidity and Moisture: Measurement and Control in 
%       Science and Industry. A. Wexler, Ed., Reinhold, New York.
%       - refer to page 79
%
%   List, R.J. (1968) Smithsonian Meteorological Tables, 6th ed. 
%       Smithsonian Institute Press, DC. 
%       - refer to page 343 in the 6th edition
%
% PURPOSE:
% 
%   Calculate the latent heat of vaporization
%
% INPUT:
%
%   T = temperature in C
%
% OUTPUT:
%
%   Le = latent heat of vaporization [J kg^-1]
%       - use for evaporation of liquid water
%   Lei = latent heat of vaporization + latent heat of fusion [J kg^-1]
%       - use for sublimation of ice
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For evaporation...
Le=(2.501-.00237*T)*1e6;

% For sublimation additional energy is required to support two phase 
% changes. Note that latent heat of fusion (LeF) is not considered
% strongly dependent on temperature. 
LeF = 333.55*1e3; % [J kg^-1]
Lei = Le + LeF; 
