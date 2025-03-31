
function Nu = calc_nusselt(Pr,Re)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% AUTHORS:
%
%   Christopher Cox (NOAA) christopher.j.cox@noaa.gov
%
% REFERENCES:
%
%   Lienhard, J.H. (2020) Heat transfer in flat-plate boundary layers: A 
%       correlation for laminar, transitional, and turbulent flow. J.
%       Heat Transfer, 142, 061805. https://doi.org/10.1115/1.4046795
%
% PURPOSE:
% 
%   Calculate the Nusselt number using Lienhard's 2020 parameterization.
%   Generally, Nu = C * Re^m * Pr^n.
%   This is for Nusselt >= 0.6.
%
% INPUT:
%
%   (may be vectors)
%   Pr = Prandtl number (>= 0.6, else return error)
%   Re = Reynold's number
%
% OUTPUT:
%
%   Nu = Nusselt number [dimensionless]
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Test that the flow qualifies for this equation
i_prlow = find(Pr < 0.6);
if length(i_prlow) == length(Pr)
    disp('Prandtl number too small to calculate Nusselt number. Aborting.');
    return
elseif ~isempty(i_prlow) & ~isscalar(Pr)
    disp([num2str(length(i_prlow) / length(Pr) * 100),'% of Prandtl numbers too small to calculate Nusselt number. NaN will be returned at these indices.'])
    Pr(i_prlow) = NaN;
end

% Calculation
Cf = 0.455 ./ (log(0.06.*Re)).^2; % Eq. (7)

Nu = (Re.*Pr.*(Cf./2)) ./ (1+12.7.*(Pr.^(2/3)-1).*sqrt(Cf./2)); % Eq. (6)