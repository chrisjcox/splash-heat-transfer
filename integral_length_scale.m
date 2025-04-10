
function [L,DN] = integral_length_scale(u,dn,dt,integ_time)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% AUTHORS:
%
%   Christopher Cox (NOAA) christopher.j.cox@noaa.gov
%
% REFERENCES:
%
%   Pope, S.B. (2000) Turbulent Flows, Cambridge University Press
%
%   Schrader, P. (1993) Computing the statistical stability of integral
%       length scale measurements by autoregressive simulation. J. Wind
%       Eng. Indust. Aerodyn., 46-47, 487-496,
%       https://doi.org/10.1016/0167-6105(93)903160-G
%
%   zum Berge, K. (2021) A two-day case study: ... Bound.-Lay. Meteorol.,
%       180, 53-78, https://doi.org/10.1007/s10546-021-00608-2
%
% PURPOSE:
% 
%   Calculate the turbulent integral length scale, which is defined as the
%   integral of the autocorrelation function of fluctuations in the 
%   horizontal velocity (e.g., Pope, 2000). Here, we multiply by the wind
%   speed, following the approach of Schrader (1993) and zum Berge (2021).
%
%   There is ambiguity in the data preparation of the fluctuations and
%   assumes as a default that the autocorrelation function decays to 0 in 
%   less than 10 min. Therefore, this is a practical implementation.
%
% INPUT:
%
%   u = wind velocity in m s^-1
%   dn = time scale in Matlab datenums
%   integ_time = time scale for integration in min
%   dt = sample rate of u n Hz
%
% OUTPUT:
%
%   L = L11 = longitudinal length scale of turbulence
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Smooth with boxcar convolution using integration period (10 min)

U = movmean(u,dt*60*integ_time); % dt [s^-1] * 60 [s] * integ_time  [min]

% Calculate fluctuations about mean wind speed.
r = u - U;

% Make a time vector 
DN = datenum(floor(dn(1)):1/(1440/integ_time):ceil(dn(end)));

% Calculate
L = NaN.*DN;
for h = 1:length(DN)

    % subset the period to integrate
    ind = find(dn >= DN(h) & dn < DN(h)+2/((1440/integ_time))); % the 2 doubles the length to ensure robust autocorr
    if isempty(ind); continue; end
    
    % calculate autocorrelation function
    [rho_a lags] = xcorr(r(ind)',r(ind)',dt*60*integ_time);
    
    % normalize
    rho_a = rho_a(lags>=0)'./max(rho_a);
    
    % get a time vector
    time_vector = (1:length(rho_a))./dt; % [s]
    
    % zum Berge Eq. (5)
    L(h) = abs( nanmean(u(ind)).* trapz(time_vector(1:dt*60*integ_time),rho_a(1:dt*60*integ_time)) ); % [m]

end