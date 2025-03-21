
function [T,x] = diff1d(Tinf,Twall,tstar,xht,dx,ustar,A,rho,F,Uref,sources)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUTHORS:
%
%   Christopher Cox (NOAA) christopher.j.cox@noaa.gov
%
%   Appreciate input from Chris Fairall (NWRA)
%
% PURPOSE:
%
% Solves 1d heat diffusion equation,
%
% ∂T/∂t = c*∂^2T/∂x^2    (1)
% 
% using a forward difference scheme. 
% 
% Cases -
% (case a) c is thermal diffusivity for a laminar solution.
% (case b) c is eddy diffusivity for a turbulent solution.
% (case c) c is the sum of viscous, turbulent, and form drag stresses.
%
% Case c: this code becomes a transient solution for the steady-state model
% proposed by Dutsch et al. (2025). Thus,
%
% ∂T/∂t = (𝛎/Smo*∂^2T/∂x^2  + Km/Stu*∂^2T/∂x^2 + ɑ_c*ustar*xstar*exp(-Az))
%
% The problem is framed as the development of a thermal boundary layer for
% a turbulent fluid of one temperature adjacent to a plate of another
% temperature. Wall temperture and bulk fluid temperature cannot change.
% The plate need not be flat (case c).
%
% INPUT:
%
% tstar   = time in seconds to run the model
% xht     = depth of fluid. Practical limit [m]. The physics assume infinite.
% dx      = discretized incriment of depth [m]
% ustar   = friction velocity [m/s]
% A       = height of form drag surface features (e.g., wave height)
% rho     = air density
% F       = surface sensible heat flux (W/m^2)
% Uref    = wind speed
% sources: 1 = visc, 2 = visc+turb, 3 = visc+turb+form
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Some constants
nu = 1.3e-5                  ; % m2/s kinematic viscocity, air
kappa = 0.4                  ; % von karman
D = 22e-6                    ; % mass diffusivity of air
Sx = nu/D                    ; % Schmidt number 
St = 1                       ; % turbulent Schmidt number: 1 makes Reynolds analogy, but it could be larger
cp = 1005                    ; % specific heat, air
m = 1                        ; % Fairall et al. (2000)
lmda = 12                    ; % Reichardt (1951) 
Us = 0                       ; % define wind speed at the surface as 0 m/s
%ustar = sqrt((1.03e-3 + 0.04e-3*Uref^1.48) / Uref^0.21)*Uref;

% Frame the problem
dt = 0.2 * dx^2              ; % define dt for stable solution wrt dx
nx = xht/dx;
x = (1:nx)*dx                ; % [m]
A = 1/A                      ; % We actually need the inverse
nt = 1e6                     ; % Specify large number of steps. Code will run until we hit tstar

% Initialize
T = ones(1, nx)*Tinf         ; % Bulk fluid temperature
T(1) = Twall                 ; % Wall temperature

% Preliminary calculations
xstar = -F/cp/(ustar*rho); % MO scaling parameter derived Dutsch et al,  Eq. (16)
delt = lmda*nu/ustar         ; % dissipation length scale

% Solve for alpha_c using Dutsch et al. Eqs 10, 12, 13 then calculate Km
alpha_vector = -1:0.001:1;
dU = NaN*alpha_vector;
for k = 1:length(alpha_vector)
    Km = (kappa*ustar.*x .* sqrt(1-alpha_vector(k).*exp(-A.*x))) ./ (1+(delt./x).^m); % Eq 10
    dUdx = (ustar^2.*(1-alpha_vector(k).*exp(-A.*x))) ./ (nu + Km); % Eq 12
    dU(k) = (ustar*trapz(x,dUdx)); % Eq 13
    %dU(k) = trapz(x,dUdx);
end
% find the alpha that best fits Eq 13
ii = find(abs(dU-(Uref-Us))<mode(diff(alpha_vector))*2); 
if length(ii) > 1; ii = ii(1); end
if isempty(ii) % happens at low wind speed
    alpha_c = 0;
else
    alpha_c = abs(alpha_vector(ii)); 
end
alpha_c = 0.11
% Recalc Km using optimal alpha_c
Km = (kappa*ustar.*x .* sqrt(1-alpha_c.*exp(-A.*x))) ./ (1+(delt./x).^m); 
alpha_h = alpha_c * 0.3; % Mueller and Veron (2010)

disp('Solving..')
for n = 2:nt  % for time

    Tn = T; % temporary copy
    
    for i = 2:nx-1 % for height (distance from wall)

        viscous_term   = (nu./Sx .* dt / dx^2 * (Tn(i+1) - 2 * Tn(i) + Tn(i-1)));
        turbulent_term = (Km(i)./St .* dt / dx^2 * (Tn(i+1) - 2 * Tn(i) + Tn(i-1)));
        formdrag_term  = (alpha_h*ustar*xstar.*exp(-(A)*x(i)));

        % What components are we including?
        switch sources
            case 1
                T(i) = Tn(i) + viscous_term;
            case 2
                T(i) = Tn(i) + viscous_term + turbulent_term;
            case 3
                T(i) = Tn(i) + viscous_term + turbulent_term + formdrag_term;
        end

    end

    if n*dt >= tstar && n < nt
        disp(['Simulation run to ',num2str(tstar),'s. Done.'])
        break; 
    elseif n == nt
        disp('Need more time. Done.')
    end

end