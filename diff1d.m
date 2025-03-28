
function [T,x,ustar,alpha_c] = diff1d(Tinf,Twall,tstar,xht,dx,A,rho,F,Uref,sources,ustar,alpha_c)

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
%   Required:
% tstar   = time in seconds to run the model
% xht     = depth of fluid. Practical limit [m]. The physics assume infinite.
% dx      = discretized incriment of depth [m]
% A       = height of form drag surface features (e.g., wave height)
% rho     = air density
% F       = surface sensible heat flux (W/m^2)
% Uref    = wind speed
% sources: 1 = visc, 2 = visc+turb, 3 = visc+turb+form
%   Optional (if ignore, use []): 
% ustar   = friction velocity [m/s]
% alpha_c = fraction of momentum

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN:
disp(' ')
disp('Preliminaries...')

% Some constants
nu = 1.33e-5                       ; % m2/s kinematic viscocity, air
kappa = 0.4                       ; % von karman
D = 22e-6                         ; % mass diffusivity of air
Sx = nu/D                         ; % Schmidt number 
St = 1                            ; % turbulent Schmidt number: 1 makes Reynolds analogy, but it could be larger
cp = 1005                         ; % specific heat, air
m = 1                             ; % Fairall et al. (2000)
lmda = 12                         ; % Reichardt (1951) 
Us = 0                            ; % define wind speed at the surface as 0 m/s

% Frame the problem
dt = 0.2 * dx^2                   ; % define dt for stable solution wrt dx
nx = xht/dx                       ;
x = (1:nx)*dx                     ; % [m]
Ainv = 1/A                        ; % We actually need the inverse
nt = 1e6                          ; % Specify large number of steps. Code will run until we hit tstar

% Initialize
T = ones(1, int32(nx))*Tinf       ; % Bulk fluid temperature
T(1) = Twall                      ; % Wall temperature

% Preliminary calculations
if isempty(ustar)
    ustar = edson2013_ustar(Uref) ; % Calulate ustar if it wasn't passed
end
disp(['ustar = ',num2str(ustar)])
xstar = -F/cp/(ustar*rho)         ; % MO scaling parameter derived Dutsch et al,  Eq. (16)
delt = lmda*nu/ustar              ; % dissipation length scale

% Solve for alpha_c
if isempty(alpha_c)
    alpha_c = alpha_solver(kappa,nu,ustar,xht,Ainv,delt,m,Uref,Us);
end
disp(['alpha_c = ',num2str(alpha_c)])
alpha_h = alpha_c * 0.3; % Mueller and Veron (2010)

% Recalc Km using optimal alpha_c
if sources < 3
    Km = kappa*ustar*x;
else
    Km = Kmcalc(kappa,ustar,x,alpha_c,Ainv,delt,m);
end

% Solve the 1d diffusion equation
disp('Solving 1D diffusion...')
for n = 2:nt  % for time

    Tn = T; % temporary copy
    
    for i = 2:nx-1 % for height (distance from wall)

        viscous_term   = (nu./Sx .* dt / dx^2 * (Tn(i+1) - 2 * Tn(i) + Tn(i-1)));
        turbulent_term = (Km(i)./St .* dt / dx^2 * (Tn(i+1) - 2 * Tn(i) + Tn(i-1)));
        formdrag_term  = (alpha_h*ustar*xstar.*exp(-(Ainv)*x(i)));

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
        disp(['Simulation run to t* = ',num2str(tstar),'s.'])
        disp('Done.')
        disp(' ')
        break; 
    elseif n == nt
        disp('Need more time.')
        disp('Done.')
        disp(' ')
    end

end

end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBS:


% Calculate ustar from Uref using Edson et al. (2013) parameterization
% https://doi.org/10.1175/JPO-D-12-0173.1

function ustar = edson2013_ustar(Uref)

    ustar = sqrt((1.03e-3 + 0.04e-3*Uref^1.48) / Uref^0.21)*Uref;

end


% Calculate Km
% Eq. (10) from Dutsch et al.

function Km = Kmcalc(kappa,ustar,x,alph,Ainv,delt,m)

    Km = (kappa*ustar.*x .* sqrt(1-alph.*exp(-Ainv.*x))) ./ (1+(delt./x).^m); % Eq 10

end


% Solve for alpha_c using Dutsch et al. Eqs 10, 12, 13

function alpha_c = alpha_solver(kappa,nu,ustar,xht,Ainv,delt,m,Uref,Us)

    disp('Solving for alpha_c...')

    % first, specify x at high resolution
    dx = 1e-6;
    nx = xht/dx;
    x = (1:nx)*dx;

    % a little slow...
    % we will be coarse with the alpha sampling and interpolate later
    % I'll look for a more robust way to speed this up later
    alpha_vector = 0:0.1:1;
    dU = NaN*alpha_vector;
    for k = 1:length(alpha_vector)
        Km = Kmcalc(kappa,ustar,x,alpha_vector(k),Ainv,delt,m);
        dUdx = (ustar^2.*(1-alpha_vector(k).*exp(-Ainv.*x))) ./ (nu + Km); % Eq 12
        dU(k) = (ustar*trapz(x,dUdx)); % Eq 13
    end

    % interp alpha results
    alpha_vector_hi = 0:0.001:1;
    dU_hi = interp1(alpha_vector,dU,alpha_vector_hi);

    % find the alpha that best fits Eq 13
    ii = find(abs(dU_hi-(Uref-Us))<mode(diff(alpha_vector_hi))*10); 
    if length(ii) > 1; ii = ii(1); end
    if isempty(ii) % happens at low wind speed
        alpha_c = 0;
    else
        alpha_c = abs(alpha_vector_hi(ii)); 
    end

end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%