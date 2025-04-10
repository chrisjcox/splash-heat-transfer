
function [T,x,ustar,alpha_c] = diff1d(Tinf,Twall,tstar,xht,dx,A,rho,F,Uref,sources,ustar,alpha_c)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% AUTHORS:
%
%   Christopher Cox (NOAA) christopher.j.cox@noaa.gov
%
%   Appreciate input from Chris Fairall (NWRA) & Lorena Barba's online
%   Navier-Stokes materials.
%
% REFERENCES:
%
%   Barba, Lorena A., and Forsyth, Gilbert F. (2018). J. Open Source Edu., 
%      1, 21, https://doi.org/10.21105/jose.00021
%
%   Dutsch et al. (submitted, 2025). J. Adv. Mod. Earth Sys.  
% 
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
% where Schmidt number scaling is replaced with the Prandtl number to
% yield diffusivity:
%
% ∂T/∂t = (𝛎/Pr*∂^2T/∂x^2  + Km/Prt*∂^2T/∂x^2 + ɑ_c*ustar*xstar*exp(-Az))
%
% The Reynold's Analogy is applied so Stu = Prt = 1.
%
% The problem is framed as the development of a thermal boundary layer for
% a turbulent fluid of one temperature adjacent to a plate of another
% temperature. Wall temperture and bulk fluid temperature cannot change.
% The plate need not be flat (case c), but currently skin drag is not
% directly parameterized.
%
% INPUT:
%   Required:
% tstar   = time in seconds to run the model
% xht     = depth of fluid. Practical limit [m]. The physics assume inf.
% dx      = discretized incriment of depth [m]
% A       = height of form drag surface features (e.g., wave height)
% rho     = air density
% F       = surface sensible heat flux (W/m^2)
% Uref    = wind speed. Should be at xht.
% sources: 1 = visc, 2 = visc+turb, 3 = visc+turb+form
%   Optional (if ignore, use []): 
% ustar   = friction velocity [m/s]
% alpha_c = fraction of momentum

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN:
%tic
fprintf('\nPreliminaries...\n')

% % % Some constants
kappa = 0.4                                                ; % von Karman
Prt = 0.85                                                 ; % turbulent Prandtl number from based on RNG method from Li (2019) (https://doi.org/10.1016/j.atmosres.2018.09.015)
cp = 1005                                                  ; % specific heat, air, J/(kg K)
m = 1                                                      ; % Fairall et al. (2000)
lmda = 12                                                  ; % Reichardt (1951) 
Us = 0                                                     ; % define wind speed at the surface as 0 m/s

% % % Frame the problem
nx = xht/dx                                                ; % number of spatial layers
x = (1:nx)*dx                                              ; % [m]
Ainv = 1/A                                                 ; % We actually need the inverse

% % % Initialize
T = ones(1, int32(nx))*Tinf                                ; % Bulk fluid temperature
TK = Tinf+273.15                                           ; % in Kelvins
T(1) = Twall                                               ; % Wall temperature

% % % Preliminary calculations
if isempty(ustar)
    ustar = edson2013_ustar(Uref)                          ; % Calulate ustar if it wasn't passed
end
fprintf(['ustar = ',num2str(ustar),'\n'])
xstar = -F/cp/(ustar*rho)                                  ; % MO scaling parameter derived Dutsch et al,  Eq. (16)
kair = calc_kair(Tinf)                                     ; % Thermal Conductivity, air. Kannuluik & Carman (1951) (https://doi.org/10.1071/CH9510305)
mu = calc_dnyamic_viscocity(TK)                            ; % Sutherland's formula for dynanmic viscocity of air
nu = mu/rho                                                ; % m2/s kinematic viscocity, air
Pr = calc_prandtl(cp,mu,kair)                              ; % Prandtl Number 
delt = lmda*nu/ustar                                       ; % dissipation length scale

% % % Solve for alpha_c

if isempty(alpha_c)
    alpha_c = alpha_solver(kappa,nu,ustar,xht,Ainv,delt,m,Uref,Us);
end

fprintf(['alpha_c = ',num2str(alpha_c),'\n'])
alpha_h = alpha_c * 0.3; % Mueller and Veron (2010)

% % % Recalc Km using optimal alpha_c

if sources < 3
    Km = kappa*ustar*x;
else
    Km = Kmcalc(kappa,ustar,x,alpha_c,Ainv,delt,m);
end


% % % Define dt for stability

% 1) For this equation, diffusivity is our analog to fluid velocity
%    We will use the max value.
Dvisc = nu/Pr;
Dturb = max(Km);
Dform = max(calc_Dform(alpha_h,ustar,xstar,Ainv,x));
Dmax  = Dvisc+Dturb+Dform;

% 2) Sigma is Courant–Friedrichs–Lewy number. 0.2-0.5 are typical.
%    This is probably unconventional, but to make the code run more
%    efficiently, I will choose the largest value of dt that satisfies the
%    CFL criterion that sigma <= 1
sigma = 1e1; % arbitrary number > 1
dt = 1e1; % 10s is too coarse
while sigma > 1
    dt = dt/1.01;
    sigma = Dmax*2*dt/dx^2;
end
fprintf(['sigma = ',num2str(sigma),'\n'])
fprintf(['dt = ',num2str(dt),'\n'])

% 3) Calculate number of time steps
nt = tstar/dt;


% % % Solve the 1d diffusion equation
fprintf('Solving 1D diffusion...\n')

for n = 2:nt  % for time

    Tn = T; % temporary copy
    
    for i = 2:nx-1 % for height (distance from wall)

        viscous_term   = (nu./Pr .* dt / dx^2 * (Tn(i+1) - 2 * Tn(i) + Tn(i-1)));
        turbulent_term = (Km(i)./Prt .* dt / dx^2 * (Tn(i+1) - 2 * Tn(i) + Tn(i-1)));
        formdrag_term  = calc_Dform(alpha_h,ustar,xstar,Ainv,x(i));

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

end
fprintf(['Simulation run to t* = ',num2str(tstar),'s.\nDone.\n\n'])
%toc
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBS:


% Calculate the form drag term

function Dform = calc_Dform(alpha_h,ustar,xstar,Ainv,x)

    Dform = (alpha_h*ustar*xstar.*exp(-(Ainv)*x));

end

% Calculate ustar from Uref using Edson et al. (2013) parameterization
% https://doi.org/10.1175/JPO-D-12-0173.1

function ustar = edson2013_ustar(Uref)

    ustar = sqrt((1.03e-3 + 0.04e-3*Uref^1.48) / Uref^0.21)*Uref;

end


% Calculate Km
% Eq. (10) from Dutsch et al.
% Will return a alph by x matrix if alph is a vector 

function Km = Kmcalc(kappa,ustar,x,alph,Ainv,delt,m)

    Km = (kappa*ustar.*x .* sqrt(1-alph'.*exp(-Ainv.*x))) ./ (1+(delt./x).^m); % Eq 10

end


% Solve for alpha_c using Dutsch et al. Eqs 10, 12, 13

function alpha_c = alpha_solver(kappa,nu,ustar,xht,Ainv,delt,m,Uref,Us)

    disp('Solving for alpha_c...')

    % first, specify x at high resolution
    dx = 1e-6;
    nx = xht/dx;
    x = (1:nx)*dx;

    % calculate for a coarse range of alpha to save time
    alpha_vector = 0:0.25:1;
    Km = Kmcalc(kappa,ustar,x,alpha_vector,Ainv,delt,m);
    dUdx = (ustar^2.*(1-alpha_vector'.*exp(-Ainv.*x))) ./ (nu + Km); % Eq 12
    dU = (ustar*sum(dUdx,2)*dx);

    % interp alpha results
    diff_hi = 1e-5;
    alpha_vector_hi = 0:diff_hi:1;
    dU_hi = interp1(alpha_vector,dU,alpha_vector_hi);

    % find the alpha that best fits Eq 13
    ii = find(abs(dU_hi-(Uref-Us))<diff_hi*2); 

    if length(ii) > 1; ii = ii(1); end
    if isempty(ii) % happens at low wind speed
        alpha_c = 0;
    else
        alpha_c = abs(alpha_vector_hi(ii)); 
    end

end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%