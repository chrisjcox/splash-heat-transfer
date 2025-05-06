
%% (0) Preliminaries

% Add paths, load Avery Picnic data, 
load('/Users/ccox/Documents/Projects/2021/splash/_level2/a50dbl.mat')
a50 = data10;
a50.temp(a50.temp_qc ~= 0) = NaN;
a50.skin_temp_surface(a50.skin_temp_surface_qc ~= 0) = NaN;
a50.atmos_pressure(a50.atmos_pressure_qc ~= 0) = NaN;
a50.wspd_vec_mean(a50.wspd_vec_mean_qc ~= 0) = NaN;
a50snow_depth = a50.snow_depth;
a50snow_depth(a50.snow_depth_qc ~= 0) = NaN;
a50snow_depth = a50snow_depth.*0.01;
a50snow_depth = naninterp(a50snow_depth);

% Case study is April (14)18-19(22), 2023, in particular when Tair > 0
% Define the case study 
time1 = datenum([2023 4 18 23 19 0]);
time2 = datenum([2023 4 18 23 39 0]);

% Flag warm and cold regimes in general: must be snowy
ind_wrmreg = find(a50.temp > 0 & a50.snow_flag == 1 & a50snow_depth > 0.2);
ind_cldreg = find(a50.temp < 0 & a50.snow_flag == 1 & a50snow_depth > 0.2);

% Define unniversal constants/variables
K0 = 273.15;
cp = 1005                                                            ; % specific heat, air [m^2 s^-2 K^-1]
kair = calc_kair(a50.temp)                                           ; % thermal conductivity of air [W m^-1 C^-1]
mu = calc_dnyamic_viscocity(a50.temp+K0)                             ; % dynamic viscosity of air [kg m^-1 s^-1]
[rho_dry, rho] = calc_rho_air(a50.temp+K0,a50.rh,a50.atmos_pressure) ; % air density [kg m^-3]
eta = 2.178e-5                                                       ; % mass diffusivity of water vapor in air [m^2 s^-1] (Massman, 1998; https://doi.org/10.1016/S1352-2310(97)00391-9)
nu = 1.6e-5                                                          ; % m2/s kinematic viscocity, air


%% (1) ---- First calculate the integral length scale ----
%   Takes a long time, so I saved it.

%get_integralLscale
load /Users/ccox/Documents/Projects/2021/splash/science/sublimation/sos_case_study/intLens
a50.L = interp1_noerror_incNaN(intLens.DN,intLens.L,a50.dn);

%% (2) ---- Convective heat transfer equation ----

% Reynolds Number
Re = calc_reynolds_number(rho,a50.wspd_vec_mean,a50.L,mu);
% Prandtl Number
Pr = calc_prandtl(cp,mu,kair);
% Nusselt Number for boundary layer turbulent flow
Nu = calc_nusselt(Pr,Re);
% Convective heat transfer equation and coefficient
[Q,h] = convective_heat(a50.temp+K0,a50.skin_temp_surface+K0,a50.L,Nu,kair,1);
% Load some precalculated bulk values
load /Users/ccox/Documents/Projects/2021/splash/science/sublimation/sos_case_study/offsettest.mat

% Make a plot for the case study
ind = find(a50.dn >= datenum([2023 4 18 13 0 0]) & a50.dn < datenum([2023 4 19 13 0 0]));
clf;
vals = [hs_50(ind) hs_50p2(ind) hs_50p4(ind) hs_50p10(ind)];
[fillhandle3,msg]=jbfill(a50.dn(ind)',max(vals'),min(vals'),'k','k',[],0.25); hold on;
hold on;
tmpHs = a50.Hs(ind); tmpHs(tmpHs > 20) = NaN;
plot(a50.dn(ind),tmpHs,'linewidth',2); hold on;
tmpQ = Q; tmpQ(abs(tmpQ) > 100) = NaN;
plot(a50.dn(ind),tmpQ(ind),'linewidth',2);
grid on; hline(0);
ylim([-100 25]);
xlim([a50.dn(ind(1)) a50.dn(ind(end))]);
legend('Bulk range','H_s from EC','Q','location','southwest');
datetick('x','mmmdd HH','keepticks','keeplimits');
set(gca,'fontsize',16); ylabel('[W m^-^2]');
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0 0 20 18]);
print('-dpng','-r300','/Users/ccox/Documents/Projects/2021/splash/science/sublimation/pub/figs/Figure_9.png');

rmse_Q = sqrt(nanmean((a50.Hs(ind)-tmpQ(ind)).^2));
rmse_B = sqrt(nanmean((a50.Hs(ind)-nanmean(vals,2)).^2));


%% (3) ---- Surface Renewal Theory ----

% thermal diffusivity of air [m^2 s^-1]
alpha = calc_thermal_diffusivity(kair,rho,cp);
% estimate t* using SRT and h from (2)
tstar = calc_tstar(h,alpha,rho,cp,1);

% Clayson Eq. (8) H
H = calc_srt_flux_claysonEq8(a50.temp,a50.skin_temp_surface,rho,cp,alpha,tstar,'H');

% Clayson Eq. (8) E
es = calc_Pw_Pws_Buck(a50.skin_temp_surface,a50.atmos_pressure,'water',1996); % surface vapor pressure, assume it is saturated there
[Qs,r] = calc_q(es,a50.atmos_pressure); % to specific humidity
q = (a50.mixing_ratio./1000) ./ (1+(a50.mixing_ratio./1000)); % specific humidity
[Le, Lei] = calc_Le(a50.skin_temp_surface); 
E = Le.*rho.*(Qs-q).*(eta./tstar).^0.5;

cmap=lines;
ind3 = find(a50.dn >= datenum([2023 4 18 0 0 0]) & a50.dn <= datenum([2023 4 20 0 0 0]) & a50.snow_flag == 1 & a50.temp > 0);
clf;
h1=plot(H(ind_cldreg),E(ind_cldreg),'.','color',[0.75 0.75 0.75],'markersize',12); hold on;
h2=plot(H(ind_wrmreg),E(ind_wrmreg),'.','color',cmap(1,:),'markersize',12);
h3=plot(H(ind3),E(ind3),'.','color',cmap(3,:),'markersize',12);
grid on; vline(0); hline(0);
plot(-1.*[-125:125],-125:125,'k');
xlim([-125 125]); ylim([-125 125]);
legend([h1 h2 h3],'Cold Regime', ...
    'Warm Regime', ...
    '18-19 April 2023', ...
    'location','southwest');
set(gca,'xtick',[-120:20:120],'ytick',[-120:20:120]);
xlabel('Sensible Heat Flux ({\it H_s}) [Wm^-^2]');
ylabel('Latent Heat Flux ({\it H_l}) [Wm^-^2]');
set(gca,'fontsize',16); 
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0 0 22 20]);
print('-dpng','-r300','/Users/ccox/Documents/Projects/2021/splash/science/sublimation/pub/figs/Figure_12/Figure_12c.png');


ind_wrmregFig12c = find(a50.temp > 0 & a50.snow_flag == 1 & a50snow_depth > 0.2 & abs(H) > 5);
nancorrcoef(H(ind_wrmregFig12c),E(ind_wrmregFig12c))
nanmean(H(ind_wrmregFig12c)+E(ind_wrmregFig12c))
nanmean(H(ind_wrmregFig12c)) / nanmean(E(ind_wrmregFig12c))

%% (4) ---- Test theory at SOS tower, 4/18 case ----

time1 = datenum([2023 4 18 23 24 0]);
time2 = datenum([2023 4 18 23 34 0]);

cd /Users/ccox/Documents/Projects/2021/splash/science/sublimation/sos_case_study/sos_data/

ind50 = find(a50.dn >= datenum([2023 4 18 13 0 0]) & a50.dn < datenum([2023 4 19 13 0 0]));

% Estimate u* for the SOS tower. 
% Found to be similar to ASFS-50 for the case (ASFS-50 u* ~ 0.41)
sosf = rd_netcdf('/Users/ccox/Documents/Projects/2021/splash/science/sublimation/sos_case_study/sos_data/isfs_sos_qc_geo_tiltcor_5min_v2_20230418.nc');
sosf.dn = [datenum([2023 4 18 0 0 0]):5/1440:datenum([2023 4 18 23 59 0])]';
indsos = find(sosf.dn >= time1 & sosf.dn <= time2);

% Get the meteorology
d = dir('isfs_sos_qc_geo_tiltcor_hr_20230418*');
d = [d; dir('isfs_sos_qc_geo_tiltcor_hr_20230419*')];
data.dn = [];
data.t1m = [];
data.t2m = [];
data.t3m = [];
data.t4m = [];
data.t5m = [];
data.t6m = [];
data.t7m = [];
data.t8m = [];
data.t9m = [];
data.t10m = [];
data.t11m = [];
data.t12m = [];
data.t13m = [];
data.t14m = [];
data.t15m = [];
data.t16m = [];
data.t17m = [];
data.t18m = [];
data.t19m = [];
data.t20m = [];
data.pr = [];
data.w10m = [];
data.u = [];
data.w1m = [];
data.w3m = [];
data.tac1m = [];
data.rh10m = [];
data.rh4m = [];
for k = 1:length(d)
    k
    file = rd_netcdf(d(k).name,'base_time','time','T_1m_c','T_2m_c','T_3m_c','T_4m_c','T_5m_c','T_6m_c','T_7m_c','T_8m_c','T_9m_c','T_10m_c', ...
        'T_11m_c','T_12m_c','T_13m_c','T_14m_c','T_15m_c','T_16m_c','T_17m_c','T_18m_c','T_19m_c','T_20m_c', ...
        'P_c_c','u_10m_c','v_10m_c','w_10m_c','u_1m_c','v_1m_c','w_1m_c','u_3m_c','v_3m_c','w_3m_c','Tirga_1m_c','RH_10m_c','RH_4m_c');
    data.dn = [data.dn; datenum([1970 1 1 0 0 0]) + double(file.base_time)/86400 + file.time./86400];
    file.T_1m_c(file.T_1m_c > 10^6) = NaN;
    file.T_2m_c(file.T_2m_c > 10^6) = NaN;
    file.T_3m_c(file.T_3m_c > 10^6) = NaN;
    file.T_4m_c(file.T_4m_c > 10^6) = NaN;
    file.T_5m_c(file.T_5m_c > 10^6) = NaN;
    file.T_6m_c(file.T_6m_c > 10^6) = NaN;
    file.T_7m_c(file.T_7m_c > 10^6) = NaN;
    file.T_8m_c(file.T_8m_c > 10^6) = NaN;
    file.T_9m_c(file.T_9m_c > 10^6) = NaN;
    file.T_10m_c(file.T_10m_c > 10^6) = NaN;
    file.T_11m_c(file.T_11m_c > 10^6) = NaN;
    file.T_12m_c(file.T_12m_c > 10^6) = NaN;
    file.T_13m_c(file.T_13m_c > 10^6) = NaN;
    file.T_14m_c(file.T_14m_c > 10^6) = NaN;
    file.T_15m_c(file.T_15m_c > 10^6) = NaN;
    file.T_16m_c(file.T_16m_c > 10^6) = NaN;
    file.T_17m_c(file.T_17m_c > 10^6) = NaN;
    file.T_18m_c(file.T_18m_c > 10^6) = NaN;
    file.T_19m_c(file.T_19m_c > 10^6) = NaN;
    file.T_20m_c(file.T_20m_c > 10^6) = NaN;
    file.u_10m_c(file.u_10m_c > 10^6) = NaN;
    file.v_10m_c(file.v_10m_c > 10^6) = NaN;
    file.w_10m_c(file.w_10m_c > 10^6) = NaN;
    file.u_1m_c(file.u_1m_c > 10^6) = NaN;
    file.v_1m_c(file.v_1m_c > 10^6) = NaN;
    file.w_1m_c(file.w_1m_c > 10^6) = NaN;
    file.u_3m_c(file.u_3m_c > 10^6) = NaN;
    file.v_3m_c(file.v_3m_c > 10^6) = NaN;
    file.w_3m_c(file.w_3m_c > 10^6) = NaN;
    file.Tirga_1m_c(file.Tirga_1m_c > 10^6) = NaN;
    data.t1m = [data.t1m; file.T_1m_c];
    data.t2m = [data.t2m; file.T_2m_c];
    data.t3m = [data.t3m; file.T_3m_c];
    data.t4m = [data.t4m; file.T_4m_c];
    data.t5m = [data.t5m; file.T_5m_c];
    data.t6m = [data.t6m; file.T_6m_c];
    data.t7m = [data.t7m; file.T_7m_c];
    data.t8m = [data.t8m; file.T_8m_c];
    data.t9m = [data.t9m; file.T_9m_c];
    data.t10m = [data.t10m; file.T_10m_c];
    data.t11m = [data.t11m; file.T_11m_c];
    data.t12m = [data.t12m; file.T_12m_c];
    data.t13m = [data.t13m; file.T_13m_c];
    data.t14m = [data.t14m; file.T_14m_c];
    data.t15m = [data.t15m; file.T_15m_c];
    data.t16m = [data.t16m; file.T_16m_c];
    data.t17m = [data.t17m; file.T_17m_c];
    data.t18m = [data.t18m; file.T_18m_c];
    data.t19m = [data.t19m; file.T_19m_c];
    data.t20m = [data.t20m; file.T_20m_c];
    data.rh10m = [data.rh10m; file.RH_10m_c];
    data.rh4m = [data.rh4m; file.RH_4m_c];

    tmp = sqrt(file.u_10m_c.^2+file.v_10m_c.^2+file.w_10m_c.^2);
    tmp = reshape(tmp,20*3600,1);
    data.u = [data.u; tmp];
    tmp = reshape(file.Tirga_1m_c,20*3600,1);
    data.tac1m = [data.tac1m; tmp];
    data.w10m = [data.w10m; sqrt(nanmean(file.u_10m_c).^2+nanmean(file.v_10m_c).^2+nanmean(file.w_10m_c).^2)'];
    data.w1m = [data.w1m; sqrt(nanmean(file.u_1m_c).^2+nanmean(file.v_1m_c).^2+nanmean(file.w_1m_c).^2)'];
    data.w3m = [data.w3m; sqrt(nanmean(file.u_3m_c).^2+nanmean(file.v_3m_c).^2+nanmean(file.w_3m_c).^2)'];
    data.pr = [data.pr; file.P_c_c(1,:)'];
end

data.dn20 = floor(data.dn(1)):1/1728000:ceil(data.dn(end)); data.dn20 = data.dn20(1:end-1)';
ind_twr20 = find(data.dn20 >= time1 & data.dn20 < time2);

cd /Users/ccox/Documents/Projects/2021/splash/science/sublimation/sos_case_study/

ind_twr = find(data.dn >= time1 & data.dn <= time2);
t1m = data.t1m(ind_twr);
t2m = data.t2m(ind_twr);
t3m = data.t3m(ind_twr);
t4m = data.t4m(ind_twr);
t5m = data.t5m(ind_twr);
t6m = data.t6m(ind_twr);
t7m = data.t7m(ind_twr);
t8m = data.t8m(ind_twr);
t9m = data.t9m(ind_twr);
t10m = data.t10m(ind_twr);
t11m = data.t11m(ind_twr);
t12m = data.t12m(ind_twr);
t13m = data.t13m(ind_twr);
t14m = data.t14m(ind_twr);
t15m = data.t15m(ind_twr);
t16m = data.t16m(ind_twr);
t17m = data.t17m(ind_twr);
t18m = data.t18m(ind_twr);
t19m = data.t19m(ind_twr);
t20m = data.t20m(ind_twr);
u = data.w3m(ind_twr);
u10 = data.w10m(ind_twr);
pr = data.pr(ind_twr);
rh4m = data.rh4m(ind_twr);
rh10m = data.rh10m(ind_twr);

% A couple things from the slow file
%ustar = mean([nanmean(sqrt(sosf.u_w__10m_c(indsos).^2)) nanmean(sqrt(sosf.u_w__10m_ue(indsos).^2)) nanmean(sqrt(sosf.u_w__10m_d(indsos).^2)) nanmean(sqrt(sosf.u_w__10m_uw(indsos).^2))]);
F = mean([sosf.w_tc__3m_c(indsos);sosf.w_tc__3m_ue(indsos);sosf.w_tc__3m_d(indsos);sosf.w_tc__3m_uw(indsos)]).*nanmean(rho).*cp; %nanmean(a50.Hs(ind50));
ustar = mean([nanmean(sqrt(sosf.u_w__3m_c(indsos).^2)) nanmean(sqrt(sosf.u_w__3m_ue(indsos).^2)) nanmean(sqrt(sosf.u_w__3m_d(indsos).^2)) nanmean(sqrt(sosf.u_w__3m_uw(indsos).^2))]);

% Calculate some things
%
% Integral length scale. I don't actually know optimal integration time. If
% the value is too low then the autocorrelation function doesn't converge
% but if it is too high then advective processes are treated as turbulence.
% Same issue as the integtration window for the turbulence. This can be
% seen by lookng at the full two data sos data set as a f(integ_time), like
% so:
% integ_time = 20:30;
% L_sos = NaN(length(integ_time),length(data.dn));
% for k = 1:length(integ_time)
%     k
%     [L_tmp,L_dn] = integral_length_scale(naninterp(data.u),data.dn20,20,integ_time(k)); 
%     L_tmp = interp1(L_dn,L_tmp,data.dn);
%     L_sos(k,:) = L_tmp;
% end
% plot(integ_time,nanmean(L_sos,2));
% % The function increase to around 15 min, then plateaus ~20-40 min as expected, then increases again.
% L_sosm = nanmean(L_sos(:,ind_twr),1)';% nanmean(L_sos(20:22,ind_twr))';
% nanmean(L_sosm)
[L_tmp,L_dn] = integral_length_scale(naninterp(data.u),data.dn20,20,10); 
L_sosm = nanmean(L_tmp);

Tref = nanmean([nanmean(t6m) nanmean(t7m) nanmean(t8m) nanmean(t9m) nanmean(t10m) ...
             nanmean(t11m) nanmean(t12m) nanmean(t13m) nanmean(t14m) nanmean(t15m) nanmean(t16m) nanmean(t17m) nanmean(t18m) nanmean(t19m) nanmean(t20m)]);
[rho_dry_sos, rho_sos] = calc_rho_air(Tref+K0,rh10m,pr);
kair_sos = calc_kair(Tref);
mu_sos = calc_dnyamic_viscocity(Tref+K0);
Pr_sos = calc_prandtl(cp,mu_sos,kair_sos);
Re_sos = calc_reynolds_number(rho_sos,u10,L_sosm,mu_sos);
Nu_sos = calc_nusselt(Pr_sos,Re_sos); 

% Calculate convective/conductive heat equation
[Q,h] = convective_heat(Tref+K0,Tref.*0+K0,L_sosm,Nu_sos,kair_sos,1);
%F = nanmean(Q);

% 1D transient heat conduction
xht = 6;
dx = 1e-2;
x = 0:dx:xht; % [m]

% Calculate tstar
alpha = calc_thermal_diffusivity(kair_sos,rho_sos,cp);
tstar = calc_tstar(h,alpha,rho_sos,cp,1); nanmean(tstar)
% Distribution of t* following Haussecker Eq. 13
t = 0:0.1:14;
P = P_t_tstar(nanmean(tstar),t);

% Now make a distribution of many random values that meets the criteria
% It doesn't work well to use the mean t*
n_sims =5000;
exp_values = exprnd(nanmean(tstar), n_sims, 1);

clf;
hb = histogram(exp_values,0:0.1:20); hb.Normalization = 'pdf';
hold on; grid on;
plot(t,P,'linewidth',2); 
xlabel('{\it t} [s]'); ylabel('P({\it t},{\it t*})');
set(gca,'xtick',0:2:14);
xlim([0 14]);
legend('N=1000 samples','Theoretical Distribution');
set(gca,'fontsize',16);
vline(nanmean(tstar));
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0 0 20 18]);
print('-dpng','-r300','probs.png');


% now I can calculate a stable profile using many t that are distributed in
% this theoretical way. I can then extract a temperature distribution at
% any height. I will do so and extract the temeprature at the height of the
% 1 m and 2 m, etc on the tower and compare.

% similarly, this should help me determine whether the observed fluctuations are tied
% to SRT as I should be able to find out how frequently the stable layer
% reaches the sensor height.

Uref = nanmean(u); %nanmean(a50.wspd_vec_mean(ind));
Tsfc = 0;
T = Tref;
A = 1;

Tcalcs_all = NaN(xht/dx,n_sims);
Tcalcs_turb = NaN(xht/dx,n_sims);
Tcalcs_visc = NaN(xht/dx,n_sims);

tic
for k = 1:n_sims

    disp([num2str((round(k/length(exp_values)*10000))/100),'% Complete...']);
        
    if k == 1 % calculate alpha_c
        [Tout,x,ustar_ret,alpha_c_ret] = diff1d(T,Tsfc,exp_values(k),xht,dx,A,nanmean(rho_sos),F,Uref,9,3,ustar,0.2,0);
    else % use the initial calc to save time
        [Tout,x,ustar_ret,alpha_c_ret] = diff1d(T,Tsfc,exp_values(k),xht,dx,A,nanmean(rho_sos),F,Uref,9,3,ustar,alpha_c_ret,0);
    end
    Tcalcs_all(:,k) = Tout;

    [Tout_visc,x,ustar_ret,alpha_c_ret] = diff1d(T,Tsfc,exp_values(k),xht,dx,A,nanmean(rho_sos),F,Uref,9,1,ustar,alpha_c_ret,0);
    Tcalcs_visc(:,k) = Tout_visc;

    [Tout_turb,x,ustar_ret,alpha_c_ret] = diff1d(T,Tsfc,exp_values(k),xht,dx,A,nanmean(rho_sos),F,Uref,9,2,ustar,alpha_c_ret,0);
    Tcalcs_turb(:,k) = Tout_turb;

end
toc
%save Tcalcs_Apr15.mat Tcalcs_all Tcalcs_visc Tcalcs_turb exp_values
load Tcalcs_Apr15.mat



% 1d semi-infinite solid
x = 0:dx:xht; % [m]
sis1d = (Tsfc + (T - Tsfc) .* erf(x./(2*sqrt(nanmean(alpha).*nanmean(tstar)))));
h3=plot(nanmean(Tcalcs_visc'),x(1:end-1),'color',cmap(2,:),'linewidth',2); hold on;
plot(sis1d,x,'color',cmap(1,:),'linewidth',2)




ii = find(Tcalcs_all(30,:) == -999 | Tcalcs_turb(30,:) == -999);
Tcalcs_turb(:,ii) = NaN;
Tcalcs_visc(:,ii) = NaN;
Tcalcs_all(:,ii) = NaN;

clf;
cmap=lines;
cmap2 = (colormap(sky(n_sims)));
[sortdata sortdata1] = sort(exp_values,'descend');
for k = 1:n_sims
    h0=plot(Tcalcs_all(:,sortdata1(k)),x,'color',cmap2(k,:),'linewidth',2); hold on; 
    hold on;
end
h1=plot(nanmean(Tcalcs_all'),x,'c','linewidth',2); hold on;
h2=plot(prctile(Tcalcs_all',25),x,'c--','linewidth',2); hold on;
h2=plot(prctile(Tcalcs_all',75),x,'c--','linewidth',2); hold on;

h3=plot(nanmean(Tcalcs_visc'),x,'color',cmap(2,:),'linewidth',2); hold on;
h4=plot(nanmean(Tcalcs_turb'),x,'color',cmap(3,:),'linewidth',2); hold on;

% 0.84 m of snow reported beneath the SOS Tower "c" (same as temps) at 4/18/23 23:30 UTC
h5=errorbar([nanmean(t1m) nanmean(t2m) nanmean(t3m) nanmean(t4m) nanmean(t5m) nanmean(t6m) nanmean(t7m) nanmean(t8m) nanmean(t9m) nanmean(t10m)], ...
         [1:10]-0.84, ...
         [0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05], ...
         [0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05], ...
         [nanstd(t1m) nanstd(t2m) nanstd(t3m) nanstd(t4m) nanstd(t5m) nanstd(t6m) nanstd(t7m) nanstd(t8m) nanstd(t9m) nanstd(t10m)], ...
         [nanstd(t1m) nanstd(t2m) nanstd(t3m) nanstd(t4m) nanstd(t5m) nanstd(t6m) nanstd(t7m) nanstd(t8m) nanstd(t9m) nanstd(t10m)], ...
         'k.','markersize',24);
xlabel('Temperature [C]'); ylabel('Height above snow [m]');
grid on;
legend([h0 h1 h2 h3 h4 h5],'Modeled Profiles','Model Mean','Model 25th/75th percentile','Viscous Only','No Form Drag','Measured','location','northwest');
colormap(cmap2);        % Use the custom colormap (cmap2)
c = colorbar;           % Add a colorbar
c.Ticks = linspace(0, 10, length(t)); % Set tick marks to match the number of lines
c.TickLabels = fliplr((0:10:length(t)).*0.1);     % Label ticks with corresponding line numbers
c.Label.String = '{\itt} [s]';
%set(gca,'yscale','log');
set(gca,'fontsize',16);
ylim([0 3.5]); xlim([0 8]);
set(gca,'ytick',0:.25:3.5,'xtick',0:1:8);
%set(gca,'ytick',[0 0.1 0.25 0.5 0.75 1:6],'xtick',0:1:8);
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0 0 20 18]);
print('-dpng','-r300','/Users/ccox/Documents/Projects/2021/splash/science/sublimation/pub/figs/Figure_11lin.png');






cd /Users/ccox/Documents/Projects/2021/splash/Simba_NOAA0101/data_processing/sledd/
ts=get_avery_raw; cd /Users/ccox/Documents/Projects/2021/splash/science/sublimation/sos_case_study/

ind = find(floor(ts.dn) == datenum([2023 4 18 0 0 0]));
rind = find(ts.range > 1.18);
tsimba = (ts.t(:,ind(end))+273.15);

Tsfc = 273.15;
T = nanmean(t4m)+273.15;
clf;
cmap=lines;
%Tcalc1 = Tsfc + (T - Tsfc) .* erf(x ./ (2*sqrt(alpha.*tstar)));
%plot(mean(Tcalc1-273.15),x); set(gca,'yscale','log'); hold on;
%Tcalc2 = T2 + (T1 - T2) .* erf(x ./ (2*sqrt(alpha*tstar)));
Tcalc2max = Tsfc + (T - Tsfc) .* exp( ((-(h.*A)./(rho.*cp.*x))).*prctile(tstar,90) );
Tcalc2min = Tsfc + (T - Tsfc) .* exp( ((-(h.*A)./(rho.*cp.*x))).*prctile(tstar,10) );
Tcalc2mean = Tsfc + (T - Tsfc) .* exp( ((-(h.*A)./(rho.*cp.*x))).*nanmean(tstar) );
%plot(mean(Tcalc2max)-273.15,x); set(gca,'yscale','lin'); hold on;
%plot(mean(Tcalc2min)-273.15,x); set(gca,'yscale','lin'); hold on;
h1=plot(nanmean(Tcalc2mean)-273.15,x,'k','color',cmap(1,:),'linewidth',2); hold on; 
h2=plot(nanmean(Tcalc2min)-273.15,x,'k--','color',cmap(1,:),'linewidth',1); 
h2=plot(nanmean(Tcalc2max)-273.15,x,'k--','color',cmap(1,:),'linewidth',1);
%plot(tsimba-273.15,ts.range-1.18,'color',cmap(2,:),'linewidth',1)
set(gca,'yscale','lin');



h3=errorbar([mean(t1m) mean(t2m) mean(t3m) mean(t4m)], ...
         [0.11 1.11 2.11 3.11], ...
         [std(t1m) std(t2m) std(t3m) std(t4m)], ...
         [0.04 0.04 0.04 0.04], ...
         'both','.','markersize',24,'color',cmap(2,:));
ylim([0 3.5]);
xlabel('Temperature [C]'); ylabel('Height above snow [m]');
grid on;
legend([h1 h2 h3],'Mean Modeled Profile','10th/90th percentile','Measured','location','northwest');

set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0 0 20 18]);
print -dpng -r300 /Users/ccox/Documents/Projects/2021/splash/science/sublimation/pub/figs/Figure_9.png



























Lall = interp1_noerror_incNaN(intLens.DN,intLens.L,a50.dn);
Qs = NaN.*a50.dn;
Prs = NaN.*a50.dn;
for k = 1:length(a50.dn)

Rsp = 287; % gas constant
T = a50.temp(k)+273.15; % air temp for the case  study
Tsfc = a50.skin_temp_surface(k)+273.15; 
p = a50.atmos_pressure(k) * 100; % air pressure for  the case study
rho = p / (Rsp*T); % air density for the case study

u = a50.wspd_vec_mean(k); % wind speed  mean(a50.wspd_w_std(ind)); 

% sutherland's formula for dynanmic viscocity of air
mu = 1.716e-5*(T/273.15)^(3/2)*((273.15+111)/(T+111));

% give me a range of characteristic lengths (m)
L = Lall(k);

% thermal conductivity of air; formula suggested by chatgpt...
kair = 0.0257*(1+0.0031*(T-300));

% Reynolds Number
Re = (rho.*u.*L) ./ mu;
% Prandtl Number
cp = 1005; % specific heat, J/(kg K)
Pr = (cp * mu) / kair;
Prs(k) = Pr;
% Nusselt Number for boundary layer turbulent flow; e.g., here:
% http://www.mhtl.uwaterloo.ca/courses/ece309_mechatronics/lectures/pdffiles/chapter6.pdf
% Nu = C * Re^m * Pr^n, but what are the appropriate coefficients?
% Lienhard (2020): https://doi.org/10.1115/1.4046795
Cf = 0.455 ./ (log(0.06.*Re)).^2;
Nu = (Re.*Pr.*(Cf./2)) ./ (1+12.7.*(Pr.^(2/3)-1).*sqrt(Cf./2));
% back out the % convective heat transfer coeficient, h, from Nu (https://en.wikipedia.org/wiki/Nusselt_number)
h = (Nu * kair) ./ L;
% convective heat transfer equation, Q in W/m2
Q = h .* 1 .* (T-Tsfc); % 1 is A is 1 m2 

Qs(k) = Q;

end





clf;
plot(L,Q,'linewidth',2); hold on;
ylabel('Heat Transfer [W m^-^2]'); xlabel('L [m]'); grid on; set(gca,'xscale','log');
plot(mean(Li(ind_Li)),abs(mean(data10.Hs(ind))),'.','markersize',16);
set(gca,'fontsize',16);
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0 0 20 18]);


exp_values = [];
for k = 1:length(tstar)

    t = 0:0.1:10;
    lambda = 1/tstar(k);
    % by this equation, tstar is understood as the e-folding time of a log-linear probability distribution of t
    P = lambda .* exp(-lambda.*t); % Eq. 13  % trapz(t,P) = 1
    
    % Now make a distribution of many random values that meets the criteria
    exp_values = [exp_values; exprnd(1/lambda, 20, 1)];

end
t = 0:0.1:10;
lambda = 1/nanmean(tstar);
% by this equation, tstar is understood as the e-folding time of a log-linear probability distribution of t
P = lambda .* exp(-lambda.*t); % Eq. 13  % trapz(t,P) = 1

% Now make a distribution of many random values that meets the criteria
exp_values = exprnd(1/lambda, 12000, 1);

