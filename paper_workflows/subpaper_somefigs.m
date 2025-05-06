
%% load data
load('/Users/ccox/Documents/Projects/2021/splash/_level2/a30dbl.mat')
a30 = data10; a30m = datam;
load('/Users/ccox/Documents/Projects/2021/splash/_level2/a50dbl.mat')
a50 = data10; a50m = datam;
clearallbut a30 a50 a30m a50m

a30m.snow_depth(a30m.snow_depth_qc ~= 0) = NaN;
a30.azimuth(a30.azimuth_qc ~= 0) = NaN;
a30.sr50_dist(a30.sr50_dist_qc ~= 0) = NaN;
a30.atmos_pressure(a30.atmos_pressure_qc ~= 0) = NaN;
a30.temp(a30.temp_qc ~= 0) = NaN;
a30.rh(a30.rh_qc ~= 0) = NaN;
a30.brightness_temp_surface(a30.brightness_temp_surface_qc ~= 0) = NaN;
a30.skin_temp_surface(a30.skin_temp_surface_qc ~= 0) = NaN;
a30.subsurface_heat_flux_A(a30.subsurface_heat_flux_A_qc ~= 0) = NaN;
a30.subsurface_heat_flux_B(a30.subsurface_heat_flux_B_qc ~= 0) = NaN;
a30.wspd_u_mean(a30.wspd_u_mean_qc ~= 0) = NaN;
a30.wspd_v_mean(a30.wspd_v_mean_qc ~= 0) = NaN;
a30.wspd_w_mean(a30.wspd_w_mean_qc ~= 0) = NaN;
a30.wspd_vec_mean(a30.wspd_vec_mean_qc ~= 0) = NaN;
a30.wdir_vec_mean(a30.wdir_vec_mean_qc ~= 0) = NaN;
a30.h2o_licor(a30.h2o_licor_qc ~= 0) = NaN;
a30.co2_licor(a30.co2_licor_qc ~= 0) = NaN;
a30.down_long_hemisp(a30.down_long_hemisp_qc ~= 0) = NaN;
a30.down_short_hemisp(a30.down_short_hemisp_qc ~= 0) = NaN;
a30.up_long_hemisp(a30.up_long_hemisp_qc ~= 0) = NaN;
a30.up_short_hemisp(a30.up_short_hemisp_qc ~= 0) = NaN;
a30.s_soil(a30.s_soil_qc ~= 0) = NaN;
a30.c_soil(a30.c_soil_qc ~= 0) = NaN;
a30.k_eff_soil(a30.k_eff_soil_qc ~= 0) = NaN;
a30.soil_vwc_5cm(a30.soil_vwc_5cm_qc ~= 0) = NaN;
a30.soil_vwc_10cm(a30.soil_vwc_10cm_qc ~= 0) = NaN;
a30.soil_vwc_20cm(a30.soil_vwc_20cm_qc ~= 0) = NaN;
a30.soil_vwc_30cm(a30.soil_vwc_30cm_qc ~= 0) = NaN;
a30.soil_vwc_40cm(a30.soil_vwc_40cm_qc ~= 0) = NaN;
a30.soil_vwc_50cm(a30.soil_vwc_50cm_qc ~= 0) = NaN;
a30.soil_t_5cm(a30.soil_t_5cm_qc ~= 0) = NaN;
a30.soil_t_10cm(a30.soil_t_10cm_qc ~= 0) = NaN;
a30.soil_t_20cm(a30.soil_t_20cm_qc ~= 0) = NaN;
a30.soil_t_30cm(a30.soil_t_30cm_qc ~= 0) = NaN;
a30.soil_t_40cm(a30.soil_t_40cm_qc ~= 0) = NaN;
a30.soil_t_50cm(a30.soil_t_50cm_qc ~= 0) = NaN;
a30.Hl(a30.turbulence_qc ~= 0) = NaN;
a30.Hs(a30.turbulence_qc ~= 0) = NaN;
a30.bulk_Hl(a30.bulk_qc ~= 0) = NaN;
a30.bulk_Hs(a30.bulk_qc ~= 0) = NaN;
a30.Cd(a30.turbulence_qc ~= 0) = NaN;
a30.ustar(a30.turbulence_qc ~= 0) = NaN;
a30.zeta_level_n(a30.turbulence_qc ~= 0) = NaN;
a30.mixing_ratio(a30.mixing_ratio_qc ~= 0) = NaN;
a30.atmos_pressure(a30.atmos_pressure_qc ~= 0) = NaN;
a30.bulk_ustar(a30.bulk_qc ~= 0) = NaN;
a30.bulk_tstar(a30.bulk_qc ~= 0) = NaN;

a50m.snow_depth(a50m.snow_depth_qc ~= 0) = NaN;
a50.azimuth(a50.azimuth_qc ~= 0) = NaN;
a50.sr50_dist(a50.sr50_dist_qc ~= 0) = NaN;
a50.atmos_pressure(a50.atmos_pressure_qc ~= 0) = NaN;
a50.temp(a50.temp_qc ~= 0) = NaN;
a50.rh(a50.rh_qc ~= 0) = NaN;
a50.brightness_temp_surface(a50.brightness_temp_surface_qc ~= 0) = NaN;
a50.skin_temp_surface(a50.skin_temp_surface_qc ~= 0) = NaN;
a50.subsurface_heat_flux_A(a50.subsurface_heat_flux_A_qc ~= 0) = NaN;
a50.subsurface_heat_flux_B(a50.subsurface_heat_flux_B_qc ~= 0) = NaN;
a50.wspd_u_mean(a50.wspd_u_mean_qc ~= 0) = NaN;
a50.wspd_v_mean(a50.wspd_v_mean_qc ~= 0) = NaN;
a50.wspd_w_mean(a50.wspd_w_mean_qc ~= 0) = NaN;
a50.wspd_vec_mean(a50.wspd_vec_mean_qc ~= 0) = NaN;
a50.wdir_vec_mean(a50.wdir_vec_mean_qc ~= 0) = NaN;
a50.h2o_licor(a50.h2o_licor_qc ~= 0) = NaN;
a50.co2_licor(a50.co2_licor_qc ~= 0) = NaN;
a50.down_long_hemisp(a50.down_long_hemisp_qc ~= 0) = NaN;
a50.down_short_hemisp(a50.down_short_hemisp_qc ~= 0) = NaN;
a50.up_long_hemisp(a50.up_long_hemisp_qc ~= 0) = NaN;
a50.up_short_hemisp(a50.up_short_hemisp_qc ~= 0) = NaN;
a50.s_soil(a50.s_soil_qc ~= 0) = NaN;
a50.c_soil(a50.c_soil_qc ~= 0) = NaN;
a50.k_eff_soil(a50.k_eff_soil_qc ~= 0) = NaN;
a50.soil_vwc_5cm(a50.soil_vwc_5cm_qc ~= 0) = NaN;
a50.soil_vwc_10cm(a50.soil_vwc_10cm_qc ~= 0) = NaN;
a50.soil_vwc_20cm(a50.soil_vwc_20cm_qc ~= 0) = NaN;
a50.soil_vwc_30cm(a50.soil_vwc_30cm_qc ~= 0) = NaN;
a50.soil_vwc_40cm(a50.soil_vwc_40cm_qc ~= 0) = NaN;
a50.soil_vwc_50cm(a50.soil_vwc_50cm_qc ~= 0) = NaN;
a50.soil_t_5cm(a50.soil_t_5cm_qc ~= 0) = NaN;
a50.soil_t_10cm(a50.soil_t_10cm_qc ~= 0) = NaN;
a50.soil_t_20cm(a50.soil_t_20cm_qc ~= 0) = NaN;
a50.soil_t_30cm(a50.soil_t_30cm_qc ~= 0) = NaN;
a50.soil_t_40cm(a50.soil_t_40cm_qc ~= 0) = NaN;
a50.soil_t_50cm(a50.soil_t_50cm_qc ~= 0) = NaN;
a50.Hl(a50.turbulence_qc ~= 0) = NaN;
a50.Hs(a50.turbulence_qc ~= 0) = NaN;
a50.bulk_Hl(a50.bulk_qc ~= 0) = NaN;
a50.bulk_Hs(a50.bulk_qc ~= 0) = NaN;
a50.Cd(a50.turbulence_qc ~= 0) = NaN;
a50.ustar(a50.turbulence_qc ~= 0) = NaN;
a50.zeta_level_n(a50.turbulence_qc ~= 0) = NaN;
a50.mixing_ratio(a50.mixing_ratio_qc ~= 0) = NaN;
a50.atmos_pressure(a50.atmos_pressure_qc ~= 0) = NaN;
a50.bulk_ustar(a50.bulk_qc ~= 0) = NaN;
a50.bulk_tstar(a50.bulk_qc ~= 0) = NaN;

load /Users/ccox/Documents/Projects/2021/splash/_barr_data/barr2024.mat

addpath(genpath('/Users/ccox/Documents/Projects/2021/splash/science/sublimation/mfiles/'));
addpath(genpath('/Users/ccox/Documents/Projects/2021/mosaic/science/Cox_etal_bulk_sensitivity/'));

cd /Users/ccox/Documents/Projects/2021/splash/science/sublimation/figures/

% if you want...
%a50.dn = a50.dn - 7/24; % to mst
%a30.dn = a30.dn - 7/24; % to mst


%% Z0


a30snow_depth = a30.snow_depth;
a30snow_depth(a30.snow_depth_qc ~= 0) = NaN;
a30snow_depth = a30snow_depth.*0.01;                               % convert to meters
a30snow_depth = naninterp(a30snow_depth); %snow_depth.rolling(600, min_periods=10).mean();  % fill nans for bulk calc only

a50snow_depth = a50.snow_depth;
a50snow_depth(a50.snow_depth_qc ~= 0) = NaN;
a50snow_depth = a50snow_depth.*0.01;                               % convert to meters
a50snow_depth = naninterp(a50snow_depth); %snow_depth.rolling(600, min_periods=10).mean();  % fill nans for bulk calc only

tmp = interp1(a30.dn,a30snow_depth,a50.dn);
a50snow_depth(1:1.6244e4) = tmp(1:1.6244e4).*0.88;


% limits to screening z0
limMax = 10^-1;
limMin = 10^-10;
wMin   = 0;


[z030,zT30] = calc_z0(naninterp(4.62 - a30snow_depth./100),a30.Cd,a30.zeta_level_n,3);
z030(z030 > limMax) = NaN;
z030(z030 < limMin) = NaN;
z030(a30.wspd_vec_mean < wMin) = NaN;

[z050,zT50] = calc_z0(naninterp(4.62 - a50snow_depth./100),a50.Cd,a50.zeta_level_n,3);
z050(z050 > limMax) = NaN;
z050(z050 < limMin) = NaN;
z050(a50.wspd_vec_mean < wMin) = NaN;


% Q=a30.mixing_ratio./1000; 
% Q = Q./(Q+1);
% es=(1.0003+4.18e-6.*a30.atmos_pressure).*6.1115.*exp(22.452.*a30.skin_temp_surface./(a30.skin_temp_surface+272.55)); % saturation vapor pressure
% Qs=es.*622./(a30.atmos_pressure-.378.*es)./1000;  
% zT30 = calc_zT(naninterp(4.62 - a30snow_depth./100),naninterp(2.6 - a30snow_depth./100),a30.Cd,a30.ustar,a30.Wq_csp,(Qs-Q),a30.zeta_level_n,3);
% zT30(zT30 > 1e1) = NaN;
% zT30(zT30 < 1e-8) = NaN;
% zT30(a30.wspd_vec_mean < wMin) = NaN;
% 
% Q=a50.mixing_ratio./1000; 
% Q = Q./(Q+1);
% es=(1.0003+4.18e-6.*a50.atmos_pressure).*6.1115.*exp(22.452.*a50.skin_temp_surface./(a50.skin_temp_surface+272.55)); % saturation vapor pressure
% Qs=es.*622./(a50.atmos_pressure-.378.*es)./1000;  
% zT50 = calc_zT(naninterp(4.62 - a50snow_depth./100),naninterp(2.6 - a50snow_depth./100),a50.Cd,a50.ustar,a50.Wq_csp,(Qs-Q),a50.zeta_level_n,3);
% zT50(zT50 > 1e1) = NaN;
% zT50(zT50 < 1e-8) = NaN;
% zT50(a50.wspd_vec_mean < wMin) = NaN;

mu = 1.716e-5 .* ((a30.temp+273.15)./273.15).^3/2 .* (273.15+110.4)./((a30.temp+273.15)+110.4);
nu = mu ./ (a30.atmos_pressure.*100./(287.1.*(a30.temp+273.15).*(1+0.61*a30.mixing_ratio)));
[zT30, zQ30] = andreas_zt(a30.ustar, z030, a30.temp+273.15, nu);
mu = 1.716e-5 .* ((a50.temp+273.15)./273.15).^3/2 .* (273.15+110.4)./((a50.temp+273.15)+110.4);
nu = 1.81e-5 ./ (a50.atmos_pressure.*100./(287.1.*(a50.temp+273.15).*(1+0.61*a50.mixing_ratio)));
[zT50, zQ50] = andreas_zt(a50.ustar, z050, a50.temp+273.15, nu);

[yy mm30 u u u u] = datevec(a30.dn);
[yy mm50 u u u u] = datevec(a50.dn);
z030m = [];
z050m = [];
zT30m = [];
zT50m = [];
zQ30m = [];
zQ50m = [];
for h = 1:12
    a = z030(mm30==h); a(isnan(a))=[];
    z030m = [z030m median(a)];
    a = z050(mm50==h); a(isnan(a))=[];
    z050m = [z050m median(a)];
    a = z030(mm30==h); a(isnan(a))=[];
    zT30m = [zT30m median(a)];
    a = zT50(mm50==h); a(isnan(a))=[];
    zT50m = [zT50m median(a)];
    a = zQ30(mm30==h); a(isnan(a))=[];
    zQ30m = [zQ30m median(a)];
    a = zT50(mm50==h); a(isnan(a))=[];
    zQ50m = [zQ50m median(a)];
end

vals = NaN(11,length(a30.dn));
for k = 1:length(a30.dn)    
    vals(1,k) = a30.wspd_vec_mean(k);
    vals(2,k) = a30.skin_temp_surface(k);
    vals(3,k) = a30.temp(k);
    vals(4,k) = a30.mixing_ratio(k)/1000;
    vals(5,k) = 600;
    vals(6,k) = a30.atmos_pressure(k);
    vals(7,k) = 4.62-a30snow_depth(k);
    vals(8,k) = 2.89-a30snow_depth(k);
    vals(9,k) = 2.60-a30snow_depth(k);
    vals(10,k) = a30.rh(k);
    vals(11,k) = a30.soil_vwc_10cm(k);
end

zT302 = movmedian(naninterp(zT30),144*3); zT302(isnan(zT302)) = nanmedian(zT302);
zQ302 = movmedian(naninterp(zQ30),144*3); zQ302(isnan(zQ302)) = nanmedian(zQ302);
z0302 = movmedian(naninterp(z030),144*3); zT302(isnan(zT302)) = nanmedian(zT302);
hs_30 = NaN(length(a30.dn),1);
hl_30 = NaN(length(a30.dn),1);
for k = 1:length(a30.dn)
    tmp = vals(:,k);
    %y=cor_ice_A10_specifyZogs_splash(tmp,z0302(k),ztzq,ztzq,1,a30.snow_flag(k),'asfs30');
    y=cor_ice_A10_specifyZogs_splash(tmp,z0302(k),zT302(k),zQ302(k),1,a30.snow_flag(k),'asfs30');
    hs_30(k) = y(1);
    hl_30(k) = y(2);
end

hs_30z0 = NaN(length(a30.dn),1);
hl_30z0 = NaN(length(a30.dn),1);
for k = 1:length(a30.dn)
    tmp = vals(:,k);
    %y=cor_ice_A10_specifyZogs_splash(tmp,z0302(k),ztzq,ztzq,1,a30.snow_flag(k),'asfs30');
    y=cor_ice_A10_specifyZogs_splash(tmp,z0302(k),1e-4,1e-4,1,a30.snow_flag(k),'asfs30');
    hs_30z0(k) = y(1);
    hl_30z0(k) = y(2);
end



vals = NaN(11,length(a50.dn));
for k = 1:length(a50.dn)    
    vals(1,k) = a50.wspd_vec_mean(k);
    vals(2,k) = a50.skin_temp_surface(k);
    vals(3,k) = a50.temp(k);
    vals(4,k) = a50.mixing_ratio(k)/1000;
    vals(5,k) = 600;
    vals(6,k) = a50.atmos_pressure(k);
    vals(7,k) = 4.62-a50snow_depth(k);
    vals(8,k) = 2.89-a50snow_depth(k);
    vals(9,k) = 2.60-a50snow_depth(k);
    vals(10,k) = a50.rh(k);
    vals(11,k) = a50.soil_vwc_10cm(k);
end

zT502 = movmedian(naninterp(zT50),144*3); zT502(isnan(zT502)) = nanmedian(zT502);
zQ502 = movmedian(naninterp(zQ50),144*3); zQ502(isnan(zQ502)) = nanmedian(zQ502);
z0502 = movmedian(naninterp(z050),144*3); z0502(isnan(z0502)) = nanmedian(z0502);
hs_50 = NaN(length(a50.dn),1);
hl_50 = NaN(length(a50.dn),1);
for k = 1:length(a50.dn)
    tmp = vals(:,k);
    %y=cor_ice_A10_specifyZogs_splash(tmp,z0502(k),ztzq,ztzq,1,a50.snow_flag(k),'asfs50');
    y=cor_ice_A10_specifyZogs_splash(tmp,z0502(k),zT502(k),zQ502(k),1,a50.snow_flag(k),'asfs50');
    hs_50(k) = y(1);
    hl_50(k) = y(2);
end

hs_50z0 = NaN(length(a50.dn),1);
hl_50z0 = NaN(length(a50.dn),1);
for k = 1:length(a50.dn)
    tmp = vals(:,k);
    %y=cor_ice_A10_specifyZogs_splash(tmp,z0502(k),ztzq,ztzq,1,a50.snow_flag(k),'asfs50');
    y=cor_ice_A10_specifyZogs_splash(tmp,z0502(k),1e-4,1e-4,1,a50.snow_flag(k),'asfs50');
    hs_50z0(k) = y(1);
    hl_50z0(k) = y(2);
end



% do a test where the heights are offset to increase 2, 4, 10 m



vals = NaN(11,length(a50.dn));
for k = 1:length(a50.dn)    
    vals(1,k) = a50.wspd_vec_mean(k);
    vals(2,k) = a50.skin_temp_surface(k);
    vals(3,k) = a50.temp(k);
    vals(4,k) = a50.mixing_ratio(k)/1000;
    vals(5,k) = 600;
    vals(6,k) = a50.atmos_pressure(k);
    vals(7,k) = 4.62-a50snow_depth(k);
    vals(8,k) = 2.89-a50snow_depth(k) + 2;
    vals(9,k) = 2.60-a50snow_depth(k) + 2;
    vals(10,k) = a50.rh(k);
    vals(11,k) = a50.soil_vwc_10cm(k);
end

hs_50p2 = NaN(length(a50.dn),1);
hl_50p2 = NaN(length(a50.dn),1);
for k = 1:length(a50.dn)
    tmp = vals(:,k);
    %y=cor_ice_A10_specifyZogs_splash(tmp,z0502(k),ztzq,ztzq,1,a50.snow_flag(k),'asfs50');
    y=cor_ice_A10_specifyZogs_splash(tmp,z0502(k),zT502(k),zQ502(k),1,a50.snow_flag(k),'asfs50');
    hs_50p2(k) = y(1);
    hl_50p2(k) = y(2);
end

vals = NaN(11,length(a50.dn));
for k = 1:length(a50.dn)    
    vals(1,k) = a50.wspd_vec_mean(k);
    vals(2,k) = a50.skin_temp_surface(k);
    vals(3,k) = a50.temp(k);
    vals(4,k) = a50.mixing_ratio(k)/1000;
    vals(5,k) = 600;
    vals(6,k) = a50.atmos_pressure(k);
    vals(7,k) = 4.62-a50snow_depth(k);
    vals(8,k) = 2.89-a50snow_depth(k) + 4;
    vals(9,k) = 2.60-a50snow_depth(k) + 4;
    vals(10,k) = a50.rh(k);
    vals(11,k) = a50.soil_vwc_10cm(k);
end

hs_50p4 = NaN(length(a50.dn),1);
hl_50p4 = NaN(length(a50.dn),1);
for k = 1:length(a50.dn)
    tmp = vals(:,k);
    %y=cor_ice_A10_specifyZogs_splash(tmp,z0502(k),ztzq,ztzq,1,a50.snow_flag(k),'asfs50');
    y=cor_ice_A10_specifyZogs_splash(tmp,z0502(k),zT502(k),zQ502(k),1,a50.snow_flag(k),'asfs50');
    hs_50p4(k) = y(1);
    hl_50p4(k) = y(2);
end

vals = NaN(11,length(a50.dn));
for k = 1:length(a50.dn)    
    vals(1,k) = a50.wspd_vec_mean(k);
    vals(2,k) = a50.skin_temp_surface(k);
    vals(3,k) = a50.temp(k);
    vals(4,k) = a50.mixing_ratio(k)/1000;
    vals(5,k) = 600;
    vals(6,k) = a50.atmos_pressure(k);
    vals(7,k) = 4.62-a50snow_depth(k);
    vals(8,k) = 2.89-a50snow_depth(k) + 10;
    vals(9,k) = 2.60-a50snow_depth(k) + 10;
    vals(10,k) = a50.rh(k);
    vals(11,k) = a50.soil_vwc_10cm(k);
end

hs_50p10 = NaN(length(a50.dn),1);
hl_50p10 = NaN(length(a50.dn),1);
for k = 1:length(a50.dn)
    tmp = vals(:,k);
    %y=cor_ice_A10_specifyZogs_splash(tmp,z0502(k),ztzq,ztzq,1,a50.snow_flag(k),'asfs50');
    y=cor_ice_A10_specifyZogs_splash(tmp,z0502(k),zT502(k),zQ502(k),1,a50.snow_flag(k),'asfs50');
    hs_50p10(k) = y(1);
    hl_50p10(k) = y(2);
end

vals = NaN(11,length(a50.dn));
for k = 1:length(a50.dn)    
    vals(1,k) = a50.wspd_vec_mean(k);
    vals(2,k) = a50.skin_temp_surface(k);
    vals(3,k) = a50.temp(k);
    vals(4,k) = a50.mixing_ratio(k)/1000;
    vals(5,k) = 600;
    vals(6,k) = a50.atmos_pressure(k);
    vals(7,k) = 4.62-a50snow_depth(k);
    vals(8,k) = 2.89-a50snow_depth(k) + 20;
    vals(9,k) = 2.60-a50snow_depth(k) + 20;
    vals(10,k) = a50.rh(k);
    vals(11,k) = a50.soil_vwc_10cm(k);
end

hs_50p20 = NaN(length(a50.dn),1);
hl_50p20 = NaN(length(a50.dn),1);
for k = 1:length(a50.dn)
    tmp = vals(:,k);
    %y=cor_ice_A10_specifyZogs_splash(tmp,z0502(k),ztzq,ztzq,1,a50.snow_flag(k),'asfs50');
    y=cor_ice_A10_specifyZogs_splash(tmp,z0502(k),zT502(k),zQ502(k),1,a50.snow_flag(k),'asfs50');
    hs_50p20(k) = y(1);
    hl_50p20(k) = y(2);
end

save /Users/ccox/Documents/Projects/2021/splash/science/sublimation/sos_case_study/offsettest.mat hs_50p20 hl_50p20 hs_50p10 hl_50p10 hs_50p4 hl_50p4 hs_50p2 hl_50p2 hs_50 hl_50


%% Figure 1

clf;
m_proj('lambert','long',[-130 -60],'lat',[25 55]);
m_coast('patch',[0.75 0.75 0.75]);
m_elev('contourf',[00:500:6000]);
m_grid('box','fancy','tickdir','in');
colormap(flipud(gray));
m_line(-106.9879,38.9562,'marker','.','markersize',40,'color','y','markerfacecolor','k','linewidth',2);
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0 0 8 6]);
print -dpng -r300 Figure_1inset.png

%% Figure 2


% Figure 1a
% albedo plot
a = csvread('/Users/ccox/Documents/Projects/2021/splash/outreach/AYP.csv');
k = csvread('/Users/ccox/Documents/Projects/2021/splash/outreach/KPA.csv');

adn = datenum([a(:,1) a(:,2) a(:,3) a(:,3).*0 a(:,3).*0 a(:,3).*0]);
kdn = datenum([k(:,1) k(:,2) k(:,3) k(:,3).*0 k(:,3).*0 k(:,3).*0]);
aalb = a(:,5); aalb(aalb < -1) = NaN;
kalb = k(:,5); kalb(kalb < -1) = NaN;

clf;
plot(kdn,kalb,'.-','linewidth',1,'markersize',10); hold on; 
plot(adn,aalb,'.-','linewidth',1,'markersize',10); 
set(gca,'xtick',datenum([2021 9 20 0 0 0]):20:datenum([2023 7 1 0 0 0]));
ylim([0 1]); xlim([datenum([2021 9 20 0 0 0]) datenum([2023 7 1 0 0 0])]); grid on;
ylabel('\alpha');
datetick('x','mmmdd yyyy','keepticks','keeplimits');
legend('KPA','AYP','location','north');
set(gca,'fontsize',16);
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0 0 50 14]);
print -dpng -r300 /Users/ccox/Documents/Projects/2021/splash/science/sublimation/pub/figs/Figure_2/Figure_2a.png


%% Figures 3, 4: Calculate Sublimation

[sub_mass_ec_ayp22,sub_mass_bulk_ayp22,sub_mass_be_ayp22] = calc_sub_precip(a50.temp,a50.skin_temp_surface,a50.bulk_Hl,a50.Hl,a50.snow_flag,barr.dn_wy22,a50.dn,'');
[sub_mass_ec_kpa22,sub_mass_bulk_kpa22,sub_mass_be_kpa22] = calc_sub_precip(a30.temp,a30.skin_temp_surface,a30.bulk_Hl,a30.Hl,a30.snow_flag,barr.dn_wy22,a30.dn,'');

sub_mass_ayp22_fill_ec = sub_mass_ec_ayp22;
sub_mass_ayp22_fill_ec(isnan(sub_mass_ayp22_fill_ec)) = sub_mass_ec_kpa22(isnan(sub_mass_ayp22_fill_ec));

sub_mass_kpa22_fill_ec = sub_mass_ec_kpa22;
sub_mass_kpa22_fill_ec(isnan(sub_mass_kpa22_fill_ec)) = sub_mass_ec_ayp22(isnan(sub_mass_kpa22_fill_ec));

sub_mass_ayp22_fill_bulk = sub_mass_bulk_ayp22;
sub_mass_ayp22_fill_bulk(isnan(sub_mass_ayp22_fill_bulk)) = sub_mass_bulk_kpa22(isnan(sub_mass_ayp22_fill_bulk));

sub_mass_kpa22_fill_bulk = sub_mass_bulk_kpa22;
sub_mass_kpa22_fill_bulk(isnan(sub_mass_kpa22_fill_bulk)) = sub_mass_bulk_ayp22(isnan(sub_mass_kpa22_fill_bulk));

sub_mass_kpa22_fill_ec(barr.dn_wy22 > datenum([2022 5 4 0 0 0])) = NaN;    
sub_mass_kpa22_fill_bulk(barr.dn_wy22 > datenum([2022 5 4 0 0 0])) = NaN;    

clf;
subplot(3,3,1)
plot(barr.dn_wy22,cumsum(sub_mass_ayp22_fill_ec),'r','linewidth',2); % this is what was sublimated
hold on;
plot(barr.dn_wy22,cumsum(sub_mass_ayp22_fill_bulk),'r:','linewidth',2); % this is what was sublimated
plot(barr.dn_wy22,cumsum(sub_mass_kpa22_fill_ec),'b','linewidth',2); % this is what was sublimated
plot(barr.dn_wy22,cumsum(sub_mass_kpa22_fill_bulk),'b:','linewidth',2); % this is what was sublimated
ylabel('Cumulative Sublimation [kg/m^2]');
ttl=title('a)                     Constant {\it z_0}, {\it z_Q}');ttl.Units = 'Normalize';ttl.Position(1) = 0;ttl.HorizontalAlignment = 'left';
legend('AYP EC Sublimated Mass','AYP Bulk Sublimated Mass','KPA EC Sublimated Mass','KPA Bulk Sublimated Mass','location','northwest');
grid on;
set(gca,'xtick',datenum([2021 12 1 0 0 0]):10:datenum([2022 5 15 0 0 0]));
xlim([datenum([2021 12 1 0 0 0]) datenum([2022 5 15 0 0 0])]); ylim([0 50]);
datetick('x','mmmdd','keepticks','keeplimits');
set(gca,'fontsize',12);
subplot(3,3,4)
plot(barr.dn_wy22,sub_mass_ayp22_fill_ec,'r','linewidth',1); % this is what was sublimated
hold on;
plot(barr.dn_wy22,sub_mass_ayp22_fill_bulk,'r:','linewidth',1); % this is what was sublimated
ylabel('Daily Sublimation [kg/m^2]');
legend('AYP EC Sublimated Mass','AYP Bulk Sublimated Mass','location','northwest');
ttl=title('b)');ttl.Units = 'Normalize';ttl.Position(1) = 0;ttl.HorizontalAlignment = 'left';
grid on;
set(gca,'xtick',datenum([2021 12 1 0 0 0]):10:datenum([2022 5 15 0 0 0]));
xlim([datenum([2021 12 1 0 0 0]) datenum([2022 5 15 0 0 0])]); ylim([0 2.5]);
datetick('x','mmmdd','keepticks','keeplimits');
set(gca,'fontsize',12);
subplot(3,3,7)
plot(barr.dn_wy22,sub_mass_kpa22_fill_ec,'b','linewidth',1); % this is what was sublimated
hold on;
plot(barr.dn_wy22,sub_mass_kpa22_fill_bulk,'b:','linewidth',1); % this is what was sublimated
ylabel('Daily Sublimation [kg/m^2]');
legend('KPA EC Sublimated Mass','KPA Bulk Sublimated Mass','location','northwest');
ttl=title('c)');ttl.Units = 'Normalize';ttl.Position(1) = 0;ttl.HorizontalAlignment = 'left';
grid on;
set(gca,'xtick',datenum([2021 12 1 0 0 0]):10:datenum([2022 5 15 0 0 0]));
xlim([datenum([2021 12 1 0 0 0]) datenum([2022 5 15 0 0 0])]); ylim([0 2.5]);
datetick('x','mmmdd','keepticks','keeplimits');
set(gca,'fontsize',12);



[sub_mass_ec_ayp22,sub_mass_bulk_ayp22,sub_mass_be_ayp22] = calc_sub_precip(a50.temp,a50.skin_temp_surface,hl_50z0,a50.Hl,a50.snow_flag,barr.dn_wy22,a50.dn,'');
[sub_mass_ec_kpa22,sub_mass_bulk_kpa22,sub_mass_be_kpa22] = calc_sub_precip(a30.temp,a30.skin_temp_surface,hl_30z0,a30.Hl,a30.snow_flag,barr.dn_wy22,a30.dn,'');

sub_mass_ayp22_fill_ec = sub_mass_ec_ayp22;
sub_mass_ayp22_fill_ec(isnan(sub_mass_ayp22_fill_ec)) = sub_mass_ec_kpa22(isnan(sub_mass_ayp22_fill_ec));
sub_mass_kpa22_fill_ec = sub_mass_ec_kpa22;
sub_mass_kpa22_fill_ec(isnan(sub_mass_kpa22_fill_ec)) = sub_mass_ec_ayp22(isnan(sub_mass_kpa22_fill_ec));
sub_mass_ayp22_fill_bulk = sub_mass_bulk_ayp22;
sub_mass_ayp22_fill_bulk(isnan(sub_mass_ayp22_fill_bulk)) = sub_mass_bulk_kpa22(isnan(sub_mass_ayp22_fill_bulk));
sub_mass_kpa22_fill_bulk = sub_mass_bulk_kpa22;
sub_mass_kpa22_fill_bulk(isnan(sub_mass_kpa22_fill_bulk)) = sub_mass_bulk_ayp22(isnan(sub_mass_kpa22_fill_bulk));
sub_mass_kpa22_fill_ec(barr.dn_wy22 > datenum([2022 5 4 0 0 0])) = NaN;    
sub_mass_kpa22_fill_bulk(barr.dn_wy22 > datenum([2022 5 4 0 0 0])) = NaN;    

subplot(3,3,2)
plot(barr.dn_wy22,cumsum(sub_mass_ayp22_fill_ec),'r','linewidth',2); % this is what was sublimated
hold on;
plot(barr.dn_wy22,cumsum(sub_mass_ayp22_fill_bulk),'r:','linewidth',2); % this is what was sublimated
plot(barr.dn_wy22,cumsum(sub_mass_kpa22_fill_ec),'b','linewidth',2); % this is what was sublimated
plot(barr.dn_wy22,cumsum(sub_mass_kpa22_fill_bulk),'b:','linewidth',2); % this is what was sublimated
ylabel('Cumulative Sublimation [kg/m^2]');
ttl=title('d)                Calculated {\it z_0}, Constant {\it z_Q}');ttl.Units = 'Normalize';ttl.Position(1) = 0;ttl.HorizontalAlignment = 'left';
legend('AYP EC Sublimated Mass','AYP Bulk Sublimated Mass','KPA EC Sublimated Mass','KPA Bulk Sublimated Mass','location','northwest');
grid on;
set(gca,'xtick',datenum([2021 12 1 0 0 0]):10:datenum([2022 5 15 0 0 0]));
xlim([datenum([2021 12 1 0 0 0]) datenum([2022 5 15 0 0 0])]); ylim([0 50]);
datetick('x','mmmdd','keepticks','keeplimits');
set(gca,'fontsize',12);
subplot(3,3,5)
plot(barr.dn_wy22,sub_mass_ayp22_fill_ec,'r','linewidth',1); % this is what was sublimated
hold on;
plot(barr.dn_wy22,sub_mass_ayp22_fill_bulk,'r:','linewidth',1); % this is what was sublimated
ylabel('Daily Sublimation [kg/m^2]');
legend('AYP EC Sublimated Mass','AYP Bulk Sublimated Mass','location','northwest');
ttl=title('e)');ttl.Units = 'Normalize';ttl.Position(1) = 0;ttl.HorizontalAlignment = 'left';
grid on;
set(gca,'xtick',datenum([2021 12 1 0 0 0]):10:datenum([2022 5 15 0 0 0]));
xlim([datenum([2021 12 1 0 0 0]) datenum([2022 5 15 0 0 0])]); ylim([0 2.5]);
datetick('x','mmmdd','keepticks','keeplimits');
set(gca,'fontsize',12);
subplot(3,3,8)
plot(barr.dn_wy22,sub_mass_kpa22_fill_ec,'b','linewidth',1); % this is what was sublimated
hold on;
plot(barr.dn_wy22,sub_mass_kpa22_fill_bulk,'b:','linewidth',1); % this is what was sublimated
ylabel('Daily Sublimation [kg/m^2]');
legend('KPA EC Sublimated Mass','KPA Bulk Sublimated Mass','location','northwest');
ttl=title('f)');ttl.Units = 'Normalize';ttl.Position(1) = 0;ttl.HorizontalAlignment = 'left';
grid on;
set(gca,'xtick',datenum([2021 12 1 0 0 0]):10:datenum([2022 5 15 0 0 0]));
xlim([datenum([2021 12 1 0 0 0]) datenum([2022 5 15 0 0 0])]); ylim([0 2.5]);
datetick('x','mmmdd','keepticks','keeplimits');
set(gca,'fontsize',12);


[sub_mass_ec_ayp22,sub_mass_bulk_ayp22,sub_mass_be_ayp22] = calc_sub_precip(a50.temp,a50.skin_temp_surface,hl_50,a50.Hl,a50.snow_flag,barr.dn_wy22,a50.dn,'');
[sub_mass_ec_kpa22,sub_mass_bulk_kpa22,sub_mass_be_kpa22] = calc_sub_precip(a30.temp,a30.skin_temp_surface,hl_30,a30.Hl,a30.snow_flag,barr.dn_wy22,a30.dn,'');

sub_mass_ayp22_fill_ec = sub_mass_ec_ayp22;
sub_mass_ayp22_fill_ec(isnan(sub_mass_ayp22_fill_ec)) = sub_mass_ec_kpa22(isnan(sub_mass_ayp22_fill_ec));
sub_mass_kpa22_fill_ec = sub_mass_ec_kpa22;
sub_mass_kpa22_fill_ec(isnan(sub_mass_kpa22_fill_ec)) = sub_mass_ec_ayp22(isnan(sub_mass_kpa22_fill_ec));
sub_mass_ayp22_fill_bulk = sub_mass_bulk_ayp22;
sub_mass_ayp22_fill_bulk(isnan(sub_mass_ayp22_fill_bulk)) = sub_mass_bulk_kpa22(isnan(sub_mass_ayp22_fill_bulk));
sub_mass_kpa22_fill_bulk = sub_mass_bulk_kpa22;
sub_mass_kpa22_fill_bulk(isnan(sub_mass_kpa22_fill_bulk)) = sub_mass_bulk_ayp22(isnan(sub_mass_kpa22_fill_bulk));
sub_mass_kpa22_fill_ec(barr.dn_wy22 > datenum([2022 5 4 0 0 0])) = NaN;    
sub_mass_kpa22_fill_bulk(barr.dn_wy22 > datenum([2022 5 4 0 0 0])) = NaN;    

subplot(3,3,3)
plot(barr.dn_wy22,cumsum(sub_mass_ayp22_fill_ec),'r','linewidth',2); % this is what was sublimated
hold on;
plot(barr.dn_wy22,cumsum(sub_mass_ayp22_fill_bulk),'r:','linewidth',2); % this is what was sublimated
plot(barr.dn_wy22,cumsum(sub_mass_kpa22_fill_ec),'b','linewidth',2); % this is what was sublimated
plot(barr.dn_wy22,cumsum(sub_mass_kpa22_fill_bulk),'b:','linewidth',2); % this is what was sublimated
ylabel('Cumulative Sublimation [kg/m^2]');
ttl=title('g)                     Calculated {\it z_0}, {\it z_Q}');ttl.Units = 'Normalize';ttl.Position(1) = 0;ttl.HorizontalAlignment = 'left';
legend('AYP EC Sublimated Mass','AYP Bulk Sublimated Mass','KPA EC Sublimated Mass','KPA Bulk Sublimated Mass','location','northwest');
grid on;
set(gca,'xtick',datenum([2021 12 1 0 0 0]):10:datenum([2022 5 15 0 0 0]));
xlim([datenum([2021 12 1 0 0 0]) datenum([2022 5 15 0 0 0])]); ylim([0 50]);
datetick('x','mmmdd','keepticks','keeplimits');
set(gca,'fontsize',12);
subplot(3,3,6)
plot(barr.dn_wy22,sub_mass_ayp22_fill_ec,'r','linewidth',1); % this is what was sublimated
hold on;
plot(barr.dn_wy22,sub_mass_ayp22_fill_bulk,'r:','linewidth',1); % this is what was sublimated
ylabel('Daily Sublimation [kg/m^2]');
legend('AYP EC Sublimated Mass','AYP Bulk Sublimated Mass','location','northwest');
ttl=title('h)');ttl.Units = 'Normalize';ttl.Position(1) = 0;ttl.HorizontalAlignment = 'left';
grid on;
set(gca,'xtick',datenum([2021 12 1 0 0 0]):10:datenum([2022 5 15 0 0 0]));
xlim([datenum([2021 12 1 0 0 0]) datenum([2022 5 15 0 0 0])]); ylim([0 2.5]);
datetick('x','mmmdd','keepticks','keeplimits');
set(gca,'fontsize',12);
subplot(3,3,9)
plot(barr.dn_wy22,sub_mass_kpa22_fill_ec,'b','linewidth',1); % this is what was sublimated
hold on;
plot(barr.dn_wy22,sub_mass_kpa22_fill_bulk,'b:','linewidth',1); % this is what was sublimated
ylabel('Daily Sublimation [kg/m^2]');
legend('KPA EC Sublimated Mass','KPA Bulk Sublimated Mass','location','northwest');
ttl=title('i)');ttl.Units = 'Normalize';ttl.Position(1) = 0;ttl.HorizontalAlignment = 'left';
grid on;
set(gca,'xtick',datenum([2021 12 1 0 0 0]):10:datenum([2022 5 15 0 0 0]));
xlim([datenum([2021 12 1 0 0 0]) datenum([2022 5 15 0 0 0])]); ylim([0 2.5]);
datetick('x','mmmdd','keepticks','keeplimits');
set(gca,'fontsize',12);

set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0 0 45 30]);
print -dpng -r300 /Users/ccox/Documents/Projects/2021/splash/science/sublimation/pub/figs/Figure_3.png




[sub_mass_ec_ayp22,sub_mass_bulk_ayp22,sub_mass_be_ayp22] = calc_sub_precip(a50.temp,a50.skin_temp_surface,hl_50,a50.Hl,a50.snow_flag,barr.dn_wy22,a50.dn,'');
[sub_mass_ec_ayp23,sub_mass_bulk_ayp23,sub_mass_be_ayp23] = calc_sub_precip(a50.temp,a50.skin_temp_surface,hl_50,a50.Hl,a50.snow_flag,barr.dn_wy23,a50.dn,'');

[sub_mass_ec_kpa22,sub_mass_bulk_kpa22,sub_mass_be_kpa22] = calc_sub_precip(a30.temp,a30.skin_temp_surface,hl_30,a30.Hl,a30.snow_flag,barr.dn_wy22,a30.dn,'');
[sub_mass_ec_kpa23,sub_mass_bulk_kpa23,sub_mass_be_kpa23] = calc_sub_precip(a30.temp,a30.skin_temp_surface,hl_30,a30.Hl,a30.snow_flag,barr.dn_wy23,a30.dn,'');

sub_mass_ayp22_fill = sub_mass_be_ayp22;
sub_mass_ayp22_fill(isnan(sub_mass_ayp22_fill)) = sub_mass_be_kpa22(isnan(sub_mass_ayp22_fill));

sub_mass_kpa22_fill = sub_mass_be_kpa22;
sub_mass_kpa22_fill(isnan(sub_mass_kpa22_fill)) = sub_mass_be_ayp22(isnan(sub_mass_kpa22_fill));

sub_mass_ayp23_fill = sub_mass_be_ayp23;
sub_mass_ayp23_fill(isnan(sub_mass_ayp23_fill)) = sub_mass_be_kpa23(isnan(sub_mass_ayp23_fill));

sub_mass_kpa23_fill = sub_mass_be_kpa23;
sub_mass_kpa23_fill(isnan(sub_mass_kpa23_fill)) = sub_mass_be_ayp23(isnan(sub_mass_kpa23_fill));

sub_mass_kpa22_fill(barr.dn_wy22 > datenum([2022 5 4 0 0 0])) = NaN;    
sub_mass_kpa23_fill(barr.dn_wy23 > datenum([2023 5 16 0 0 0])) = NaN; 

% I'm setting a period in Nov 2022 to 0 sublimation for when south facing slopes became snow free
sub_mass_kpa23_fill_orig = sub_mass_kpa23_fill;
sub_mass_kpa23_fill(barr.dn_wy23 >= datenum([2022 11 13 0 0 0]) & barr.dn_wy23 <= datenum([2022 11 28 0 0 0])) = 0;

clf;
subplot(2,1,1)
plot(barr.dn_wy22,nancumsum(barr.armpluv_new_swe_wy22.*1000),'k-','color',[1 0.7 0],'linewidth',2); % this is what went in
hold on;
plot(barr.dn_wy22,cumsum(sub_mass_ayp22_fill),'r','linewidth',2); % this is what was sublimated
plot(barr.dn_wy22,cumsum(sub_mass_kpa22_fill),'b','linewidth',2); % this is what was sublimated
plot(barr.dn_wy22,(cumsum(sub_mass_kpa22_fill)+cumsum(sub_mass_ayp22_fill))./2,'k--','linewidth',2);
ylabel('Sublimated SWE [kg/m^2]');
legend('SAIL Precip Mass PUVIO','AYP BE Sublimated Mass','KPA BE Sublimated Mass','Mean Sublimated Mass','location','northwest');
grid on;
set(gca,'xtick',datenum([2021 10 1 0 0 0]):10:datenum([2022 6 1 0 0 0]));
xlim([datenum([2021 10 1 0 0 0]) datenum([2022 6 1 0 0 0])]); ylim([0 700]);
datetick('x','mmmdd','keepticks','keeplimits');
title('2021-2022');
set(gca,'fontsize',12);
subplot(2,1,2)
plot(barr.dn_wy23,cumsum(barr.armpluv_new_swe_wy23.*1000),'k-','color',[1 0.7 0],'linewidth',2); % this is what went in
hold on;
plot(barr.dn_wy23,cumsum(sub_mass_ayp23_fill),'r','linewidth',2); % this is what was sublimated
plot(barr.dn_wy23,cumsum(sub_mass_kpa23_fill),'b','linewidth',2); % this is what was sublimated
plot(barr.dn_wy23,(cumsum(sub_mass_kpa23_fill)+cumsum(sub_mass_ayp23_fill))./2,'k--','linewidth',2);
plot(barr.dn_wy23,cumsum(sub_mass_kpa23_fill_orig),'color',[0.5 0.5 0.5],'linewidth',1); % this is what was sublimated is I include late Nov at KPA
ylabel('Sublimated SWE [kg/m^2]');
legend('SAIL Precip Mass PUVIO','AYP BE Sublimated Mass','KPA BE Sublimated Mass','Mean Sublimated Mass','KPA all','location','northwest');
grid on;
set(gca,'xtick',datenum([2022 10 1 0 0 0]):10:datenum([2023 6 1 0 0 0]));
xlim([datenum([2022 10 1 0 0 0]) datenum([2023 6 1 0 0 0])]); ylim([0 700]);
datetick('x','mmmdd','keepticks','keeplimits');
title('2022-2023');
set(gca,'fontsize',12);
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0 0 23 34]);
print -dpng -r300 /Users/ccox/Documents/Projects/2021/splash/science/sublimation/pub/figs/Figre_4/Figure_4ab.png


clf;
subplot(2,1,1)
plot(barr.dn_wy22,cumsum(sub_mass_ayp22_fill),'r','linewidth',2); % this is what was sublimated
hold on;
plot(barr.dn_wy22,cumsum(sub_mass_kpa22_fill),'b','linewidth',2); % this is what was sublimated
plot(barr.dn_wy22,(cumsum(sub_mass_kpa22_fill)+cumsum(sub_mass_ayp22_fill))./2,'k--','linewidth',2);
ylabel('Sublimated SWE [kg/m^2]');
legend('AYP BE Sublimated Mass','KPA BE Sublimated Mass','Mean Sublimated Mass','location','northwest');
grid on;
set(gca,'xtick',datenum([2021 10 1 0 0 0]):10:datenum([2022 6 1 0 0 0]));
xlim([datenum([2021 10 1 0 0 0]) datenum([2022 6 1 0 0 0])]); ylim([0 70]);
datetick('x','mmmdd','keepticks','keeplimits');
title('2021-2022');
set(gca,'fontsize',12);
subplot(2,1,2)
plot(barr.dn_wy23,cumsum(sub_mass_ayp23_fill),'r','linewidth',2); % this is what was sublimated
hold on;
plot(barr.dn_wy23,cumsum(sub_mass_kpa23_fill),'b','linewidth',2); % this is what was sublimated
plot(barr.dn_wy23,(cumsum(sub_mass_kpa23_fill)+cumsum(sub_mass_ayp23_fill))./2,'k--','linewidth',2);
plot(barr.dn_wy23,cumsum(sub_mass_kpa23_fill_orig),'color',[0.5 0.5 0.5],'linewidth',1); % this is what was sublimated is I include late Nov at KPA
ylabel('Sublimated SWE [kg/m^2]');
legend('AYP BE Sublimated Mass','KPA BE Sublimated Mass','Mean Sublimated Mass','KPA all','location','northwest');
grid on;
set(gca,'xtick',datenum([2022 10 1 0 0 0]):10:datenum([2023 6 1 0 0 0]));
xlim([datenum([2022 10 1 0 0 0]) datenum([2023 6 1 0 0 0])]); ylim([0 70]);
datetick('x','mmmdd','keepticks','keeplimits');
title('2022-2023');
set(gca,'fontsize',12);
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0 0 23 34]);
print -dpng -r300 /Users/ccox/Documents/Projects/2021/splash/science/sublimation/pub/figs/Figure_4/Figure_4cd.png

% peak snow depth
% March 10, 2022
% March 22, 2023


[sub_mass_ec_ayp22c,sub_mass_bulk_ayp22c,sub_mass_be_ayp22c] = calc_sub_precip(a50.temp,a50.skin_temp_surface,hl_50,a50.Hl,a50.snow_flag,barr.dn_wy22,a50.dn,'cold');
[sub_mass_ec_kpa22c,sub_mass_bulk_kpa22c,sub_mass_be_kpa22c] = calc_sub_precip(a30.temp,a30.skin_temp_surface,hl_30,a30.Hl,a30.snow_flag,barr.dn_wy22,a30.dn,'cold');

[sub_mass_ec_ayp22w,sub_mass_bulk_ayp22w,sub_mass_be_ayp22w] = calc_sub_precip(a50.temp,a50.skin_temp_surface,hl_50,a50.Hl,a50.snow_flag,barr.dn_wy22,a50.dn,'warm');
[sub_mass_ec_kpa22w,sub_mass_bulk_kpa22w,sub_mass_be_kpa22w] = calc_sub_precip(a30.temp,a30.skin_temp_surface,hl_30,a30.Hl,a30.snow_flag,barr.dn_wy22,a30.dn,'warm');

sub_mass_ayp22_fillc = sub_mass_be_ayp22c;
sub_mass_ayp22_fillc(isnan(sub_mass_ayp22_fillc)) = sub_mass_be_kpa22c(isnan(sub_mass_ayp22_fillc));

sub_mass_kpa22_fillc = sub_mass_be_kpa22c;
sub_mass_kpa22_fillc(isnan(sub_mass_kpa22_fillc)) = sub_mass_be_ayp22c(isnan(sub_mass_kpa22_fillc));

sub_mass_ayp22_fillw = sub_mass_be_ayp22w;
sub_mass_ayp22_fillw(isnan(sub_mass_ayp22_fillw)) = sub_mass_be_kpa22w(isnan(sub_mass_ayp22_fillw));

sub_mass_kpa22_fillw = sub_mass_be_kpa22w;
sub_mass_kpa22_fillw(isnan(sub_mass_kpa22_fillw)) = sub_mass_be_ayp22w(isnan(sub_mass_kpa22_fillw));

sub_mass_ayp22_fillc(isnan(sub_mass_ayp22_fillc)) = 0;
sub_mass_ayp22_fillw(isnan(sub_mass_ayp22_fillw)) = 0;
sub_mass_kpa22_fillc(isnan(sub_mass_kpa22_fillc)) = 0;
sub_mass_kpa22_fillw(isnan(sub_mass_kpa22_fillw)) = 0;

sub_mass_kpa22_fillc(barr.dn_wy22 > datenum([2022 5 4 0 0 0])) = NaN;    
sub_mass_kpa22_fillw(barr.dn_wy22 > datenum([2022 5 4 0 0 0])) = NaN;    




[sub_mass_ec_ayp23c,sub_mass_bulk_ayp23c,sub_mass_be_ayp23c] = calc_sub_precip(a50.temp,a50.skin_temp_surface,hl_50,a50.Hl,a50.snow_flag,barr.dn_wy23,a50.dn,'cold');
[sub_mass_ec_kpa23c,sub_mass_bulk_kpa23c,sub_mass_be_kpa23c] = calc_sub_precip(a30.temp,a30.skin_temp_surface,hl_30,a30.Hl,a30.snow_flag,barr.dn_wy23,a30.dn,'cold');

[sub_mass_ec_ayp23w,sub_mass_bulk_ayp23w,sub_mass_be_ayp23w] = calc_sub_precip(a50.temp,a50.skin_temp_surface,hl_50,a50.Hl,a50.snow_flag,barr.dn_wy23,a50.dn,'warm');
[sub_mass_ec_kpa23w,sub_mass_bulk_kpa23w,sub_mass_be_kpa23w] = calc_sub_precip(a30.temp,a30.skin_temp_surface,hl_30,a30.Hl,a30.snow_flag,barr.dn_wy23,a30.dn,'warm');
sub_mass_be_kpa23w(barr.dn_wy23 >= datenum([2022 11 13 0 0 0]) & barr.dn_wy23 <= datenum([2022 11 28 0 0 0])) = 0;
sub_mass_be_kpa23c(barr.dn_wy23 >= datenum([2022 11 13 0 0 0]) & barr.dn_wy23 <= datenum([2022 11 28 0 0 0])) = 0;

sub_mass_ayp23_fillc = sub_mass_be_ayp23c;
sub_mass_ayp23_fillc(isnan(sub_mass_ayp23_fillc)) = sub_mass_be_kpa23c(isnan(sub_mass_ayp23_fillc));

sub_mass_kpa23_fillc = sub_mass_be_kpa23c;
sub_mass_kpa23_fillc(isnan(sub_mass_kpa23_fillc)) = sub_mass_be_ayp23c(isnan(sub_mass_kpa23_fillc));

sub_mass_ayp23_fillw = sub_mass_be_ayp23w;
sub_mass_ayp23_fillw(isnan(sub_mass_ayp23_fillw)) = sub_mass_be_kpa23w(isnan(sub_mass_ayp23_fillw));

sub_mass_kpa23_fillw = sub_mass_be_kpa23w;
sub_mass_kpa23_fillw(isnan(sub_mass_kpa23_fillw)) = sub_mass_be_ayp23w(isnan(sub_mass_kpa23_fillw));

sub_mass_ayp23_fillc(isnan(sub_mass_ayp23_fillc)) = 0;
sub_mass_ayp23_fillw(isnan(sub_mass_ayp23_fillw)) = 0;
sub_mass_kpa23_fillc(isnan(sub_mass_kpa23_fillc)) = 0;
sub_mass_kpa23_fillw(isnan(sub_mass_kpa23_fillw)) = 0;


sub_mass_kpa23_fillc(barr.dn_wy23 > datenum([2023 5 16 0 0 0])) = NaN;    
sub_mass_kpa23_fillw(barr.dn_wy23 > datenum([2023 5 16 0 0 0])) = NaN;    


clf;
subplot(2,1,1)
plot(barr.dn_wy22,cumsum(sub_mass_ayp22_fillc),'r--','linewidth',2); % this is what was sublimated
hold on;
plot(barr.dn_wy22,cumsum(sub_mass_ayp22_fillw),'r','linewidth',2); % this is what was sublimated
plot(barr.dn_wy22,cumsum(sub_mass_kpa22_fillc),'b--','linewidth',2); % this is what was sublimated
plot(barr.dn_wy22,cumsum(sub_mass_kpa22_fillw),'b','linewidth',2); % this is what was sublimated
ylabel('Sublimated SWE [kg/m^2]');
legend('AYP BE Sublimated Mass (cold regime)','AYP BE Sublimated Mass (warm regime)', ...
       'KPA BE Sublimated Mass (cold regime)','KPA Mean Sublimated Mass (warm regime)','location','northwest');
grid on;
set(gca,'xtick',datenum([2021 10 1 0 0 0]):10:datenum([2022 6 1 0 0 0]));
xlim([datenum([2021 10 1 0 0 0]) datenum([2022 6 1 0 0 0])]); ylim([0 70]);
datetick('x','mmmdd','keepticks','keeplimits');
title('2021-2022');
set(gca,'fontsize',12);
subplot(2,1,2)
plot(barr.dn_wy23,cumsum(sub_mass_ayp23_fillc),'r--','linewidth',2); % this is what was sublimated
hold on;
plot(barr.dn_wy23,cumsum(sub_mass_ayp23_fillw),'r','linewidth',2); % this is what was sublimated
plot(barr.dn_wy23,cumsum(sub_mass_kpa23_fillc),'b--','linewidth',2); % this is what was sublimated
plot(barr.dn_wy23,cumsum(sub_mass_kpa23_fillw),'b','linewidth',2); % this is what was sublimated
ylabel('Sublimated SWE [kg/m^2]');
legend('AYP BE Sublimated Mass (cold regime)','AYP BE Sublimated Mass (warm regime)', ...
       'KPA BE Sublimated Mass (cold regime)','KPA Mean Sublimated Mass (warm regime)','location','northwest');
grid on;
set(gca,'xtick',datenum([2022 10 1 0 0 0]):10:datenum([2023 6 1 0 0 0]));
xlim([datenum([2022 10 1 0 0 0]) datenum([2023 6 1 0 0 0])]); ylim([0 70]);
datetick('x','mmmdd','keepticks','keeplimits');
title('2022-2023');
set(gca,'fontsize',12);
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0 0 23 34]);
print -dpng -r300 /Users/ccox/Documents/Projects/2021/splash/science/sublimation/pub/figs/Figure_4/Figure_4ef.png


%% Figure 5 


Hl = a50.Hl; 
Hs = a50.Hs;
ind1 = find(a50.temp < 0 & a50.snow_flag == 1 & a50snow_depth > 0.2);
ind2 = find(a50.temp > 0 & a50.snow_flag == 1 & a50snow_depth > 0.2);
ind3 = find(a50.dn >= datenum([2023 4 18 0 0 0]) & a50.dn <= datenum([2023 4 20 0 0 0]) & a50.snow_flag == 1 & a50.temp > 0);

cmap = lines;
clf;
h1=plot(Hs(ind1),Hl(ind1),'.','color',[0.75 0.75 0.75],'markersize',12); hold on;
xlim([-125 125]); ylim([-125 125]);
grid on; hline(0); vline(0);
plot(-1.*[-125:125],-125:125,'k');
legend(h1,'Cold Regime', ...
    'location','southwest');
set(gca,'xtick',[-120:20:120],'ytick',[-120:20:120]);
xlabel('Sensible Heat Flux ({\it H_s}) [Wm^-^2]');
ylabel('Latent Heat Flux ({\it H_l}) [Wm^-^2]');
set(gca,'fontsize',16); 
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0 0 22 20]);
print -dpng -r300 /Users/ccox/Documents/Projects/2021/splash/science/sublimation/pub/figs/Figure_5/Figure_5a.png

clf;
h1=plot(Hs(ind1),Hl(ind1),'.','color',[0.75 0.75 0.75],'markersize',12); hold on;
h2=plot(Hs(ind2),Hl(ind2),'.','color',cmap(1,:),'markersize',12);
h3=plot(Hs(ind3),Hl(ind3),'.','color',cmap(3,:),'markersize',12);
xlim([-125 125]); ylim([-125 125]);
grid on; hline(0); vline(0);
plot(-1.*[-125:125],-125:125,'k');
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
print -dpng -r300 /Users/ccox/Documents/Projects/2021/splash/science/sublimation/pub/figs/Figure_5/Figure_5b.png

ind2Fig12 = find(a50.temp > 0 & a50.snow_flag == 1 & a50snow_depth > 0.2 & abs(Hs) > 5);
nancorrcoef(Hs(ind2Fig12),Hl(ind2Fig12))
nanmean(Hs(ind2Fig12)+Hl(ind2Fig12))
nanmean(Hs(ind2Fig12)) / nanmean(Hl(ind2Fig12))




Hl = a50.bulk_Hl; 
Hs = a50.bulk_Hs;
ind1 = find(a50.temp < 0 & a50.snow_flag == 1 & a50snow_depth > 0.2);
ind2 = find(a50.temp > 0 & a50.snow_flag == 1 & a50snow_depth > 0.2);
ind3 = find(a50.dn >= datenum([2023 4 18 0 0 0]) & a50.dn <= datenum([2023 4 20 0 0 0]) & a50.snow_flag == 1 & a50.temp > 0);

clf;
h1=plot(Hs(ind1),Hl(ind1),'.','color',[0.75 0.75 0.75],'markersize',12); hold on;
h2=plot(Hs(ind2),Hl(ind2),'.','color',cmap(1,:),'markersize',12);
h3=plot(Hs(ind3),Hl(ind3),'.','color',cmap(3,:),'markersize',12);
xlim([-125 125]); ylim([-125 125]);
grid on; hline(0); vline(0);
plot(-1.*[-125:125],-125:125,'k');
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
print -dpng -r300 /Users/ccox/Documents/Projects/2021/splash/science/sublimation/pub/figs/Figure_12/Figure_12b.png

ind2Fig12 = find(a50.temp > 0 & a50.snow_flag == 1 & a50snow_depth > 0.2 & abs(Hs) > 5);
nancorrcoef(Hs(ind2Fig12),Hl(ind2Fig12))
nanmean(Hs(ind2Fig12)+Hl(ind2Fig12))
nanmean(Hs(ind2Fig12)) / nanmean(Hl(ind2Fig12))

%% Figure 6b 


dn2 = a50.dn - 7/24;
ind = find(dn2 >= datenum([2023 4 16 0 0 0]) & dn2 <= datenum([2023 4 21 0 0 0]));

Hlf6 = a50.Hl(ind); Hlf6(Hlf6<-20|Hlf6>98)=NaN;
clf;
plot(dn2(ind),a50.bulk_Hl(ind),'linewidth',1); hold on;
plot(dn2(ind),Hlf6,'linewidth',2); 
plot(dn2(ind),a50.temp(ind),'linewidth',2);
hline(0); grid on; ylim([-20 100]);
set(gca,'xtick',[datenum([2023 4 16 0 0 0]):1:datenum([2023 4 21 0 0 0])])
datetick('x','mm/dd','keepticks','keeplimits');
ylabel('H_l [Wm^-^2], Temp. [C]'); 
legend('H_l (bulk)','H_l EC','Air Temp.');
xlim([datenum([2023 4 16 0 0 0]) datenum([2023 4 21 0 0 0])]);
xlabel('Time [MST]')
set(gca,'fontsize',24);
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0 0 40 14]);
print -dpng -r300 /Users/ccox/Documents/Projects/2021/splash/science/sublimation/pub/figs/Figure_6/Figure_6b.png


%% Figure 13

clf; 
yyaxis left
scatter(dp,rats_l0,n_l0/5,hl_l0); hold on;
scatter(dp,rats_g0,n_g0/5,hl_g0,'filled'); grid on; 
xlim([-30 5]); ylim([-3 6]); box on;
h=colorbar; h.Label.String = 'H_l [Wm^-^2]'; colormap(jet); caxis([0 30]); grid on;
set(gcf,'PaperPositionMode','manual'); hline(-1); vline(0); 
ylabel('Bowen Ratio (H_s / H_l)'); xlabel('T_d [C]');
yyaxis right
ind_all = find(a50.temp > 0 & a50.snow_flag == 1 & a50snow_depth > 0.2 & abs(Hs) > 1 & abs(Hl) > 1);
histN(a50.dew_point(ind_all),-30:0.5:5,'line','color',cmap(2,:));
ylabel('%');
ax = gca;
ax.YAxis(2).Color = cmap(2,:);
ax.YAxis(1).Color = 'k';
set(gca,'fontsize',14)
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0 0 20 16]);
print -dpng -r300 /Users/ccox/Documents/Projects/2021/splash/science/sublimation/pub/figs/Figure_13.png