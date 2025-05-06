

addpath '/Users/ccox/Documents/Projects/2021/splash/Simba_NOAA0101/data_processing/sledd'
cd /Users/ccox/Documents/Projects/2021/splash/science/sublimation/sos_case_study/
load('/Users/ccox/Documents/Projects/2021/splash/_level2/a50dbl.mat');
a50 = datam;
load('/Users/ccox/Documents/Projects/2021/splash/_level2/a30dbl.mat')
a30 = datam;

% Process the tower -----------
s = rd_netcdf('/Users/ccox/Documents/Projects/2021/splash/science/sublimation/sos_case_study/sos_data/isfs_sos_qc_geo_tiltcor_hr_20230418_23.nc');
s.dn = datenum([1970 1 1 0 0 0]) + double(s.base_time)/86400 + s.time./86400;
si = find(s.dn >= datenum([2023 4 18 23 24 0]) & s.dn < datenum([2023 4 18 23 34 0]));

s.P_10m_c(s.P_10m_c > 10^5) = NaN;
s.P_20m_c(s.P_20m_c > 10^5) = NaN;
s.P_10m_c = nanmean(s.P_10m_c);
s.P_20m_c = nanmean(s.P_20m_c);

pr = repmat(1:20,3600,1).*NaN;
for k = 1:length(s.P_10m_c)
    pr(k,:) = interp1([10 20],[s.P_10m_c(:,k);s.P_20m_c(:,k)],1:20,'linear','extrap');
end

t = [s.T_1m_c s.T_2m_c s.T_3m_c s.T_4m_c ...
    s.T_5m_c s.T_6m_c s.T_7m_c s.T_8m_c ...
    s.T_9m_c s.T_10m_c s.T_11m_c s.T_12m_c ...
    s.T_13m_c s.T_14m_c s.T_15m_c s.T_16m_c ...
    s.T_17m_c s.T_18m_c s.T_19m_c s.T_20m_c];
t(t > 10^5) = NaN;


rh = [s.RH_1m_c s.RH_2m_c s.RH_3m_c s.RH_4m_c ...
    s.RH_5m_c s.RH_6m_c s.RH_7m_c s.RH_8m_c ...
    s.RH_9m_c s.RH_10m_c s.RH_11m_c s.RH_12m_c ...
    s.RH_13m_c s.RH_14m_c s.RH_15m_c s.RH_16m_c ...
    s.RH_17m_c s.RH_18m_c s.RH_19m_c s.RH_20m_c];
rh(rh > 10^5) = NaN;

PREF = nanmean(pr(si,1));
twr_theta = ((t+273.15).* (PREF./pr).^0.286)-273.15;

% -----------------------------


% Process Avery met -----------
% the height of this measurements is 1.32 m

ind = find(a50.dn >= datenum([2023 4 18 23 24 0]) & a50.dn < datenum([2023 4 18 23 34 0]));
temp_vaisala_50 = a50.temp(ind);
theta_vaisala_50 = ((temp_vaisala_50+273.15).*(PREF./a50.atmos_pressure(ind)).^0.286)-273.15;

ind = find(a30.dn >= datenum([2023 4 18 23 24 0]) & a30.dn < datenum([2023 4 18 23 34 0]));
temp_vaisala_30 = a30.temp(ind);
%theta_vaisala_30 = ((temp_vaisala_30+273.15).*((a30.atmos_pressure(ind)+0.01)./a30.atmos_pressure(ind)).^0.286)-273.15;
theta_vaisala_30 = ((temp_vaisala_30+273.15).*(PREF./a30.atmos_pressure(ind)).^0.286)-273.15;

% -----------------------------


% Process sonde ---------------

x = rd_netcdf('/Users/ccox/Documents/Projects/2021/splash/science/sublimation/sos_case_study/data_radiosonde/gucsondewnpnM1.b1.20230418.232900.cdf');
s_ht = x.alt-x.alt(1); s_ht(1) = 2; % make first level 2 m gnd station
%s_theta = ((x.tdry+273.15).*(x.pres(1)./x.pres).^0.286)-273.15;
s_theta = ((x.tdry+273.15).*(PREF./x.pres).^0.286)-273.15;
wvout1 = humidRS80(x.rh,x.tdry+273.15,1,9,x.pres);

% -----------------------------


% Process SIMBA ---------------

ts=get_avery_raw; cd /Users/ccox/Documents/Projects/2021/splash/science/sublimation/sos_case_study/

ind = find(floor(ts.dn) == datenum([2023 4 18 0 0 0]));
rind = find(ts.range > 1.18);
% correct the in air sensors
t_vaisala = interp1(a50.dn,a50.temp,ts.dn);
t_sfc = interp1(a50.dn,a50.skin_temp_surface,ts.dn);
sza = interp1(datam.dn,datam.zenith_true,ts.dn);
simbt = mean(ts.t(1:10,:))'; simbt(7672:7739) = NaN;
[yyyy mm dd HH u u] = datevec(ts.dn);
i_tcor = find(ts.dn > datenum([2023 1 1 0 0 0]) & ts.dn < datenum([2023 4 24 0 0 0]) & simbt-t_vaisala < 6 & simbt-t_vaisala > -4);
pfits = polyfit(sza(i_tcor),simbt(i_tcor)-t_vaisala(i_tcor),2);
tcor = polyval(pfits,sza);
%ts.t(1:24,:) = ts.t(1:24,:) - tcor(ind(4));
pr_simb = data10.atmos_pressure(findnearest(data10.dn,datenum([2023 4 18 0 0 0])));
solar_off = (ts.t(1,ind(4))-mean(temp_vaisala_50));
%ts.t(1:24,:) = ts.t(1:24,:) - solar_off;

rh_str = ts.range.*0 + 100;
rh_str(rind) = interp1([0 max(ts.range(rind)-ts.range(rind(end)))],[100 nanmean(rh(si,1))],ts.range(rind)-ts.range(rind(end)));
%ts_theta = ((ts.t(:,ind(end))+273.15).*(pr_simb./pr_simb).^0.286)-273.15;
ts_theta = (ts.t(:,ind(end))+273.15);
ts_theta(rind) = (ts_theta(rind) .* (PREF./pr_simb).^0.286);
ts_theta = ts_theta - 273.15;

wvout3 = humidRS80(ts.t(:,ind(end)).*0+100,ts.t(:,ind(end))+273.15,1,9,ts.t(:,ind(end)).*0+x.pres(1));

% Snow
% Avery: 112 cm (154.2 cm at pinger) ts.range(find(ts.flag(:,ind(end)) == 3,1,'last'))
sdoff = 1-0.84+(1.12-0.84); % offset to account for snow depth at twr (first term) and diff in snow depth between tower and Avery

% -----------------------------



% Get cloud height ------------

c = rd_netcdf('/Users/ccox/Documents/Projects/2021/splash/science/sublimation/sos_case_study/data_radiosonde/gucceil10mM1.b1.20230418.000011.nc');
c.dn = datenum([1970 1 1 0 0 0]) + double(c.base_time)/86400 + c.time_offset./86400;
ic = findnearest(c.dn,datenum([2023 4 18 23 29 0]));

% -----------------------------


cmap = lines;
clf;

ti = tiledlayout(1,1);
ax1 = axes(ti);

toff = nanmean(theta_vaisala_30)-nanmean(theta_vaisala_50);

h1=plot(s_theta,s_ht,'linewidth',2,'color',cmap(3,:)); hold on;
h2=plot(nanmean(twr_theta(si,:)),[1:20]+sdoff,'k.-','linewidth',2,'color',cmap(6,:),'markersize',14);
h3=plot(ts_theta,ts.range,'linewidth',2,'color',cmap(2,:));

plot(nanmean(theta_vaisala_30),1.32+1.18,'ko','color',cmap(1,:),'markersize',10,'markerfacecolor',cmap(1,:));
%plot(nanmean(theta_vaisala_50),1.32+1.18,'ko','color',cmap(2,:),'markersize',10,'markerfacecolor',cmap(2,:));

set(gca,'yscale','log');
xlim([-1 10]); ylim([0.02 5000]);

ax=get(gca);
curtick = ax.YTick;
set(gca, 'YTickLabel', cellstr(num2str(curtick(:)))); 
xlabel('{\it \theta} [C]'); ylabel('Height above soil [m]');
hline(ts.range(rind(end)))
hline(c.first_cbh(ic))
vline(0)
hline(1000)
grid on; 
set(gca,'xtick',-2:11,'fontsize',22)

set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0 0 30 30]);
print('-dpng','-r300','/Users/ccox/Documents/Projects/2021/splash/science/sublimation/pub/figs/Figure_7/Figure_7orig.png');




% 
% 
% 
% cmap = lines;
% clf;
% 
% h1=plot(s_theta,x.alt-x.alt(1)+3,'linewidth',2,'color',cmap(5,:)); hold on;
% h2=plot(nanmean(twr_theta(si,:)),[1:20]+sdoff,'linewidth',2,'color',cmap(1,:));
% h3=plot(ts_theta+solar_off,ts.range,'linewidth',2,'color',cmap(2,:));
% 
% plot(nanmean(theta_vaisala_50),1.32+1.18,'ko','color',cmap(2,:),'markersize',10,'markerfacecolor',cmap(2,:));
% plot(nanmean(theta_vaisala_30),1.32+1.18,'ko','color',cmap(1,:),'markersize',10,'markerfacecolor',cmap(1,:));
% 
% set(gca,'yscale','log');
% xlim([-1 10]); ylim([0.02 5000]);
% 
% ax=get(gca);
% curtick = ax.YTick;
% set(gca, 'YTickLabel', cellstr(num2str(curtick(:)))); 
% xlabel('{\it \theta} [C]'); ylabel('Height above soil [m]');
% hline(ts.range(rind(end)))
% hline(c.first_cbh(ic))
% vline(0)
% hline(1000)
% grid on; 
% set(gca,'xtick',-2:11,'fontsize',22)
% 
% 
% 
% 
% 
% cmap = lines;
% clf;
% 
% ti = tiledlayout(1,1);
% ax1 = axes(ti);
% 
% toff = nanmean(theta_vaisala_30)-nanmean(theta_vaisala_50);
% 
% h1=plot(s_theta,s_ht,'linewidth',2,'color',cmap(2,:)); hold on;
% h2=plot(nanmean(twr_theta(si,:))-toff,[1:20]+sdoff,'linewidth',2,'color',cmap(2,:));
% h2a=plot(nanmean(twr_theta(si,:)),[1:20]+sdoff,'k--','linewidth',1,'color',cmap(1,:));
% h3=plot(ts_theta,ts.range,'linewidth',2,'color',cmap(2,:));
% 
% plot(nanmean(theta_vaisala_30),1.32+1.18,'ko','color',cmap(1,:),'markersize',10,'markerfacecolor',cmap(1,:));
% plot(nanmean(theta_vaisala_50),1.32+1.18,'ko','color',cmap(2,:),'markersize',10,'markerfacecolor',cmap(2,:));
% 
% set(gca,'yscale','log');
% xlim([-1 10]); ylim([0.02 5000]);
% 
% ax=get(gca);
% curtick = ax.YTick;
% set(gca, 'YTickLabel', cellstr(num2str(curtick(:)))); 
% xlabel('{\it \theta} [C]'); ylabel('Height above soil [m]');
% hline(ts.range(rind(end)))
% hline(c.first_cbh(ic))
% vline(0)
% hline(1000)
% grid on; 
% set(gca,'xtick',-2:11,'fontsize',22)
% ax1.XColor = cmap(2,:);
% 
% 
% ax2 = axes(ti);
% 
% 
% wvout2 = humidRS80(rh,t+273.15,1,9,s.P_10m_c);
% rh2 = humidRS80(wvout2,t-toff+273.15,8,4,s.P_10m_c);
% wvout2 = humidRS80(rh2,t-toff+273.15,1,9,s.P_10m_c);
% 
% q_sonde = wvout1 ./ (1+wvout1);
% h4=plot(q_sonde,x.alt-x.alt(1)+3,'linewidth',2,'color',cmap(5,:)); ylim([0 5000]); % +2 bc "surface" will be 2 m above 1 m of snow
% hold on;
% sdoff = ts.range(find(ts.flag(:,ind(end)) == 3,1,'last'))-0.84; % offset to account for diff in snow depth between tower and Avery (second term)
% 
% q = wvout2 ./ (1+wvout2);
% h5=plot(nanmean(q(si,:)),[1:20]+sdoff,'linewidth',2,'color',cmap(5,:));
% 
% q_sfc = wvout3 ./ (1+wvout3);
% q_sfc(1:find(ts.flag(:,ind(end)) == 3,1,'last')) = NaN;
% h6=plot(q_sfc,ts.range,'linewidth',2,'color',cmap(5,:));
% 
% set(gca,'yscale','log');
% hline(ts.range(find(ts.flag(:,ind(end)) == 3,1,'last')))
% ax=get(gca);
% curtick = ax.YTick;
% set(gca, 'YTickLabel', cellstr(num2str(curtick(:)))); 
% xlabel('{\it q} [g/kg]'); ylabel('Height above soil [m]');
% hline(c.first_cbh(ic))
% plot([0.5 0.9],[1000 1000],'k--','linewidth',1); 
% plot([0 0],[10^-2 5000],'k--','linewidth',2); 
% xlim([0.5 0.9]); ylim([0.02 5000]);
% grid on; 
% legend([h1 h4],'Potential Temperature (\theta)','Specific Humidty','location','south');
% set(gca,'xtick',0.5:0.1:0.9);
% set(gca,'fontsize',22);
% 
% ax2.XAxisLocation = 'top';
% ax2.YAxisLocation = 'right';
% ax2.Color = 'none';
% ax2.XColor = cmap(5,:);
% ax2.YColor = 'none';
% ax1.Box = 'on';
% ax2.Box = 'on';
% ax1.XTickLabelRotation = 0;
% ax2.XTickLabelRotation = 0;
% 
% set(gcf,'PaperPositionMode','manual');
% set(gcf,'PaperUnits','centimeters');
% set(gcf,'PaperPosition',[0 0 30 30]);
% print('-dpng','-r300','prof_v2024aug13.png');
% 
% 
% 
