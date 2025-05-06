days = [15:31 1:30];
cmapr = flipud(cbrewer('seq', 'Reds',length(days), 'linear'));
cmapb = flipud(cbrewer('seq', 'Blues',length(days), 'linear'));
Ta = [];
dT = [];
for j = 1:length(days)
    j
    theday = num2str(days(j));
    if length(theday) == 1; theday = ['0',theday]; end
    if j <= 17
        file = rd_netcdf(['isfs_sos_qc_geo_tiltcor_hr_202303',theday,'_23.nc'],'T_10m_c','T_3m_c');
    else
        file = rd_netcdf(['isfs_sos_qc_geo_tiltcor_hr_202304',theday,'_23.nc'],'T_10m_c','T_3m_c');
    end
    file.T_10m_c(file.T_10m_c > 1e6) = NaN;
    file.T_3m_c(file.T_3m_c > 1e6) = NaN;
    Ta = [Ta; mean(file.T_10m_c,'omitnan')];
    dT = [dT; mean(file.T_10m_c,'omitnan')-mean(file.T_3m_c,'omitnan')];
end    
[a b] = sort(Ta);


dtprof = [];
for j = 1:length(days)
    j
    
    theday = num2str(days(j));
    if length(theday) == 1; theday = ['0',theday]; end
    if j <= 17
        file = rd_netcdf(['isfs_sos_qc_geo_tiltcor_hr_202303',theday,'_23.nc'],'tc_2m_c','tc_3m_c','tc_5m_c','tc_10m_c', ...
        'tc_15m_c','tc_20m_c','w_2m_c','w_3m_c','w_5m_c','w_10m_c','w_15m_c','w_20m_c','T_10m_c','T_3m_c');
    else
        file = rd_netcdf(['isfs_sos_qc_geo_tiltcor_hr_202304',theday,'_23.nc'],'tc_2m_c','tc_3m_c','tc_5m_c','tc_10m_c', ...
        'tc_15m_c','tc_20m_c','w_2m_c','w_3m_c','w_5m_c','w_10m_c','w_15m_c','w_20m_c','T_10m_c','T_3m_c');
    end

    data1.dn20 = floor(data.dn(1)):1/1728000:ceil(data.dn(end)); data.dn20 = data.dn20(1:end-1)';
    
    file.T_10m_c(file.T_10m_c > 1e6) = NaN;
    file.T_3m_c(file.T_3m_c > 1e6) = NaN;

    file.tc_2m_c(file.tc_2m_c > 1e6) = NaN;
    file.tc_3m_c(file.tc_3m_c > 1e6) = NaN;
    file.tc_5m_c(file.tc_5m_c > 1e6) = NaN;
    file.tc_10m_c(file.tc_10m_c > 1e6) = NaN;
    file.tc_15m_c(file.tc_15m_c > 1e6) = NaN;
    file.tc_20m_c(file.tc_20m_c > 1e6) = NaN;
    
    file.w_2m_c(file.w_2m_c > 1e6) = NaN;
    file.w_3m_c(file.w_3m_c > 1e6) = NaN;
    file.w_5m_c(file.w_5m_c > 1e6) = NaN;
    file.w_10m_c(file.w_10m_c > 1e6) = NaN;
    file.w_15m_c(file.w_15m_c > 1e6) = NaN;
    file.w_20m_c(file.w_20m_c > 1e6) = NaN;
    
    data1.tc2 = file.tc_2m_c; data1.tc2 = reshape(data1.tc2,20*3600,1);
    data1.tc3 = file.tc_3m_c; data1.tc3 = reshape(data1.tc3,20*3600,1);
    data1.tc5 = file.tc_5m_c; data1.tc5 = reshape(data1.tc5,20*3600,1);
    data1.tc10 = file.tc_10m_c; data1.tc10 = reshape(data1.tc10,20*3600,1);
    data1.tc15 = file.tc_15m_c; data1.tc15 = reshape(data1.tc15,20*3600,1);
    data1.tc20 = file.tc_20m_c; data1.tc20 = reshape(data1.tc20,20*3600,1);
    
    data1.w2 = file.w_2m_c; data1.w2 = reshape(data1.w2,20*3600,1);
    data1.w3 = file.w_3m_c; data1.w3 = reshape(data1.w3,20*3600,1);
    data1.w5 = file.w_5m_c; data1.w5 = reshape(data1.w5,20*3600,1);
    data1.w10 = file.w_10m_c; data1.w10 = reshape(data1.w10,20*3600,1);
    data1.w15 = file.w_15m_c; data1.w15 = reshape(data1.w15,20*3600,1);
    data1.w20 = file.w_20m_c; data1.w20 = reshape(data1.w20,20*3600,1);
    
    
    % cmap = flipud(cbrewer('seq', 'Greens',20, 'linear'));
    % clf;
    % histN(data1.tc3-median(data1.tc3,'omitnan'),-2:0.01:2,'line','color',cmap(3,:)); hold on
    % histN(data1.tc5-median(data1.tc5,'omitnan'),-2:0.01:2,'line','color',cmap(5,:)); 
    % histN(data1.tc10-median(data1.tc10,'omitnan'),-2:0.01:2,'line','color',cmap(10,:));
    % histN(data1.tc15-median(data1.tc15,'omitnan'),-2:0.01:2,'line','color',cmap(15,:));
    % histN(data1.tc20-median(data1.tc20,'omitnan'),-2:0.01:2,'line','color',cmap(20,:));
    
    
    hts = [2, 3, 5, 10, 15, 20]-0.84;
    thevarstc = [var(data1.tc2-median(data1.tc2,'omitnan'),'omitnan') ...
     var(data1.tc3-median(data1.tc3,'omitnan'),'omitnan') ...
     var(data1.tc5-median(data1.tc5,'omitnan'),'omitnan') ...   
     var(data1.tc10-median(data1.tc10,'omitnan'),'omitnan') ...
     var(data1.tc15-median(data1.tc15,'omitnan'),'omitnan') ...
     var(data1.tc20-median(data1.tc20,'omitnan'),'omitnan')];
    thevarsw = [var(data1.w2-median(data1.w2,'omitnan'),'omitnan') ...
     var(data1.w3-median(data1.w3,'omitnan'),'omitnan') ...
     var(data1.w5-median(data1.w5,'omitnan'),'omitnan') ...   
     var(data1.w10-median(data1.w10,'omitnan'),'omitnan') ...
     var(data1.w15-median(data1.w15,'omitnan'),'omitnan') ...
     var(data1.w20-median(data1.w20,'omitnan'),'omitnan')];
    
    
    % clf;
    % cmap1 = colormap(lines);
    % t = tiledlayout(1,1);
    % ax1 = axes(t);
    % plot(thevarstc,hts,'color',cmap1(1,:),'linewidth',2);
    % ax1.Color = 'none';
    % ax1.XColor = cmap1(1,:);
    % xlim([0 0.4])
    % xlabel('\sigma^2 T [C]');
    % set(ax1,'fontsize',14);
    % axis square
    % ax2 = axes(t);
    % plot(thevarsw,hts,'color',cmap1(2,:),'linewidth',2);
    % ax2.XAxisLocation = 'top';
    % ax2.Color = 'none';
    % ax2.XColor = cmap1(2,:);
    % xlim([0 4])
    % xlabel('\sigma^2 w [C]'); ylabel('Height [m]');
    % axis square
    % set(ax2,'fontsize',14);
    
    prof = thevarstc./max(thevarstc);
    if abs(mean(file.T_10m_c,'omitnan')-mean(file.T_3m_c,'omitnan')) < 0.12 %mean(file.T_10m_c,'omitnan') > -0.1
        plot(prof,hts,'k','color',cmapr(40,:),'linewidth',2); hold on;
        plot(thevarsw./max(thevarsw),hts,'k','color',cmapb(40,:),'linewidth',2);
    else
        plot(prof,hts,'k','color',cmapr(10,:),'linewidth',2); hold on;
        plot(thevarsw./max(thevarsw),hts,'k','color',cmapb(10,:),'linewidth',2);
    end
    
    dtprof = [dtprof; prof(1)-prof(end)];

end


