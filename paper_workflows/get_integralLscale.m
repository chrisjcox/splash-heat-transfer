
% Calculate for the whole experiment

cd /PSL/Observations/Campaigns/SPLASH/asfs50/2_level_product/zenodo_archive_v1/

d = dir('sledwind10hz*.nc');

intLens.DN = [];
intLens.L = [];

for k = 1:length(d)-1

    disp([num2str(round(k/length(d)*10000)/100),'% Complete...']);
    x1 = rd_netcdf(d(k).name);
    x2 = rd_netcdf(d(k+1).name);
    dn1 = double(x1.base_time)./86400 + x1.time_offset./86400000 + datenum([1970 1 1 0 0 0]);
    dn2 = double(x2.base_time)./86400 + x2.time_offset./86400000 + datenum([1970 1 1 0 0 0]);
    dn = [dn1;dn2];
    x1.metek_u(x1.metek_u < -100) = NaN;
    x1.metek_u(x1.metek_v < -100) = NaN;
    x1.metek_u(x1.metek_w < -100) = NaN;
    x2.metek_u(x2.metek_u < -100) = NaN;
    x2.metek_u(x2.metek_v < -100) = NaN;
    x2.metek_u(x2.metek_w < -100) = NaN;
    u1 = sqrt(x1.metek_u.^2+x1.metek_v.^2+x1.metek_w.^2); 
    u2 = sqrt(x2.metek_u.^2+x2.metek_v.^2+x2.metek_w.^2); 
    try
        u = naninterp([u1;u2]);
    catch
        continue
    end

    % Calculate
    Ltmp = integral_length_scale(u,dn);

    intLens.DN = [intLens.DN DN];
    intLens.L = [intLens.L Ltmp];

end

clearallbut intLens
cd /home/ccox/
save intLens