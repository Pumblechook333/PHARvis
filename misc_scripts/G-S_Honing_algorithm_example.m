
clear

raytrace_params.origin_lat = -77.8464;             % latitude of the start point of rays
raytrace_params.origin_lon = 166.6683;            % longitude of the start point of rays
raytrace_params.origin_alt = 4; % altitude of the start point of rays
raytrace_params.tol=[1e-12 0.005 1];

raytrace_params.nhops=2;% number of hops attempted
raytrace_params.init_delta_elev=1;% elevation angle perturb
raytrace_params.init_delta_bear=10;% bearing angle perturb
raytrace_params.fr=10.3; % frequency in MHz
Newton_params.elev_tol=18; % lowest elevation angle attempted
Newton_params.iteration_num=10; % number of Newton's iteration
Newton_params.rece_dis=100; % size of receiving zone (m)

ht_start = 90;          % start height for ionospheric grid (km)
ht_inc = 10;             % height increment step length (km)
num_ht = 192;
lat_start = -90;
lat_inc = 1;
num_lat = 46;
lon_start= 0;
lon_inc = 4;
num_lon = 90;
iono_grid_parms = [lat_start, lat_inc, num_lat, lon_start, lon_inc, num_lon, ...
    ht_start, ht_inc, num_ht ];
B_ht_start = ht_start;          % start height for geomagnetic grid (km)
B_ht_inc =10;                  % height increment (km)
B_num_ht = ceil(num_ht .* ht_inc ./ B_ht_inc);
B_lat_start = lat_start;
B_lat_inc = 1.0;
B_num_lat = ceil(num_lat .* lat_inc ./ B_lat_inc);
B_lon_start = lon_start;
B_lon_inc = 4.0;
B_num_lon = ceil(num_lon .* lon_inc ./ B_lon_inc);
geomag_grid_parms = [B_lat_start, B_lat_inc, B_num_lat, B_lon_start, ...
    B_lon_inc, B_num_lon, B_ht_start, B_ht_inc, B_num_ht];
load('D:\eclipse_ionosphere\NoEclips_inosphere_linear_final_with_B_6_54.mat')

lat=h5read('D:\Doppler_newton_method\Dec_4\RRI_20211204_065415_070412_lv1_13.0.0.h5','/CASSIOPE Ephemeris/Geographic Latitude (deg)');
lon=h5read('D:\Doppler_newton_method\Dec_4\RRI_20211204_065415_070412_lv1_13.0.0.h5','/CASSIOPE Ephemeris/Geographic Longitude (deg)');
alt=h5read('D:\Doppler_newton_method\Dec_4\RRI_20211204_065415_070412_lv1_13.0.0.h5','/CASSIOPE Ephemeris/Altitude (km)')*1000;% spacecraft height (m)
satellite_location.lon=lon;
satellite_location.lat=lat;
satellite_location.alt=alt;
% load('D:\Doppler_newton_method\Dec_8\Eclips_inosphere_final_with_B_Dec_8_6_30.mat')
% Bx=B_x;
% By=B_y;
% Bz=B_z;
Ne_5=zeros(size(Ne));
%  Ne(90-76:90-70,35:38,:)=Ne(90-76:90-70,35:38,:)*2;
%      Ne=Ne*factor;

raytrace_params.OX_mode=1;
[ ~,ray_O_receive]=Honing_raytracing(satellite_location,...
    Newton_params,raytrace_params,iono_grid_parms,geomag_grid_parms,Ne,Ne_5,collision_freq,...
    Bx,By,Bz);
save('D:\Doppler_newton_method\Dec_4\epop_ray_receive_noeclipse_OmodeDec4_6to7UT.mat','ray_O_receive');% reverse eclipse and non-eclipse filename
%%
for i=1:length(ray_O_receive)
    lon=ray_O_receive(i).lon;
    lon(lon<0)=lon(lon<0)+360;
    scatter3(ray_O_receive(i).lat,lon,ray_O_receive(i).height,5,ray_O_receive(i).electron_density);
    hold on
end
xlabel('Latitude');
ylabel('Longitude');
zlabel('Altitude km');
colormap(viridis)
c=colorbar;
c.Label.String = 'Ne cm^-^3';
title('SAMI3 Dec 8 6:21:14-6:29:11 O mode');
%%
i=1;

% load ionosphere and magnetic fiel information here
% grid size of ionosphere and magnetic field must match the parameters from
% line 16 to line 37 above.
%     ray_O_receive=struct('lat',[],'lon',[],'height',[],'hour',[],'minute',[],'elevs',[],'bearing',[],'ground_range',[],...
%         'group_range',[],'absorption',[],'ionofactor',[]);
%     ray_X_receive=struct('lat',[],'lon',[],'height',[],'hour',[],'minute',[],'elevs',[],'bearing',[],'ground_range',[],...
%         'group_range',[],'absorption',[],'ionofactor',[]);
for factor=1:0.1:10
    load('D:\eclipse_ionosphere\Eclips_inosphere_final_with_B_6_54.mat')
    Bx=B_x;
    By=B_y;
    Bz=B_z;
    Ne_5=zeros(size(Ne));
    %  Ne(90-76:90-70,35:38,:)=Ne(90-76:90-70,35:38,:)*2;
    Ne=Ne*factor;
    %  raytrace_params.OX_mode=1;%  O mode
    % [~,ray_O_receive{i}]=Honing_raytracing(satellite_location,...
    %     Newton_params,raytrace_params,iono_grid_parms,geomag_grid_parms,Ne,Ne_5,collision_fr,...
    %     Bx,By,Bz);
    raytrace_params.OX_mode=1;%  X mode
    [ ~,ray_O_receive{i}]=Honing_raytracing(satellite_location,...
        Newton_params,raytrace_params,iono_grid_parms,geomag_grid_parms,Ne,Ne_5,collision_fr,...
        Bx,By,Bz);
    %         ray_O_receive{i}.ionofactor=factor;


    %         raytrace_params.OX_mode=-1;%  X mode
    %         [ ~,ray_X_receive{i}]=Honing_raytracing(satellite_location,...
    %             Newton_params,raytrace_params,iono_grid_parms,geomag_grid_parms,Ne,Ne_5,collision_fr,...
    %             Bx,By,Bz);
    %         ray_X_receive{i}.ionofactor=factor;
    factor_index(i)=factor;
    i=i+1;
end

%    save('D:\Doppler_newton_method\Dec_7\ray_O_receive_epop_doppler_Dec_07_7_0_100m_threshold_V2','ray_O_receive');
%    save('D:\Doppler_newton_method\Dec_7\ray_X_receive_epop_doppler_Dec_07_7_0_100m_threshold_V2','ray_X_receive');

%%
% find the proper ionosphere for eclipse raytracing December 4 2021
% 065415_070412 UT e-POP conjunction
% ray_O_receive_cell=ray_O_receive;


figure(2)
count=1;
ionofacotr=0.1:0.1:10;
for ionofactorcount=1:length(ionofacotr)
    filname=['Doppler_eclipsediff_iono_8UT_Xmode_factor_',num2str(ionofactorcount),'.mat'];
    if exist(filname,'file')
        % for count=1:length(ray_O_receive_cell)
        load(filname);
        %     clear vars ray_O_receive phase_path timeindex
        %     ray_O_receive=ray_O_receive_cell{count};
        for i=1:length(ray_X_receive)
            phase_path(i)=ray_X_receive(i).phase_path(end);
            timeindex(i)=ray_X_receive(i).loca_num;
        end
        phasepath_diff=diff(phase_path);
        time_diff=diff(timeindex);
        timefit=timeindex(1:end-1);
        fre=10.3e6;
        c=2.99792458e8;
        lambda=c/fre;
        doppler=-1/lambda*(phasepath_diff*1000./time_diff);
        doppler(doppler>300)=nan;
        doppler(doppler<-300)=nan;
        [doppler_O_fit,S]=polyfit(timeindex(1:end-1)',doppler',8);
        [O_poly,delta]=polyval(doppler_O_fit,timeindex(1:end-1),S);
        Doppler_diff_iono(count).doppler=doppler;
        Doppler_diff_iono(count).timecount=timefit;
        Doppler_diff_iono(count).doppler_poly=O_poly;
        % timescale=hour*ones(size(O_poly));
        % cmap = colormap('jet');  % 使用 'jet' 颜色映射
        O=plot(timeindex(1:end-1),doppler,'o');
        hold on
        %     O.LineWidth=1.5;
        clear vars ray_X_receive
        factorindex_array(count)=factor_index;
        count=count+1;
    end
end
xlabel('Second');
ylabel('Hz');
title('Diff IRI raytracing Doppler');
grid on
% save('Dec_4_Xmode_8UT_different_iono_doppler_all.mat', 'Doppler_diff_iono','factorindex_array');
%%
% i=2;
clear
load('D:\eclipse_ionosphere\Eclips_inosphere_linear_final_with_B_6_54.mat')
TEC=sum(Ne,3)*10000;
figure(1)
contourf(0:4:356,Latitude_x,TEC*1e6*1e-16,20);
colormap(viridis)
h = colorbar;
set(get(h,'label'),'string','TECU','FontSize',12);%name the colorbar
hold on
load('D:\Doppler_newton_method\Dec4_6to7UT_Diff_iono\Dec4_Omode_eclipse_6to7UT_ray_path.mat');
ray_time_select=ray_index_time(ray_index_time<310);

lat=h5read('D:\Doppler_newton_method\Dec_4\RRI_20211204_065415_070412_lv1_13.0.0.h5','/CASSIOPE Ephemeris/Geographic Latitude (deg)');
lon=h5read('D:\Doppler_newton_method\Dec_4\RRI_20211204_065415_070412_lv1_13.0.0.h5','/CASSIOPE Ephemeris/Geographic Longitude (deg)');
alt=h5read('D:\Doppler_newton_method\Dec_4\RRI_20211204_065415_070412_lv1_13.0.0.h5','/CASSIOPE Ephemeris/Altitude (km)')*1000;% spacecraft height (m)
lon(lon<0)=lon(lon<0)+360;
scatter(lon,lat,20,'o','filled','y');
for i=1:length(ray_time_select)
    ray_lat=ray_O_receive_match(i).lat;
    ray_lon=ray_O_receive_match(i).lon;
    ray_alt=ray_O_receive_match(i).height;
    lat_select=ray_lat(ray_alt<300&ray_alt>200);
    lon_select=ray_lon(ray_alt<300&ray_alt>200);
    %  alt_select=ray_alt(ray_alt<300&ray_alt>200);
    scatter(ray_lon,ray_lat,20,'o','filled');
end
scatter(166.6683,-77.8464,100,'p','filled','r');
scatter(155.894,-79.147,100,'d','filled');
xlabel('Longitude');
ylabel('Latitude');
title('SAIM3 Eclipse Dec 4, 2021 6:54 UT');
set(gca,'FontSize',12);

%%
load('D:\Doppler_newton_method\Dec_4\diff_iono\Dec_4_Xmode_8UT_different_iono_doppler_all.mat');
time_matrix=nan(100,600);
doppler_matrix=nan(100,600);
for i=1:length(Doppler_diff_iono)
    time=Doppler_diff_iono(i).timecount;
    sim_doppler=Doppler_diff_iono(i).doppler;
    for time_count=1:length(time)
        time_matrix(i,time(time_count))=time(time_count);
        doppler_matrix(i,time(time_count))=sim_doppler(time_count);
    end
end


load('D:\Doppler_newton_method\Dec_4\Doppler_20211204_083615_084612_1.5second_interval_ionosphere.mat')
epop_time=t_second(~isnan(doppler_1));
epop_doppler=doppler_1(~isnan(doppler_1));
i=1;
j=1;

for T_count=1:length(epop_time)
    for time_count=1:600
        if abs(time_count-epop_time(T_count))<1
            doppler_diff=doppler_matrix(:,time_count)-epop_doppler(T_count);
            if any(~isnan(doppler_diff))
                [select_doppler,iono_index]=min(abs(doppler_diff));
                if abs(select_doppler)<6
                    ray_index_iono(i)=iono_index;
                    ray_index_time(j)=time_count;

                    j=j+1;
                    i=i+1;

                end
            end
        end
    end
end
for i=1:length(ray_index_iono)
    time_sim(i)=time_matrix(ray_index_iono(i),ray_index_time(i));
    doppler_sim(i)=doppler_matrix(ray_index_iono(i),ray_index_time(i));
end
plot(time_sim,doppler_sim,'o');
hold on
plot(t_second,doppler_1,'d');
grid on
xlabel('Time Second');
ylabel('Hz');
title('Doppler Dec 4,2021 6:54:15-7:4:12');
legend('Simulation O mode','RRI');
%  for i=1:length(ray_index_iono)
%  iono_factor_select(i)=
%  endD
figure(2)
T_sim_UT=time_sim;
plot(T_sim_UT,ray_index_iono*0.1+1);
grid on
xlabel('UT');
ylabel('Iono factor');
title('Doppler Dec 4,2021 6:54:15-7:4:12');
% save('Doppler_match_Dec_4_065415_070412.mat');
%%
factor=1:0.1:10;
load('Doppler_match_Dec_4_065415_070412.mat');
ray_O_receive_match=struct('lat',[],'lon',[],'height',[],'elevs',[],'bearing',[],'loca_num',[],'absorption',[],'group_range',[],...
    'electron_density',[],'phase_path',[],'refractive_index',[],'range',[]);
for i=1:length(ray_index_iono)
    load(['Doppler_diff_iono_factor',num2str(factor(ray_index_iono(i))),'Omode_6to7UT_eclipse.mat']);
    for timeindex=1:length(ray_O_receive)
        time(timeindex)=ray_O_receive(timeindex).loca_num;
    end
    raynumber=find(time==ray_index_time(i));
    ray_O_receive_match(i)=ray_O_receive(raynumber);
    clearvars time
end
%%
figure(3)
for i=1:190
    lon=ray_O_receive_match(i).lon;
    lon(lon<0)=lon(lon<0)+360;
    scatter3(ray_O_receive_match(i).lat,lon,ray_O_receive_match(i).height,5,ray_O_receive_match(i).electron_density);
    % drawnow
    % pause(0.2)
    hold on
end
xlabel('Latitude');
ylabel('Longitude');
zlabel('Altitude km');
colormap(viridis)
c=colorbar;
set(get(c,'label'),'string','Ne cm^-^3','FontSize',12);%name the colorbar
title('Dec 4 6:54:15-7:4:12 O mode');
%%


% December 8 6UT Omode


clear
% figure(1)
count=1;
for ionofactor=1:0.1:10

    if exist(['Doppler_diff_iono_Omode_Dec_8_6UT_factor',num2str(ionofactor),'.mat'],'file')
        % for count=1:length(ray_O_receive_cell)
        load(['Doppler_diff_iono_Omode_Dec_8_6UT_factor',num2str(ionofactor),'.mat']);
        %     clear vars ray_O_receive phase_path timeindex
        %     ray_O_receive=ray_O_receive_cell{count};
        for i=1:length(ray_O_receive)
            phase_path(i)=ray_O_receive(i).phase_path(end);
            timeindex(i)=ray_O_receive(i).loca_num;
        end
        phasepath_diff=diff(phase_path);
        time_diff=diff(timeindex);
        timefit=timeindex(1:end-1);
        fre=10.3e6;
        c=2.99792458e8;
        lambda=c/fre;
        doppler=-1/lambda*(phasepath_diff*1000./time_diff);
        doppler(doppler>300)=nan;
        doppler(doppler<-300)=nan;
        [doppler_O_fit,S]=polyfit(timeindex(1:end-1)',doppler',8);
        [O_poly,delta]=polyval(doppler_O_fit,timeindex(1:end-1),S);
        Doppler_diff_iono(count).doppler=doppler;
        Doppler_diff_iono(count).timecount=timefit;
        Doppler_diff_iono(count).doppler_poly=O_poly;
        % timescale=hour*ones(size(O_poly));
        % cmap = colormap('jet');  % 使用 'jet' 颜色映射
        O=plot(timeindex(1:end-1),doppler,'o');
        hold on
        %     O.LineWidth=1.5;
        clear vars ray_O_receive
        factorindex_array(count)=factor_index;
        count=count+1;
    end
end
%  save('D:\Doppler_newton_method\Dec_8\different_ionosphere\Dec_8_Omode_6UT_different_iono_doppler_all.mat', 'Doppler_diff_iono','factorindex_array');
%%
% i=2;
clear

figure(1)

load('D:\Doppler_newton_method\Dec_8\different_ionosphere\Dec_8_Omode_6UT_different_iono_doppler_all.mat');
time_matrix=nan(100,600);
doppler_matrix=nan(100,600);
for i=1:length(Doppler_diff_iono)
    time=Doppler_diff_iono(i).timecount;
    sim_doppler=Doppler_diff_iono(i).doppler;
    for time_count=1:length(time)
        time_matrix(i,time(time_count))=time(time_count);
        doppler_matrix(i,time(time_count))=sim_doppler(time_count);
    end
end


load('D:\Doppler_newton_method\Dec_8\Doppler_20211208_062114_062911_1.5second_interval_ionosphere.mat')
epop_time=t_second(~isnan(doppler_1));
epop_doppler=doppler_1(~isnan(doppler_1));
i=1;
j=1;

for T_count=1:length(epop_time)
    for time_count=1:600
        if abs(time_count-epop_time(T_count))<1
            doppler_diff=doppler_matrix(:,time_count)-epop_doppler(T_count);
            if any(~isnan(doppler_diff))
                [select_doppler,iono_index]=min(abs(doppler_diff));
                if abs(select_doppler)<3
                    ray_index_iono(i)=iono_index;
                    ray_index_time(j)=time_count;

                    j=j+1;
                    i=i+1;

                end
            end
        end
    end
end
for i=1:length(ray_index_iono)
    time_sim(i)=time_matrix(ray_index_iono(i),ray_index_time(i));
    doppler_sim(i)=doppler_matrix(ray_index_iono(i),ray_index_time(i));
end
plot(time_sim,doppler_sim,'o');
hold on
plot(t_second,doppler_1,'d');
grid on
xlabel('Time Second');
ylabel('Hz');
title('Doppler Dec 8,2021 6:21:14-6:29:11');
legend('Simulation O mode','RRI');
%  for i=1:length(ray_index_iono)
%  iono_factor_select(i)=
%  end
figure(2)
T_sim_UT=time_sim;
plot(T_sim_UT,ray_index_iono*0.1+1);
grid on
xlabel('UT');
ylabel('Iono factor');
title('Iono index Dec 8,2021 6:21:14-6:29:11');
save('D:\Doppler_newton_method\Dec_8\different_ionosphere\Doppler_match_Dec_8_062114_062911_3Hz.mat');
%%
factor=1:0.1:10;
load('D:\Doppler_newton_method\Dec_8\different_ionosphere\Doppler_match_Dec_8_062114_062911_3Hz.mat');
ray_O_receive_match=struct('lat',[],'lon',[],'height',[],'elevs',[],'bearing',[],'loca_num',[],'absorption',[],'group_range',[],...
    'electron_density',[],'phase_path',[],'refractive_index',[],'range',[]);
for i=1:length(ray_index_iono)
    load(['D:\Doppler_newton_method\Dec_8\different_ionosphere\Doppler_diff_iono_Omode_Dec_8_6UT_factor',num2str(factor(ray_index_iono(i))),'.mat']);
    for timeindex=1:length(ray_O_receive)
        time(timeindex)=ray_O_receive(timeindex).loca_num;
    end

    raynumber=find(time==ray_index_time(i));
    ray_O_receive_match(i)=ray_O_receive(raynumber);
    clearvars time
end
save('D:\Doppler_newton_method\Dec_8\different_ionosphere\Dec8_Omode_6UT_ray_path_3Hz.mat','ray_O_receive_match','ray_index_iono','ray_index_time');
%%
clear
load('D:\Doppler_newton_method\Dec_8\different_ionosphere\Dec8_Omode_6UT_ray_path_3Hz.mat','ray_O_receive_match','ray_index_iono','ray_index_time');
figure(11)

for i=1:length(ray_O_receive_match)
    lon=ray_O_receive_match(i).lon;
    lon(lon<0)=lon(lon<0)+360;
    scatter3(ray_O_receive_match(i).lat,lon,ray_O_receive_match(i).height,5,ray_O_receive_match(i).electron_density);
    % drawnow
    % pause(0.2)
    hold on
end
xlabel('Latitude');
ylabel('Longitude');
zlabel('Altitude km');
colormap(viridis)
c=colorbar;
c.Label.String = 'Ne cm^-^3';
% caxis([0,3e5]);
title('Dec 8 6:21:14-6:29:11 O mode');
set(get(c,'label'),'string','Ne cm^-^3','FontSize',12);%name the colorbar