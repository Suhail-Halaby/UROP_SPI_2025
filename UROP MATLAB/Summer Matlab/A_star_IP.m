clear
clc
close all

initialize_plane(1);

Konst = struct;

Konst.Ue = 15; 
Konst.alpha = 0;


Konst.time = [2025,8,25,0,0,0];


Konst.A_long = A_long;
Konst.B_long = B_long;
Konst.A_lat = A_lat;
Konst.B_lat = B_lat;

Konst.lat0 = init_lat;
Konst.long0 = init_long;
Konst.al0 = init_alt;

Konst.MG = MG;
Konst.MGs = {};
Konst.MGs{1} = MGs1;
Konst.MGs{2} = MGs2;
Konst.MGs{3} = MGs3;
Konst.MGs{4} = MGs4;
Konst.MGs{5} = MGs5;
Konst.MGs{6} = MGs6;


Konst.Elv = unique(energy_factor_lookup(:,1));
Konst.Az = unique(energy_factor_lookup(:,2));
Konst.Z = Z;

clearvars -except Konst


%% Calling A-star

% Positional Range
Lat_range  = 51.9:0.01:52.1;
Long_range = -0.1:0.01:0.1; 
[Lat_samp, Long_samp] = meshgrid( Lat_range , Long_range); 
Obs = zeros(size(Lat_samp));


% Obstacle Handilng -- feed forbidden nodes -- work in progress
Obs(10:15,10:15) = 1;   



% Domain assembly
domain_grid = cat(3,Lat_samp,Long_samp,Obs);

% Start and End Poses
start_pose = [51.92;-0.08];
end_pose =   [52.08;0.08];


[Ldist, Mdist, Ndist] = geodetic2enu(start_pose(1),start_pose(2),80, end_pose(1), end_pose(2),0,wgs84Ellipsoid);
path_len = sqrt(Ldist^2 + Mdist^2 + Ndist^2);

disp('Total (Flat) Dist (km) :')
disp(path_len/1000)


alt = 80;
evaluated = @(current,next) quasi_static_costs([current;alt],[next;alt],Konst,[current;alt],path_len);
heuristic = @(current,goal) quasi_static_costs([current;alt],[goal;alt],Konst,[current;alt],path_len);


plot_limits = [min(Lat_range),max(Lat_range),min(Long_range),max(Long_range)];

path_var = Astar_SP(domain_grid,start_pose,end_pose,evaluated,heuristic, plot_limits);

path_lat = domain_grid(1,path_var(1,:),1);
path_long = domain_grid(path_var(2,:),1,2);

plot(path_lat,path_long,'.-k',MarkerSize=14,MarkerEdgeColor='g',LineWidth=1.5)


%% Obstacle Plotting

for i = 10:15
    plot(Lat_range(i),Long_range(10:15),'*k',MarkerSize=15)
end

%% Plot Wind Feild
hold on
lat_plt  = unique(Konst.MG(:,1));
long_plt = unique(Konst.MG(:,2));

[LON, LAT] = meshgrid(long_plt, lat_plt);

MWS = Konst.MGs{3}(:,:,1);   % wind speed at t=1
MWD = Konst.MGs{4}(:,:,1);   % wind dir   at t=1

% Convert meteo direction (from which wind blows) → math convention
theta = deg2rad(MWD);

% Convert to u,v components
[U,V] = pol2cart(theta, MWS);

% Plot
plot(51.92,0.03,'.r',MarkerSize=10)

quiver(LAT, LON, U, V, 0.25, 'k');
ylabel('Longitude'); xlabel('Latitude');
legend('Closed Points','','Open Points','Optimal Path')
axis equal tight


%% Outputting the Waypoint Mission:

[path_y,path_x,path_z] = geodetic2enu(path_lat',path_long,alt,start_pose(1),start_pose(2),0,wgs84Ellipsoid);

path_z = alt*ones(size(path_z));


clearvars -except Konst path_x path_y path_z

wpts = zeros(6,length(path_x));

wpts(1:3,:) = [path_x';path_y';path_z'];

first_wp = wpts(1:3,1);

wpts(1:3,:) = wpts(1:3,:) - first_wp;

for i = 2:length(path_x)
    wpts(4:6,i) = wpts(1:3,i) - wpts(1:3,i-1);
end

wpts = wpts + zeros(6,1);

writematrix(wpts,'A_star_wpts.csv');

% Method

% specify control points
% feed in a set of gridded GPS coords, with altitudes
% preform A-star between 2 control points using the time (hour) of the first control point
% once a path is determined with the help of the quasi-steady costs, a
% waypoint mission will be created
% this is saved, and then the mission is flown with the simulink
% the time at the end control point is updated