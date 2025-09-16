function [net_cost] = quasi_static_costs(coord1,coord2,Konst,cp_ref,path_len)
%% Initializing Aircraft Model Position

[Y,X,h] = geodetic2enu(coord1(1),coord1(2),coord1(3),cp_ref(1),cp_ref(2),cp_ref(3),wgs84Ellipsoid);
point1 = [X,Y,h];

[Y,X,h] = geodetic2enu(coord2(1),coord2(2),coord2(3),cp_ref(1),cp_ref(2),cp_ref(3),wgs84Ellipsoid);
point2 = [X,Y,h];

clear X Y h


%% Straight-Line Processing

% poits consist of [X,Y,Z] locations

% X --> North
% Y --> East
% Z --> Up

p2p_Ndist = sqrt( sum((point1(1:3) - point2(1:3)).^2) );
p2p_Hdist = sqrt( sum((point1(1:2) - point2(1:2)).^2) );

p2p_clmgrad =  atan(   (point2(3) - point1(3))   / p2p_Hdist );
p2p_heading = cart2pol(  point2(1) - point1(1) , point2(2) - point1(2) );

p2p_tng = (point2 - point1) ./ vecnorm(point2 - point1);


%% Finding Midpoint Position

center_pos = (point2+point1) ./ 2;

[cen_lat,cen_long,cen_alt] = enu2geodetic(center_pos(2),center_pos(1),center_pos(3),cp_ref(1),cp_ref(2),cp_ref(3),wgs84Ellipsoid);



%% Average Wind Conditions for Track Midpoint

MG = Konst.MG; 
time = Konst.time(4);

% Low level wind

LLW_speed_rec = Konst.MGs{1};
 
LLW_speed = interp3(unique(MG(:,2)),unique(MG(:,1)),unique(MG(:,3)),LLW_speed_rec,cen_long,cen_lat,time);

LLW_dir_rec = Konst.MGs{2};
 
LLW_dir = interp3(unique(MG(:,2)),unique(MG(:,1)),unique(MG(:,3)),LLW_dir_rec,cen_long,cen_lat,time);

[LLWx, LLWy] = pol2cart(LLW_dir,LLW_speed);

% Med level wind


MLW_speed_rec = Konst.MGs{3};
 
MLW_speed = interp3(unique(MG(:,2)),unique(MG(:,1)),unique(MG(:,3)),MLW_speed_rec,cen_long,cen_lat,time);

MLW_dir_rec = Konst.MGs{4};
 
MLW_dir = interp3(unique(MG(:,2)),unique(MG(:,1)),unique(MG(:,3)),MLW_dir_rec,cen_long,cen_lat,time);

[MLWx, MLWy] = pol2cart(MLW_dir,MLW_speed);


% High Wind Level

HLW_speed_rec = Konst.MGs{5};
 
HLW_speed = interp3(unique(MG(:,2)),unique(MG(:,1)),unique(MG(:,3)),HLW_speed_rec,cen_long,cen_lat,time);

HLW_dir_rec = Konst.MGs{6};
 
HLW_dir = interp3(unique(MG(:,2)),unique(MG(:,1)),unique(MG(:,3)),HLW_dir_rec,cen_long,cen_lat,time);

[HLWx, HLWy] = pol2cart(HLW_dir,HLW_speed);


% Altitude Interpolation
Vx = [LLWx,MLWx,HLWx];
Vy = [LLWy,MLWy,HLWy];
Alts = [10,80,120];

Wind_X = interp1(Alts,Vx,cen_alt);
Wind_Y = interp1(Alts,Vy,cen_alt);

%% Resolve Wind into body-ish frame

[b_theta,~] = cart2pol(p2p_tng(1),p2p_tng(2));

U_comp =  Wind_X*cos(b_theta) + Wind_Y*sin(b_theta);
V_comp = -Wind_X*sin(b_theta) + Wind_Y*cos(b_theta);


Wind_UW = [U_comp;0];
Wind_V = V_comp;

%% Longitudinal QS Inputs

U_e = Konst.Ue;
alpha = Konst.alpha;
gamma = p2p_clmgrad; 

X_long_qs = [U_e;U_e*tan(alpha);0;gamma-alpha];

U_long_qs = -pinv(Konst.B_long)*(Konst.A_long*X_long_qs + Konst.A_long(:,1:2)*Wind_UW ); 

%VER1 = Konst.A_long*X_long_qs + Konst.A_long(:,1:2)*Wind_UW + Konst.B_long*U_long_qs;

%% Lateral QS Inputs

heading = p2p_heading;

X_lat_qs = [0;0;0;0;heading];

U_lat_qs = -pinv(Konst.B_lat)*(Konst.A_lat*X_lat_qs + Konst.A_lat(:,1)*Wind_V ); 

%VER2 = Konst.A_lat*X_lat_qs + Konst.A_lat(:,1)*Wind_V + Konst.B_lat*U_lat_qs;

%% Solar Income Model

% solar vector
time_info = Konst.time;
SI = solar_vector2(time_info,cen_lat,cen_long,cen_alt);

% casting into plane FOR
Euler = [X_long_qs(3),X_lat_qs(4),X_lat_qs(5)];
[relative_elv, relative_azimuth] = relative_sun(SI,Euler);

% solar income factor
Z = Konst.Z;

    
% relative azimuth must be in the proper range
if  relative_azimuth < 0
    relative_azimuth = 2*pi + relative_azimuth; 
end

% elevation must be clipped 
clip_elv = clip(relative_elv,min(Konst.Elv),pi/2);

SIF = interp2(Konst.Elv,Konst.Az,Z,clip_elv,relative_azimuth);

%% Defining the General Cost Function

% Contribution of Thrust Power Draw
C_thr = cost_fun(0,1,2,U_lat_qs(2));

% Solar Accumulation Factor

C_SIF = cost_fun(0,1,2,1-SIF);

net_cost = p2p_Ndist*(C_thr+C_SIF);





end


%% Defining the archetype normalized cost function 

function C =  cost_fun(threshold,limit,exponent,value)
    C = heaviside( sign(limit-threshold)*(value-threshold) ) * ...
        ( exp( exponent*(value-threshold) / (limit-threshold) )  - 1 ) / ...
        ( exp(exponent) - 1 );

    C = min(max(C,0),1);
end