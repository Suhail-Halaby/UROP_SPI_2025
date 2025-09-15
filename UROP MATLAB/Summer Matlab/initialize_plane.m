function [] = initialize_plane(qs_mode)
%% Battery
max_cap = 3600*236.80;
max_power_income = 1100*0.18*0.46208;
motor_power = 300;

%% Launch and Flight Parameters

% launch time
L_year = 2025;
L_month = 8;
L_day = 27;
L_hr = 8;

% launch location
init_alt = 60;
init_lat = 51;
init_long = 0;

% running patameters
lap_tog = 0;
laps = 10;
time_lim = 2*3600;

wind_tog = 1;
wind_type = 1;
const_wind_vector = [6;0;0];

%% Weather Function

MG = readmatrix('METEO_data.csv');


U_lat  = unique(MG(:,1));
U_lon  = unique(MG(:,2));
U_time = unique(MG(:,3));

nlat  = numel(U_lat);
nlon  = numel(U_lon);
ntime = numel(U_time);

MGs1 = nan(nlat, nlon, ntime);
MGs2 = nan(nlat, nlon, ntime);
MGs3 = nan(nlat, nlon, ntime);
MGs4 = nan(nlat, nlon, ntime);
MGs5 = nan(nlat, nlon, ntime);
MGs6 = nan(nlat, nlon, ntime);

[~, i] = ismember(MG(:,1), U_lat);
[~, j] = ismember(MG(:,2), U_lon);
[~, t] = ismember(MG(:,3), U_time);

% Convert subscripts to linear index
idx = sub2ind([nlat, nlon, ntime], i, j, t);

MGs1(idx) = MG(:,4);
MGs2(idx) = 270-MG(:,5);
MGs3(idx) = MG(:,6);
MGs4(idx) = 270-MG(:,7);
MGs5(idx) = MG(:,8);
MGs6(idx) = 270-MG(:,9);



%% System Matrices

A_long = [
            -0.0443238,            0.535343,                   0,               -9.81;
              -1.39239,            -8.47466,             11.4293,                   0;
           -0.00762359,            -4.30179,            -6.14639,                   0;
                     0,                   0,                   1,                   0 ] ;

A_lat = [
             -0.411379,           -0.649153,            -13.5366,                9.81,   0;
              -1.46667,            -16.7448,             3.78401,                   0,   0;
              0.292253,            -2.34281,           -0.384354,                   0,   0;
                     0,                   1,                   0,                   0,   0;
                     0,                   0,                   1,                   0,   0;] ;

% A = [A_long,zeros(4,5);zeros(5,4),A_lat];
%% Control Matricies

B_elv = [ 0.07638803; -11.96186; -43.11175; 0];

B_thr = [ 30/7 ; 0 ; 0 ; 0];

B_long  = [  B_elv , B_thr ];

B_ail = [ 1.178498; 35.6272; 2.449692; 0; 0];

B_rud = [ 3.77032; 0.8297701; -5.705696; 0; 0];

B_lat = [ B_ail , B_rud];

% B = [B_long,zeros(4,2);zeros(5,2),B_lat];
% 
% C = eye(9);

clear B_rud B_ail B_elv B_thr

%% Static Longitudinal TRIM

U_inf = 15;
Theta_eq = 0;

X_eq = [U_inf; 0; 0; Theta_eq];

trim_fun = @(xu) A_long*xu(1:4) + B_long*xu(5:6);
xu0 = [U_inf; 0; 0;Theta_eq;0;0];  % initial guess
xu_trim = fsolve(trim_fun, xu0);
x_trim = xu_trim(1:4);
u_trim = xu_trim(5:6);



if qs_mode == 0
    SimuPilot_Initializer();
end

%wpts = readmatrix('wpt_output.csv');
%wpts = wpts(:,2:end);

wp_rad = 50;
L1 = 100; 


%% SHADOW CASTER

energy_factor_lookup = readmatrix('Energy_Factor_Lookup.csv');

SZh = size(unique(energy_factor_lookup(:,2)),1);
SZv = size(unique(energy_factor_lookup(:,1)),1);

Z = reshape(energy_factor_lookup(:,3),SZh,SZv);

%% Back to workspace
vars = who;
for k = 1:numel(vars)
    if evalin('base', sprintf('exist(''%s'',''var'')', vars{k}))
        warning('Base workspace already has variable %s; it will be overwritten.', vars{k});
    end
    assignin('base', vars{k}, eval(vars{k}));
end



