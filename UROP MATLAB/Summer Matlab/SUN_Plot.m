clear
clc
close all

lat = deg2rad(51.5);
long = deg2rad(0);
R_e = 6.371e6;
S_vec = [cosd(232),sind(232),0];
S_vec = S_vec / norm(S_vec);

% position point on earth's surface
P = R_e*[cos(lat)*cos(long),cos(lat)*sin(long),sin(lat)];
P_quat = quaternion([0,P]);
% local frame

north = [-sin(lat)*cos(long),-sin(lat)*sin(long),cos(lat)]; 
east  = [-sin(long),cos(long),0];
up    = [cos(lat)*cos(long),cos(lat)*sin(long),sin(lat)];

frame = [north;east;up];

% position in tilted axis FOR
tilt = deg2rad(-23.5);
a = [sin(tilt),0,cos(tilt)];
north_pole = a';

q1 = quaternion([cos(tilt/2),0,sin(tilt/2),0]);

n_points = 25;
theta = linspace(0,2*pi,n_points);
P_rot = zeros(3,n_points);
F_rot = zeros(3,3,n_points);
SI = zeros(3,n_points);

% Canonical reference vectors
P0 = [0;0;1]; % North pole
U0 = [0;0;1]; % also used for Up calculation (radial from origin)


for i = 1:n_points
    q2 = quaternion([cos(theta(i)/2), sin(theta(i)/2)*a]);
    q_total = q2 * q1;  % apply q1 then q2

    % Position in ECEF
    r = rotatepoint(q_total, P)';     % 3x1 position vector
    P_rot(:,i) = r;

    % Up = normalized radial
    up = r / norm(r);

    % North = projection of rotated pole onto tangent plane
    north = north_pole - (dot(north_pole, up) * up);
    north = north / norm(north);

    % East = cross(north, up)
    east = cross(north, up);

    % Store NEU as rows
    F_rot(:,:,i) = [north'; east'; up'];

    % Sun dot prod
    SI(:,i) = [dot(S_vec,north) / norm(north) ;dot(S_vec,east) / norm(east);dot(S_vec,up) / norm(up)];
end



%% PLOTTING


[X, Y, Z] = sphere(50);
ap = [-a;[0,0,0];a];


hold on
axis equal

plot3(ap(:,1)*2*R_e,ap(:,2)*2*R_e,ap(:,3)*2*R_e,'k',LineWidth=2)
plot3(P_rot(1,:),P_rot(2,:),P_rot(3,:),'.r',MarkerSize=20)
N = reshape(F_rot(1,:,:),3,n_points)';
E = reshape(F_rot(2,:,:),3,n_points)';
U = reshape(F_rot(3,:,:),3,n_points)';
quiver3(P_rot(1,:),P_rot(2,:),P_rot(3,:),N(:,1)',N(:,2)',N(:,3)')
quiver3(P_rot(1,:),P_rot(2,:),P_rot(3,:),E(:,1)',E(:,2)',E(:,3)')
quiver3(P_rot(1,:),P_rot(2,:),P_rot(3,:),U(:,1)',U(:,2)',U(:,3)')
quiver3(0,0,0,S_vec(1),S_vec(2),S_vec(3),2*R_e,'b',LineWidth=2)
surf(X*R_e,Y*R_e,Z*R_e)
legend('Earth Axis','Points of Interest','North','East','Up','Direction of Sun')

xlabel('Earth-Centered X-Location (m)')
ylabel('Earth-Centered Y-Locationc (m)')
zlabel('Earth-Centered Z-Location (m)')

[ang,elv,rad] = cart2sph(SI(1,:),SI(2,:),SI(3,:));
ang = rad2deg(ang);
elv = rad2deg(elv);


figure

hold on
time = linspace(0,24,n_points);

plot(time,ang)
plot(time,elv)
title('Spherical Coord. Sun Progression')
legend('Azimith','Elvation')


figure

hold on
axis equal
O = zeros(1,size(SI,2));
quiver3(O,O,O,SI(1,:),SI(2,:),SI(3,:))
plot3(SI(1,:),SI(2,:),SI(3,:),'--.k',MarkerSize=12,MarkerEdgeColor='r')
 
[x,y] = meshgrid(-1:0.5:1); % Generate x and y data
z = zeros(size(x, 1)); % Generate z data
surf(x, y, z , FaceColor= 'none') % Plot the surface

plot3(cosd(ang).*cosd(elv),sind(ang).*cosd(elv),sind(elv),'.r',MarkerSize=12)
legend('Sun Vector','Visualized Pos.','Local Horizon')

xlabel('Non-dimensional North Direction')
ylabel('Non-dimensional East Direction')
zlabel('Non-dimensional Up Direction')


figure
sma = 149.598e9;
e = 0.0167;
phi = linspace(0,2*pi,100);
r_orb = @(phi,e)  sma * (1 - e^2 ) ./ (1 + e*cos(phi));
[X2,Y2] = pol2cart(phi,r_orb(phi,e));
hold on
plot(X2,Y2,'k--.',MarkerEdgeColor='r')
plot(sma*cos(phi),sma*sin(phi))
plot(0,0,'.k',MarkerSize=20)
title('Earth Orbit')


% %% Testing
% figure 
% hold on
% 
% t = [2025;3;21;0;0;0];
% hrs = 0:23;
% elv = zeros(length(hrs),1);
% ail = zeros(length(hrs),1);
% for i = 1:length(hrs)
%     time_info = datetime(t(1),t(2),t(3),hrs(i),0,0);
%     SI = solar_vector(time_info,51.2,0,0);
%    t2 = t;
%    t2(4) = hrs(i);
%     SI2 = solar_vector2(t2,51.2,0,0);
%     [ail(i),elv(i)] = cart2sph(SI(1),SI(2),SI(3));
%     [ail2(i),elv2(i)] = cart2sph(SI2(1),SI2(2),SI2(3));
% end
% ref_mid = datetime(t(1),t(2),t(3));
% [SR,SS,NOON] = sunrise(51.5,0,0,0,time_info);
% SR = datetime(SR, 'ConvertFrom', 'datenum');
% SS = datetime(SS, 'ConvertFrom', 'datenum');
% 
% 
% 
% plot(hrs,elv)
% plot(hrs,elv2,'b')
% xline(hours(SR-ref_mid),'r',LineWidth=2)
% xline(hours(SS-ref_mid),'r',LineWidth=2)
% xline(12)
% xline(6)
% xline(18)
% yline(0)
% yline(mean(elv),'r')
% 





% %% Testing
% figure
% hold on
% 
% t = datetime(2025,4,21,0,0,0);
% g = linspace(0,2*pi,100);
% gt = linspace(0,24,100);
% for i = 1:100
%     G = solar_vector(t,51.5,0,0,g(i));
%     [az(i),el(i)] = cart2sph(G(1),G(2),G(3));
% end
% 
% [SR,SS,NOON] = sunrise(51.5,0,0,0,t);
% SR = datetime(SR, 'ConvertFrom', 'datenum');
% SS = datetime(SS, 'ConvertFrom', 'datenum');
% 
% ref_mid = t;
% 
% plot(gt,el)
% xline(hours(SR-ref_mid),'r')
% xline(hours(SS-ref_mid),'r')
% xline(12)
% xline(6)
% xline(18)
% yline(0)
% yline(mean(el),'r')