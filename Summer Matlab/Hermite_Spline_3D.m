clear
clc
close all

wpts_x = [0,80,260,360,600];
wpts_y = [0,80,320,280,120];
wpts_h = [0,10,80,100,50];

hold on
plot(wpts_x,wpts_y,'.k',MarkerSize=20)

tng_x = [20,10,5,-10,-10];
tng_y = [0,20,-2,-200,1];
tng_h = [0,0,0,0,0];



% Hermite Spline
spl = @(s,Pi,Pj,Mi,Mj) (2*s^3-3*s^2+1)*Pi + (s^3-2*s^2+s)*Mi +(-2*s^3+3*s^2)*Pj +(s^3-s^2)*Mj;


A = zeros(2,100);
i = 1;

nodes = [wpts_x;wpts_y;wpts_h];
tangents = [tng_x;tng_y;tng_h];
seconds = zeros(size(nodes));
num_nodes = size(nodes,2);

mode = 1;
points = 100;
wind = [0;0;1];
[R,J0,AVC0,TNB,L_0] = Hermite_Path5(nodes,tangents,seconds,mode,points,wind);

plot(R(1,:),R(2,:),'.-r',MarkerEdgeColor='k',MarkerSize=2)


%% Path Optimisation

x0 = reshape(tangents, [], 1);  % 3N x 1 vector, necessary for fmincon
disp(x0)
x0b = reshape([tangents,seconds], [], 1);
disp(x0b)
% Optional bounds or constraints
lb = []; ub = [];  % or set limits on max tangent magnitudes if needed
minNormx = 1;
maxNormx = 1000;
minNormy = 0;
maxNormy = 10;


% Optimization
options = optimoptions('fmincon', ...
    'Display', 'iter', ...
    'Algorithm', 'sqp');
options.MaxFunctionEvaluations = 1e4;

[x_opt, fval] = fmincon(@(x) curvatureObjective(nodes,x(1:3*num_nodes), x(1+3*num_nodes:end), mode, points, wind), ...
                        x0b, [], [], [], [], [], [], @(x) unitNormConstraints(x(1:3*num_nodes),x(1+3*num_nodes:end),minNormx,maxNormx,minNormy,maxNormy) , options);
% Reshape result
tangents_opt = reshape(x_opt(1:3*num_nodes), 3, []);
seconds_opt = reshape(x_opt(1+3*num_nodes:end), 3, []);


[R_opt, J_opt, AVC_opt,TNB,L_opt] = Hermite_Path5(nodes, tangents_opt, seconds, mode, points, wind);
plot(R_opt(1,:),R_opt(2,:),'.-b',MarkerEdgeColor='k',MarkerSize=2)

disp('Cost Function Initial, Final Values')
disp(J0)
disp(J_opt)

disp('Average Radius of Curvature Initial, Final Values (m)')
disp(1/AVC0)
disp(1/AVC_opt)

disp('Average Path Length Initial, Final Values (m)')
disp(L_0)
disp(L_opt)

figure
hold on
axis equal padded
plot3(nodes(1,:),nodes(2,:),nodes(3,:),'.k',MarkerSize=20)
plot3(R(1,:),R(2,:),R(3,:),'r',LineWidth=2)
plot3(R_opt(1,:),R_opt(2,:),R_opt(3,:),'b',LineWidth=2)
legend('Nodes','Initial Path','Optimised Tangents')

% hold on 
% 
% f = 8;
% 
% figure
% axis equal padded
% for i = 1:size(R,2)
%     pause(0.05)
%     hold on
%     plot3(R_opt(1,:),R_opt(2,:),R_opt(3,:),'k',LineWidth=2)
%     quiver3(R_opt(1,i),R_opt(2,i),R_opt(3,i),f*TNB(1,i),f*TNB(2,i),f*TNB(3,i),'r',LineWidth=2)
%     quiver3(R_opt(1,i),R_opt(2,i),R_opt(3,i),f*TNB(4,i),f*TNB(5,i),f*TNB(6,i),'g',LineWidth=2)
%     quiver3(R_opt(1,i),R_opt(2,i),R_opt(3,i),f*TNB(7,i),f*TNB(8,i),f*TNB(9,i),'b',LineWidth=2)
%     quiver3(R_opt(1,i),R_opt(2,i),R_opt(3,i),-f*TNB(1,i),-f*TNB(2,i),-f*TNB(3,i),'r',LineWidth=2)
%     quiver3(R_opt(1,i),R_opt(2,i),R_opt(3,i),-f*TNB(4,i),-f*TNB(5,i),-f*TNB(6,i),'g',LineWidth=2)
%     quiver3(R_opt(1,i),R_opt(2,i),R_opt(3,i),-f*TNB(7,i),-f*TNB(8,i),-f*TNB(9,i),'b',LineWidth=2)
%     hold off
% end

%% Waypoint mission generation

% cumulative distance
N = size(R_opt,2);
dist = zeros(1,N);
for l = 2:N
    dist(l) = dist(l-1) + vecnorm(R_opt(1:3,l) - R_opt(1:3,l-1));
end

% desired number of waypoints
num_wpts = 10;
desired_d = linspace(0, dist(end), num_wpts);

% remove duplicate distances to avoid NaNs in interp1
[dist_u, ia] = unique(dist, 'stable');      % unique keeps order
idx_u = 1:N;
idx_u = idx_u(ia);                          % corresponding original indices

% interpolate to get fractional indices
frac_idx = interp1(dist_u, idx_u, desired_d, 'linear');

% round and clamp to valid integer indices
wp_ids = round(frac_idx);
wp_ids = max(1, min(N, wp_ids));  % ensure in [1, N]

% extract waypoints
R_aug = R(1:3, wp_ids);

% construct tracks
trk = zeros(size(R_aug));
for k = 2:size(R_aug,2)
    trk(:,k) = R_aug(:,k)-R_aug(:,k-1);
end

% if mode == 1
%     trk(:,end) = R_aug(:,2)-R_aug(:,end);
% end
mission = [R_aug;trk];
writematrix(mission, 'wpt_output.csv');



%% With wind
% wind = [5;0;0];
% 
% [x_opt, fval] = fmincon(@(x) curvatureObjective(x, nodes, mode, points, wind), ...
%                         x0, [], [], [], [], [], [], @(x) unitNormConstraints(x,minNorm,maxNorm) , options);
% % Reshape result
% tangents_opt = reshape(x_opt, 3, []);
% 
% 
% [R_optw, J_optw, AVC_opt] = Hermite_Path(nodes, tangents_opt, mode, points, wind);
% 
% plot3(R_optw(1,:),R_optw(2,:),R_optw(3,:),'g',LineWidth=2)
% 
% 
% 
% 
% 
function J = curvatureObjective(nodes,x,y,mode,points,wind)
    tangents = reshape(x, 3, []); 
    seconds =  reshape(y, 3, []); 
    [R, J] = Hermite_Path5(nodes,tangents,seconds, mode, points,wind);
end


function [c, ceq] = unitNormConstraints(x,y,minNormx,maxNormx,minNormy,maxNormy)
    M = reshape(x, 3, []);
    normsx = vecnorm(M);
    Q = reshape(y, 3, []);
    normsy = vecnorm(Q);
    ceq = [];  
    c = [minNormx - normsx, normsx - maxNormx,minNormy - normsy, normsy - maxNormy];  
end