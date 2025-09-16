clear 
clc
close all
hold on

[X,Y] = meshgrid(0:1:50,0:1:50);  
Obs = zeros(size(X));
Obs(1:30,30) = 1;


domain_grid = cat(3,X,Y,Obs);
start_pose = [1;1];
end_pose = [50;1];
evaluated = @(current,next) sqrt((current(1)-next(1))^2 + (current(2)-next(2))^2);
heuristic = @(current,goal) sqrt((current(1)-goal(1))^2 + (current(2)-goal(2))^2);



plot_limits = [0,50,0,50];

path_var = Astar_SP(domain_grid,start_pose,end_pose,evaluated,heuristic,plot_limits);

x = domain_grid(1,path_var(1,:),1);
y = domain_grid(path_var(2,:),1,2);

plot(x,y,'.-k',MarkerSize=8,MarkerEdgeColor='g')
axis( [0 50 0 50] )

Obs_x = [];
Obs_y = [];
for i  = 1:size(Obs,1)
    for j = 1:size(Obs,2)
        if Obs(i,j) == 1
            Obs_y = [Obs_y domain_grid(1,i,1)];
            Obs_x = [Obs_x domain_grid(j,1,2)];
        end
    end
end

plot(Obs_x,Obs_y,'k',LineWidth=2)