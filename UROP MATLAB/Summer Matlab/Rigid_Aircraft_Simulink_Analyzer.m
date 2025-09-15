%wpts(1:2,:) = wpts(1:2,:)*1e6;
modelName = 'Uncoupled_Aircraft_State_Space';
load_system(modelName);
simout = sim(modelName);

close all



%% Trajectory 

Time = simout.position.Time;
Pos = simout.position.Data;

points = 100;
theta_circ = linspace(0,2*pi,points);
x_sup = wp_rad*cos(theta_circ);
y_sup = wp_rad*sin(theta_circ);


hold on
plot3(Pos(:,1),Pos(:,2),Pos(:,4),'k',LineWidth=2)
xlabel('X (m)')
ylabel('Y (m)')
zlabel('H (m)')

plot3(wpts(1,:),wpts(2,:),wpts(3,:),'r.',MarkerSize=15)


for i = 1:size(wpts,2)
    plot3(wpts(1,i)+x_sup,wpts(2,i)+y_sup,wpts(3,i)*ones(size(x_sup)),'r')
end

hold off
axis equal
axis padded


%% State of Charge and Solar Income
figure 

Launch_hr = simout.global_time.Data(1,4);

Time = simout.SOC.Time;   % [s]
SOC  = simout.SOC.Data;

% Convert time into hours of day, wrap every 24 h
timeHours = Time/3600;

subplot(2,1,1)
plot(timeHours, 100*SOC/max_cap,'k','LineWidth',2)
title('State of Charge')
xlabel('Time of Day (hh:mm)')
ylabel('State of Charge (%)')

% Format x-axis ticks as hh:mm
xticks(0:0.5:24)                                % tick every 2 hours
xticklabels(datestr(xticks/24, 'HH:MM'))

SI = simout.power_income.Data ;

dt = 0.05;
tq = Time(1):dt:Time(end);
SIq = interp1(Time,SI,tq,'linear');

% Moving averages
k1_interval = 60/dt;     % 1 min
k2_interval = 1800/dt;    % 10 min
k3_interval = 3600/dt;   % 1 hr

MSI_k1 = movmean(SIq,k1_interval);
MSI_k2 = movmean(SIq,k2_interval);
MSI_k3 = movmean(SIq,k3_interval);

% Wrap query time to hours of day
timeHours_q = tq/3600;

subplot(2,1,2)
hold on
%plot(timeHours, SI,'k','LineWidth',2)
%plot(timeHours_q, MSI_k1,'r','LineWidth',2)
plot(timeHours_q, MSI_k2,'b','LineWidth',1.5)
plot(timeHours_q, MSI_k3,'g','LineWidth',1)
legend('PI','1min-MAvg','30min-MAvg','1hr-MAvg')

title('Solar Income')
xlabel('Time of Day (hh:mm)')
ylabel('Power Income (W)')

% Format x-axis ticks
xticks(0:0.5:24)
xticklabels(datestr(xticks/24, 'HH:MM'))



% A = [0,0,0,0,0,1;
%     1,1,1,1,1,1;
%     0,0,0,0,1,0;
%     5,4,3,2,1,0;
%     0,0,0,2,0,0;
%     20,12,6,2,0,0];
% 
% C = zeros(6);
% D = zeros(6);
% 
% for i = 1:6
%     b = zeros(6,1);
%     b(i) = 1;
%     C(i,:) = inv(A)*b;
%     for j = 1:5
%         D(i,j+1) = C(i,j)*(6-j);
%         E(i,j+1) = D(i,j)*(6-j);
%     end
% end
% 
% 
% 
% spl = @(s,Pi,Pj,Ti,Tj,Mi,Mj)   C*[s.^5;s.^4;s.^3;s.^2;s;ones(1,length(s))].*[Pi;Pj;Ti;Tj;Mi;Mj];
% 
% spl(0:.1:1,1,1,2,3,4,2);
% 
% d1 = @(s,Pi,Pj,Ti,Tj,Mi,Mj)   dot(D*[s.^5;s.^4;s.^3;s.^2;s;ones(1,length(s))],[Pi;Pj;Ti;Tj;Mi;Mj]);
% 
% d1(0:1,1,1,2,3,4,2)
% 
% for j = 1:3
%     DS(:,j) = d1(0,[1:1],[1:2],2,3,4,2)
% end


