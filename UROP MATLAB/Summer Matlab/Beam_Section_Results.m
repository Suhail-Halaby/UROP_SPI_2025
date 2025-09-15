clear
clc
close all



%% Wing 

 
% Aluminium Spar -- Outer Box
x_offs = 46;
y_offs = 5;
spar1_x = [10 10 -10 -10 ]+x_offs; 
spar1_y = [-10 10 10 -10 ]+y_offs;
[Cxso, Cyso, Ixxso, Iyyso, Ixyso, Aso] = polygon_properties(spar1_x/1000,spar1_y/1000);

% Aluminium Spar -- Inner Box
spar1_x = [8.4 8.4 -8.4 -8.4 ]+x_offs; 
spar1_y = [-8.4 8.4 8.4 -8.4 ]+y_offs;
[Cxsi, Cysi, Ixxsi, Iyysi, Ixysi,Asi] = polygon_properties(spar1_x/1000,spar1_y/1000);

% total spar properties
E_Al = 70e9;
G_Al = E_Al / 1.4;
rho_Al = 2700;
Al = [Cxso, Cyso, Ixxso, Iyyso, Ixyso, Aso]- [0, 0, Ixxsi, Iyysi, Ixysi, Asi];


% CFRP Aft Spar -- Outer Ring
radius = 4 ;  
numPoints = 100; 
theta = linspace(0, 2*pi, numPoints);

x_offc = 176;
y_offc = 3;

x_co = x_offc + (radius * cos(theta));
y_co = y_offc + (radius * sin(theta));

[Cxco, Cyco, Ixxco, Iyyco, Ixyco, Aco] = polygon_properties(x_co/1000,y_co/1000);

% CFRP Aft Spar -- Inner Ring
radius = 3.5 ;
x_ci = x_offc + (radius * cos(theta));
y_ci = y_offc + (radius * sin(theta));

[Cxci, Cyci, Ixxci, Iyyci, Ixyci, Aci] = polygon_properties(x_ci/1000,y_ci/1000);

% total spar properties
E_CF = 70e9;
G_CF = E_CF / 1.4;
rho_CF = 1700;
CF = [Cxco, Cyco, Ixxco, Iyyco, Ixyco, Aco] - [0, 0, Ixxci, Iyyci, Ixyci, Aci];


% foam wing section
chord = 250;
CC = readmatrix('naca4415.dat');
x_af = chord*CC(2:end,1);
y_af = chord*CC(2:end,2);

[Cxa, Cya, Ixxa, Iyya, Ixya, Aa] = polygon_properties(x_af/1000,y_af/1000);

% total spar properties
E_F = 5e6;
G_F = E_F / 1.4;
rho_f = 60;
Af = [Cxa, Cya, Ixxa, Iyya, Ixya, Aa] - CF - Al;


%Af = zeros(1,6)
%CF = zeros(1,6)

% processing

Cx_total = (rho_f*Af(6)*Af(1) + rho_Al*Al(6)*Al(1) + rho_CF*CF(6)*CF(1)) / ...
           (rho_f*Af(6) + rho_Al*Al(6) + rho_CF*CF(6));

Cy_total = (rho_f*Af(6)*Af(2) + rho_Al*Al(6)*Al(2) + rho_CF*CF(6)*CF(2)) / ...
           (rho_f*Af(6) + rho_Al*Al(6) + rho_CF*CF(6));



% Offsets from total centroid
dx_Al = Al(1) - Cx_total;
dy_Al = Al(2) - Cy_total;
dx_CF = CF(1) - Cx_total;
dy_CF = CF(2) - Cy_total;
dx_F = Af(1) - Cx_total;
dy_F = Af(2) - Cy_total;

% Apply parallel axis theorem for each component
% Aluminium
Ixx_Al = Al(3) + Al(6) * dy_Al^2;
Iyy_Al = Al(4) + Al(6) * dx_Al^2;
Izz_Al = Ixx_Al + Iyy_Al;

% CFRP
Ixx_CF = CF(3) + CF(6) * dy_CF^2;
Iyy_CF = CF(4) + CF(6) * dx_CF^2;
Izz_CF = Ixx_CF + Iyy_CF;

% Foam
Ixx_F = Af(3) + Af(6) * dy_F^2;
Iyy_F = Af(4) + Af(6) * dx_F^2;
Izz_F = Ixx_F + Iyy_F;

% Sectional stiffnesses (now correctly about the section centroid)
EIx = E_Al * Ixx_Al + E_CF * Ixx_CF + E_F * Ixx_F;
EIy = E_Al * Iyy_Al + E_CF * Iyy_CF + E_F * Iyy_F;
GJ  = G_Al * Izz_Al + G_CF * Izz_CF + G_F * Izz_F;


m_bar = rho_f*Af(6)+rho_Al*Al(6)+rho_CF*CF(6);  % sectional density (kg/m)
J_bar = rho_f*(sum(Af(3:4)))+rho_CF*(sum(CF(3:4)))+rho_Al*(sum(Al(3:4))); % sectional polar inertial moment

EA = E_F*Af(6)+E_Al*Al(6)+E_CF*CF(6); % EA (kgm)
GA = G_F*Af(6)+G_Al*Al(6)+G_CF*CF(6); 

GA_total = G_F*Af(6) + G_Al*Al(6) + G_CF*CF(6);
x_e = (G_F*Af(6)*Af(1) + G_Al*Al(6)*Al(1) + G_CF*CF(6)*CF(1)) / GA_total;
y_e = (G_F*Af(6)*Af(2) + G_Al*Al(6)*Al(2) + G_CF*CF(6)*CF(2)) / GA_total;





fprintf('Centroid: (%.3f mm, %.3f mm)\n', Cx_total*1000, Cy_total*1000);
fprintf('EIx = %.3e Nm^2\n', EIx);
fprintf('EIy = %.3e Nm^2\n', EIy);
fprintf('GJ  = %.3e Nm^2\n', GJ);
fprintf('EA  = %.3e Nm^-1\n', EA);
fprintf('GA  = %.3e Nm^-1\n', GA);
fprintf('m_bar = %.3f kg/m\n', m_bar);
fprintf('J_bar = %.3e kg·m\n', J_bar);
fprintf('Elastic Axis = %.3f \n', x_e*1000 / chord );


%% H-stab




% Spar 1
radius = 5 ;  
x_offc = 35;
y_offc = 0;
x_co1 = x_offc + (radius * cos(theta));
y_co1 = y_offc + (radius * sin(theta));
[Cxco1, Cyco1, Ixxco1, Iyyco1, Ixyco1, Aco1] = polygon_properties(x_co1/1000,y_co1/1000);
radius = 4.5 ;
x_ci1 = x_offc + (radius * cos(theta));
y_ci1 = y_offc + (radius * sin(theta));
[Cxci1, Cyci1, Ixxci1, Iyyci1, Ixyci1, Aci1] = polygon_properties(x_ci1/1000,y_ci1/1000);

S1 = [Cxco1, Cyco1, Ixxco1, Iyyco1, Ixyco1, Aco1] - [0, 0, Ixxci1, Iyyci1, Ixyci1, Aci1];

% Spar 2
radius = 4 ;  
x_offc = 115;
y_offc = 0;
x_co2 = x_offc + (radius * cos(theta));
y_co2 = y_offc + (radius * sin(theta));
[Cxco2, Cyco2, Ixxco2, Iyyco2, Ixyco2, Aco2] = polygon_properties(x_co2/1000,y_co2/1000);
radius = 3.5 ;
x_ci2 = x_offc + (radius * cos(theta));
y_ci2 = y_offc + (radius * sin(theta));
[Cxci2, Cyci2, Ixxci2, Iyyci2, Ixyci2, Aci2] = polygon_properties(x_ci2/1000,y_ci2/1000);

S2 = [Cxco2, Cyco2, Ixxco2, Iyyco2, Ixyco2, Aco2] - [0, 0, Ixxci2, Iyyci2, Ixyci2, Aci2];


% Spar 3
radius = 4 ;  
x_offc = 185;
y_offc = 0;
x_co3 = x_offc + (radius * cos(theta));
y_co3 = y_offc + (radius * sin(theta));
[Cxco3, Cyco3, Ixxco3, Iyyco3, Ixyco3, Aco3] = polygon_properties(x_co3/1000,y_co3/1000);
radius = 3.5 ;
x_ci3 = x_offc + (radius * cos(theta));
y_ci3 = y_offc + (radius * sin(theta));
[Cxci3, Cyci3, Ixxci3, Iyyci3, Ixyci3, Aci3] = polygon_properties(x_ci3/1000,y_ci3/1000);

S3 = [Cxco3, Cyco3, Ixxco3, Iyyco3, Ixyco3, Aco3] - [0, 0, Ixxci3, Iyyci3, Ixyci3, Aci3];


% Tail material properties (re-using from wing)
E_CF = 70e9;
G_CF = E_CF / 1.4;
rho_CF = 1700;

% Tail component list
spars = [S1; S2; S3];

% Mass-weighted centroid
Cx_tail = sum(rho_CF * spars(:,6) .* spars(:,1)) / sum(rho_CF * spars(:,6));
Cy_tail = sum(rho_CF * spars(:,6) .* spars(:,2)) / sum(rho_CF * spars(:,6));

% Parallel Axis Theorem for each spar
Ixx_total = 0; Iyy_total = 0; Izz_total = 0;
for i = 1:size(spars,1)
    dx = spars(i,1) - Cx_tail;
    dy = spars(i,2) - Cy_tail;
    A_i = spars(i,6);
    Ixx_i = spars(i,3) + A_i * dy^2;
    Iyy_i = spars(i,4) + A_i * dx^2;
    
    Ixx_total = Ixx_total + E_CF * Ixx_i;
    Iyy_total = Iyy_total + E_CF * Iyy_i;
    Izz_total = Izz_total + G_CF * (Ixx_i + Iyy_i);  % GJ
end

% EA and GA
EA_tail = E_CF * sum(spars(:,6));
GA_tail = G_CF * sum(spars(:,6));

% Sectional mass and inertia
m_bar_tail = rho_CF * sum(spars(:,6));
J_bar_tail = rho_CF * sum(spars(:,3) + spars(:,4));  % polar moment of inertia

% Elastic Axis (use G-weighted centroid)
GAi = G_CF * spars(:,6);
x_e_tail = sum(GAi .* spars(:,1)) / sum(GAi);
y_e_tail = sum(GAi .* spars(:,2)) / sum(GAi);

% Display
fprintf('\n=== H-stab Section Properties ===\n');
fprintf('Centroid: (%.3f mm, %.3f mm)\n', Cx_tail*1000, Cy_tail*1000);
fprintf('EIx = %.3e Nm^2\n', Ixx_total);
fprintf('EIy = %.3e Nm^2\n', Iyy_total);
fprintf('GJ  = %.3e Nm^2\n', Izz_total);
fprintf('EA  = %.3e Nm^-1\n', EA_tail);
fprintf('GA  = %.3e Nm^-1\n', GA_tail);
fprintf('m_bar = %.3f kg/m\n', m_bar_tail);
fprintf('J_bar = %.3e kg·m\n', J_bar_tail);
fprintf('Elastic Axis = %.3f mm from LE \n', x_e_tail*1000);




%% V-stab





% Spar 1
radius = 5 ;  
x_offc = 60;
y_offc = 0;
x_co1 = x_offc + (radius * cos(theta));
y_co1 = y_offc + (radius * sin(theta));
[Cxco1, Cyco1, Ixxco1, Iyyco1, Ixyco1, Aco1] = polygon_properties(x_co1/1000,y_co1/1000);
radius = 4.5 ;
x_ci1 = x_offc + (radius * cos(theta));
y_ci1 = y_offc + (radius * sin(theta));
[Cxci1, Cyci1, Ixxci1, Iyyci1, Ixyci1, Aci1] = polygon_properties(x_ci1/1000,y_ci1/1000);

S1 = [Cxco1, Cyco1, Ixxco1, Iyyco1, Ixyco1, Aco1] - [0, 0, Ixxci1, Iyyci1, Ixyci1, Aci1];

% Spar 2
radius = 4 ;  
x_offc = 125;
y_offc = 0;
x_co2 = x_offc + (radius * cos(theta));
y_co2 = y_offc + (radius * sin(theta));
[Cxco2, Cyco2, Ixxco2, Iyyco2, Ixyco2, Aco2] = polygon_properties(x_co2/1000,y_co2/1000);
radius = 3.5 ;
x_ci2 = x_offc + (radius * cos(theta));
y_ci2 = y_offc + (radius * sin(theta));
[Cxci2, Cyci2, Ixxci2, Iyyci2, Ixyci2, Aci2] = polygon_properties(x_ci2/1000,y_ci2/1000);

S2 = [Cxco2, Cyco2, Ixxco2, Iyyco2, Ixyco2, Aco2] - [0, 0, Ixxci2, Iyyci2, Ixyci2, Aci2];


% Spar 3
radius = 4 ;  
x_offc = 210;
y_offc = 0;
x_co3 = x_offc + (radius * cos(theta));
y_co3 = y_offc + (radius * sin(theta));
[Cxco3, Cyco3, Ixxco3, Iyyco3, Ixyco3, Aco3] = polygon_properties(x_co3/1000,y_co3/1000);
radius = 3.5 ;
x_ci3 = x_offc + (radius * cos(theta));
y_ci3 = y_offc + (radius * sin(theta));
[Cxci3, Cyci3, Ixxci3, Iyyci3, Ixyci3, Aci3] = polygon_properties(x_ci3/1000,y_ci3/1000);

S3 = [Cxco3, Cyco3, Ixxco3, Iyyco3, Ixyco3, Aco3] - [0, 0, Ixxci3, Iyyci3, Ixyci3, Aci3];


% Tail material properties (re-using from wing)
E_CF = 70e9;
G_CF = E_CF / 1.4;
rho_CF = 1700;

% Tail component list
spars = [S1; S2; S3];

% Mass-weighted centroid
Cx_tail = sum(rho_CF * spars(:,6) .* spars(:,1)) / sum(rho_CF * spars(:,6));
Cy_tail = sum(rho_CF * spars(:,6) .* spars(:,2)) / sum(rho_CF * spars(:,6));

% Parallel Axis Theorem for each spar
Ixx_total = 0; Iyy_total = 0; Izz_total = 0;
for i = 1:size(spars,1)
    dx = spars(i,1) - Cx_tail;
    dy = spars(i,2) - Cy_tail;
    A_i = spars(i,6);
    Ixx_i = spars(i,3) + A_i * dy^2;
    Iyy_i = spars(i,4) + A_i * dx^2;
    
    Ixx_total = Ixx_total + E_CF * Ixx_i;
    Iyy_total = Iyy_total + E_CF * Iyy_i;
    Izz_total = Izz_total + G_CF * (Ixx_i + Iyy_i);  % GJ
end

% EA and GA
EA_tail = E_CF * sum(spars(:,6));
GA_tail = G_CF * sum(spars(:,6));

% Sectional mass and inertia
m_bar_tail = rho_CF * sum(spars(:,6));
J_bar_tail = rho_CF * sum(spars(:,3) + spars(:,4));  % polar moment of inertia

% Elastic Axis (use G-weighted centroid)
GAi = G_CF * spars(:,6);
x_e_tail = sum(GAi .* spars(:,1)) / sum(GAi);
y_e_tail = sum(GAi .* spars(:,2)) / sum(GAi);

% Display
fprintf('\n=== V-stab Section Properties ===\n');
fprintf('Centroid: (%.3f mm, %.3f mm)\n', Cx_tail*1000, Cy_tail*1000);
fprintf('EIx = %.3e Nm^2\n', Ixx_total);
fprintf('EIy = %.3e Nm^2\n', Iyy_total);
fprintf('GJ  = %.3e Nm^2\n', Izz_total);
fprintf('EA  = %.3e Nm^-1\n', EA_tail);
fprintf('GA  = %.3e Nm^-1\n', GA_tail);
fprintf('m_bar = %.3f kg/m\n', m_bar_tail);
fprintf('J_bar = %.3e kg·m\n', J_bar_tail);
fprintf('Elastic Axis = %.3f mm from LE \n', x_e_tail*1000);








function [Cx, Cy, Ixx, Iyy, Ixy, A] = polygon_properties(x, y)
% polygon_properties - Computes centroid, second moment of area, and area of a 2D polygon
%
% Inputs:
%   x, y - vectors of polygon vertex coordinates (must be ordered CCW)
%
% Outputs:
%   Cx, Cy - centroid coordinates
%   Ixx, Iyy, Ixy - second moments of area about the origin
%   A - area of the polygon

    % Ensure inputs are column vectors
    x = x(:);
    y = y(:);
    
    % Close the polygon if not already closed
    if x(1) ~= x(end) || y(1) ~= y(end)
        x(end+1) = x(1);
        y(end+1) = y(1);
    end

    N = length(x) - 1;  % Number of edges
    A = 0;              % Area
    Cx = 0; Cy = 0;     % Centroid
    Ixx = 0; Iyy = 0; Ixy = 0;  % Second moments

    for i = 1:N
        xi = x(i);     yi = y(i);
        xi1 = x(i+1);  yi1 = y(i+1);
        
        a = xi * yi1 - xi1 * yi;  % Cross product (2*area element)
        A = A + a;
        
        Cx = Cx + (xi + xi1) * a;
        Cy = Cy + (yi + yi1) * a;

        Ixx = Ixx + (yi^2 + yi * yi1 + yi1^2) * a;
        Iyy = Iyy + (xi^2 + xi * xi1 + xi1^2) * a;
        Ixy = Ixy + (xi * yi1 + 2 * xi * yi + 2 * xi1 * yi1 + xi1 * yi) * a;
    end

    A = A / 2;
    Cx = Cx / (6 * A);
    Cy = Cy / (6 * A);

    Ixx = Ixx / 12;
    Iyy = Iyy / 12;
    Ixy = Ixy / 24;
end
