function [R,J,AVC,FS,L] = Hermite_Path5(nodes,tangents,seconds,mode,points,trg_vec)
% takes in a set of nodes and tangents, generates a spatial curve R
%
% Nodes expressed in terms of an 3xN matrix consisting of column vectors
% encoding the node position.
%
% Tangents in terms of a 3xN matrix consisting of column vectors encoding
% the forward facing tangent vectors
%
% Mode == 1 corresponds to closed pathes, Mode == 0 corresponds to open
% paths
%
% Points stipulates the number of points provided in the output curve per
% segment

if mode == 1
    Pi = nodes;
    Pj = [nodes(:,2:end),nodes(:,1)];
    Ti = tangents;
    Tj = [tangents(:,2:end),tangents(:,1)];
    Mi = seconds;
    Mj = [seconds(:,2:end),seconds(:,1)];
    nSegments = size(nodes,2);
elseif mode == 0
    Pi = nodes(:,1:end-1);
    Pj = nodes(:,2:end);
    Ti = tangents(:,1:end-1);
    Tj = tangents(:,2:end);
    Mi = seconds(:,1:end-1);
    Mj = seconds(:,2:end);
    nSegments = size(nodes,2)-1;
end


% 5th order hermite spline coefficient matrix

C = [ -6.0000,   15.0000,  -10.0000,         0,         0,    1.0000;
    6.0000,  -15.0000,   10.0000,         0,         0,         0;
   -3.0000,    8.0000,   -6.0000,         0,    1.0000,         0;
   -3.0000,    7.0000,   -4.0000,         0,         0,         0;
   -0.5000,    1.5000,   -1.5000,    0.5000,         0,         0;
    0.5000,   -1.0000,    0.5000,         0,         0,         0];

% first derivative coefficient matrix

D  = [   0,  -30.0000,   60.0000,  -30.0000,         0,         0;
         0,   30.0000,  -60.0000,   30.0000,         0,         0;
         0,  -15.0000,   32.0000,  -18.0000,         0,    1.0000;
         0,  -15.0000,   28.0000,  -12.0000,         0,         0;
         0,   -2.5000,    6.0000,   -4.5000,    1.0000,         0;
         0,    2.5000,   -4.0000,    1.5000,         0,         0];

%  second derivative coefficient matrix

E = [    0,         0, -120.0000,  180.0000,  -60.0000,         0;
         0,         0,  120.0000, -180.0000,   60.0000,         0;
         0,         0,  -60.0000,   96.0000,  -36.0000,         0;
         0,         0,  -60.0000,   84.0000,  -24.0000,         0;
         0,         0,  -10.0000,   18.0000,   -9.0000,    1.0000;
         0,         0,   10.0000,  -12.0000,    3.0000,         0];





% spline generator -- s spans 0 to 1 
spl = @(s,Pi,Pj,Ti,Tj,Mi,Mj)  dot(C*[s.^5;s.^4;s.^3;s.^2;s;ones(1,length(s))],[Pi;Pj;Ti;Tj;Mi;Mj]);  
% first derivative
d1 = @(s,Pi,Pj,Ti,Tj,Mi,Mj)   dot(D*[s.^5;s.^4;s.^3;s.^2;s;ones(1,length(s))],[Pi;Pj;Ti;Tj;Mi;Mj]);
% second derivative
d2 = @(s,Pi,Pj,Ti,Tj,Mi,Mj)   dot(E*[s.^5;s.^4;s.^3;s.^2;s;ones(1,length(s))],[Pi;Pj;Ti;Tj;Mi;Mj]);    

% curvature along segment - in m^-1
kappa = @(d1,d2) norm(cross( d1' , d2'  ) ) / norm( d1 )^3 ;
% remember, curvature is 1/R, where R is the radius of the obfuscating
% circle


% defining the curve
R  = zeros(6,nSegments*points);
FS = zeros(9,nSegments*points);
i = 1;
for k = 1:nSegments
    for s = linspace(0,1,points)
        spline = zeros(3,1);
        Ds = zeros(3,1);
        DDs = zeros(3,1);
        for j = 1:3
            spline(j,1) = spl(s,Pi(j,k),Pj(j,k),Ti(j,k),Tj(j,k),Mi(j,k),Mj(j,k));
            Ds(j,1) =     d1(s,Pi(j,k),Pj(j,k),Ti(j,k),Tj(j,k),Mi(j,k),Mj(j,k));
            DDs(j,1) =    d2(s,Pi(j,k),Pj(j,k),Ti(j,k),Tj(j,k),Mi(j,k),Mj(j,k));
        end

        curv = kappa(Ds,DDs);

        R(1:3,i) = spline;
        R(4:6,i) = Ds;
        R(7:9,i) = DDs;
        R(10,i) = curv;
        R(11,i) = norm(Ds);
        R(12,i) = -dot(Ds,trg_vec);


        % Frenet - Serret Ref. Frame

        % T vector
        FS(1:3,i) = Ds ./ norm(Ds) ;
        % N vector
        T_hat = FS(1:3,i);
        N_raw = DDs - dot(DDs, T_hat) * T_hat;
        FS(4:6,i) = N_raw ./ norm(N_raw);
        % B vector
        FS(7:9,i) = cross(FS(1:3,i),FS(4:6,i));
        
        i = i+1; 
    end
end


% finding the average curvature
delta_s = 1 / (points - 1);
L = sum(R(11,:)*delta_s);
AVC = sum(R(10,:) .* R(11,:)) * delta_s / L;
TRG = sum(R(12,:) .* R(11,:)) * delta_s / L;

C_cost = sum((R(10,:).^2).*R(11,:));


w1 = 1;
w2 = 20;
w3 = 0.8;

T1 = w1*C_cost;
T2 = w2*TRG^2;
T3 = w3*L;
J = T1+T2+T3;
