function [R,J,AVC,FS] = Hermite_Path(nodes,tangents,mode,points,trg_vec)
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
    Mi = tangents;
    Mj = [tangents(:,2:end),tangents(:,1)];
    nSegments = size(nodes,2);
elseif mode == 0
    Pi = nodes(:,1:end-1);
    Pj = nodes(:,2:end);
    Mi = tangents(:,1:end-1);
    Mj = tangents(:,2:end);
    nSegments = size(nodes,2)-1;
end

% spline generator -- s spans 0 to 1 
spl = @(s,Pi,Pj,Mi,Mj)      (2*s^3-3*s^2+1)*Pi + ...
                            (s^3-2*s^2+s)*Mi + ...
                            (-2*s^3+3*s^2)*Pj + ... 
                            (s^3-s^2)*Mj;
% first derivative
d1 = @(s,Pi,Pj,Mi,Mj)       (6*s^2-6*s)*Pi + ...
                            (3*s^2-4*s+1)*Mi + ...
                            (-6*s^2+6*s)*Pj + ... 
                            (3*s^2-2*s)*Mj;
% second derivative
d2 = @(s,Pi,Pj,Mi,Mj)       (12*s-6)*Pi + ...
                            (6*s-4)*Mi + ...
                            (-12*s+6)*Pj + ... 
                            (6*s-2)*Mj;
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
        
        Ds = d1(s,Pi(:,k),Pj(:,k),Mi(:,k),Mj(:,k));
        DDs = d2(s,Pi(:,k),Pj(:,k),Mi(:,k),Mj(:,k));
        curv = kappa(Ds,DDs);


        R(1:3,i) = spl(s,Pi(:,k),Pj(:,k),Mi(:,k),Mj(:,k));
        R(4,i) = curv;
        R(5,i) = norm(Ds);
        R(6,i) = -dot(Ds,trg_vec);


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
L = sum(rmmissing(R(end-1,:)));
AVC = sum((rmmissing(R(end-2,:)).*R(end-1,:)))/L;
w1 = 1;
w2 = 1;
J = sum(   w1*(rmmissing(R(end-2,:).^2).*R(end-1,:)) +  w2*(rmmissing(R(end,:).^2).*R(end-1,:))   )/L;