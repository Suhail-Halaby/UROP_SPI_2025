clear
clc
close all

% Vertices of each of the panels [X;Y;Z]

%% Main Wing
n_elem = 3; 

chord_dist = [0.25, 0.25, 0.20,0.10];
span = [1.9,0.1,0.1];
dihedral   = [6,8,12];
twist      = [3,3,3,3];
sweep      = [0,0,0.05,0.1];

x_root = 0.75;
z_root = 0;

symmetric = true;

plane.wing = initialize_panel(n_elem,chord_dist,span,dihedral,twist,sweep,x_root,z_root,symmetric);
clear n_elem chord_dist span dihedral twist sweep x_root z_root symmetric
%% H-stab
n_elem = 1; 

chord_dist = [0.25, 0.25];
span = 0.5;
dihedral   = [0,0];
twist      = [0,0];
sweep      = [0,0];

x_root = 1.90;
z_root = 0.4;

symmetric = true;

plane.hstab = initialize_panel(n_elem,chord_dist,span,dihedral,twist,sweep,x_root,z_root,symmetric);
clear n_elem chord_dist span dihedral twist sweep x_root z_root symmetric

%% V-stab
n_elem = 1; 

chord_dist = [0.3, 0.3];
span       = 0.4;
dihedral   = [90,0];
twist      = [0,0];
sweep      = [0,0];

x_root = 1.90;
z_root = 0;

symmetric = false;

plane.vstab = initialize_panel(n_elem,chord_dist,span,dihedral,twist,sweep,x_root,z_root,symmetric);
clear n_elem chord_dist span dihedral twist sweep x_root z_root symmetric

%% SOLAR Cell Generators

n_elem = 12; 

chord_dist = 0.152*ones(1,n_elem+1);
span       = 0.152*ones(1,n_elem);
dihedral   = 6*ones(1,n_elem);
twist      = 3*ones(1,n_elem+1);
sweep      =  zeros(1,n_elem+1);

x_root = 0.8;
z_root = 0.01;

symmetric = true;

plane.solar = initialize_panel(n_elem,chord_dist,span,dihedral,twist,sweep,x_root,z_root,symmetric);
clear n_elem chord_dist span dihedral twist sweep x_root z_root symmetric



%% Shadow Casting
%hold on
%axis equal
elv = deg2rad(-15);  
azimuth = deg2rad(175); 

[X,Y,Z] = sph2cart(azimuth,elv,1);
Sun_Vec = -[X;Y;Z];
 
area_fraction = zeros(1,size(plane.solar,1));

for i = 1 :size(plane.solar,1)
    Vrts = plane.solar{i,6};
    X_r = Vrts(1,:);
    Y_r = Vrts(2,:);
    Z_r = Vrts(3,:);

    scell = zeros(2,4);
    wing_cell = zeros(2*size(plane.wing,1),4);


    for p = 1:4
       [H,G] = project_to_plane(X_r,Y_r,Z_r,X_r(p),Y_r(p),Z_r(p),Sun_Vec);
       scell(:,p) = [H;G];
    end
    
    
    trg_panel = polyshape(scell(1,:),scell(2,:));
    solar_panel = trg_panel;
    %plot(trg_panel)
    clear H G

    for j = 1:size(plane.wing,1)
        Vrts_w = plane.wing{j,6};
        X_w = Vrts_w(1,:);
        Y_w = Vrts_w(2,:);
        Z_w = Vrts_w(3,:);
        wing_cell_temp  = zeros(2,4);
        K_list = zeros(1,4);

        for p = 1:4
            [H,G,K] = project_to_plane(X_r,Y_r,Z_r,X_w(p),Y_w(p),Z_w(p),Sun_Vec);
            wing_cell_temp (:,p) = [H;G];
            K_list(1,p) = K;
        end
        
        if max(K_list) >= 0
            w_panel = polyshape(wing_cell_temp (1,:),wing_cell_temp (2,:));
            cut_out = intersect(trg_panel,w_panel);
            trg_panel = subtract(trg_panel,cut_out);
            %plot(w_panel)
        end

    end 

    clear H G K

    for j = 1:size(plane.hstab,1)
        Vrts_w = plane.hstab{j,6};
        X_w = Vrts_w(1,:);
        Y_w = Vrts_w(2,:);
        Z_w = Vrts_w(3,:);
        h_cell = zeros(2,4);
        for p = 1:4
            [H,G,K] = project_to_plane(X_r,Y_r,Z_r,X_w(p),Y_w(p),Z_w(p),Sun_Vec);
            h_cell(:,p) = [H;G];
            K_list(:,p) = K;
        end

        if max(K_list) >= 0
            h_panel = polyshape(h_cell(1,:),h_cell(2,:));
            cut_out = intersect(trg_panel,h_panel);
            trg_panel = subtract(trg_panel,cut_out);
            %plot(h_panel)
        end

    end 

    clear H G 

    for j = 1:size(plane.vstab,1)
        Vrts_w = plane.vstab{j,6};
        X_w = Vrts_w(1,:);
        Y_w = Vrts_w(2,:);
        Z_w = Vrts_w(3,:);
        v_cell = zeros(2,4);
        for p = 1:4
            [H,G] = project_to_plane(X_r,Y_r,Z_r,X_w(p),Y_w(p),Z_w(p),Sun_Vec);
            v_cell(:,p) = [H;G];
        end
        v_panel = polyshape(v_cell(1,:),v_cell(2,:));
        %plot(v_panel)

        cut_out = intersect(trg_panel,v_panel);
        trg_panel = subtract(trg_panel,cut_out);
    end 
    
    %plot(trg_panel)

    area_fraction(i) = area(trg_panel)/area(solar_panel);
    clear  G H X_r Y_r Z_r j v_cell h_cell wing_cell_temp  scell Vrts Vrts_w X_w Y_w Z_w p K_list K cut_out
end




%clear h_panel v_panel w_panel trg_panel

%% Plotting
figure
hold on
axis equal
grid on
xlabel('X')
ylabel('Y')
zlabel('Z')
view(3)

%quiver3(0,0,0,Sun_Vec(1),Sun_Vec(2),Sun_Vec(3),2,'magenta',LineWidth=2)


for i = 1:size(plane.wing,1)
    c = plane.wing{i,6};
    fill3(c(1,:), c(2,:), c(3,:), 'cyan', 'FaceAlpha', 0.5, 'EdgeColor', 'k')    
    quiver3(c(1,:), c(2,:), c(3,:),X*ones(1,4),Y*ones(1,4),Z*ones(1,4),3,'k')
end
clear i 
for i = 1:size(plane.hstab,1)
    c = plane.hstab{i,6};
    fill3(c(1,:), c(2,:), c(3,:), 'cyan', 'FaceAlpha', 0.5, 'EdgeColor', 'k')
    quiver3(c(1,:), c(2,:), c(3,:),X*ones(1,4),Y*ones(1,4),Z*ones(1,4),3,'k')
end
clear i
for i = 1:size(plane.vstab,1)
    c = plane.vstab{i,6};
    fill3(c(1,:), c(2,:), c(3,:), 'cyan', 'FaceAlpha', 0.5, 'EdgeColor', 'k')
    quiver3(c(1,:), c(2,:), c(3,:),X*ones(1,4),Y*ones(1,4),Z*ones(1,4),3,'k')
end
clear i 
for i = 1:size(plane.solar,1)
    c = plane.solar{i,6};
    fill3(c(1,:), c(2,:), c(3,:), 'blue', 'FaceAlpha', 0.5, 'EdgeColor', 'k')
    %quiver3(c(1,:), c(2,:), c(3,:),X*ones(1,4),Y*ones(1,4),Z*ones(1,4),5,'--k')
end
clear i 
clear c X Y Z



%% PANEL GENERATION FUNCTION


function object = initialize_panel(n_elem,chord_dist,span,dihedral,twist,sweep,x_root,z_root,symmetric)
    if symmetric == true
        object = cell(2*n_elem,6);
    else
        object = cell(n_elem,6);
    end 

    for i = 1:n_elem
        object{i,1} = chord_dist(i:i+1);          % vector
        object{i,2} = span(i); 
        object{i,3} = deg2rad(dihedral(i));       % scalar
        object{i,4} = twist(i:i+1);               % vector
        object{i,5} = sweep(i:i+1);               % vector
        
    
        % corner order: front inner --> front outer --> back outer --> back
        % inner
        if i == 1
            
            cornersx = x_root+[0,sweep(i+1),sweep(i+1)+chord_dist(i+1),chord_dist(i)];
            cornersy = [0,span(i)*cosd(dihedral(i)),span(i)*cosd(dihedral(i)),0];
            cornersz = z_root+[0,span(i)*sind(dihedral(i)),span(i)*sind(dihedral(i)),0];
    
            object{i,6} = [cornersx;cornersy;cornersz];
        else 
            cornersx = cornersx(2)+[0,sweep(i+1),sweep(i+1)+chord_dist(i+1),chord_dist(i)];
            cornersy = cornersy(2)+[0,span(i)*cosd(dihedral(i)),span(i)*cosd(dihedral(i)),0];
            cornersz = cornersz(2)+[0,span(i)*sind(dihedral(i)),span(i)*sind(dihedral(i)),0];
    
            object{i,6} = [cornersx;cornersy;cornersz];
        end
    end
    
    % Symmetric Duplication 
    clear i
    if symmetric == true
        for i = n_elem+1 : 2*n_elem
            object{i,1} = object{i-n_elem,1};
            object{i,2} = object{i-n_elem,2};
            object{i,3} = object{i-n_elem,3};
            object{i,4} = object{i-n_elem,4};
            object{i,5} = object{i-n_elem,5};
        
        
            coords = object{i-n_elem,6};
        
            object{i,6} = [coords(1,:);-coords(2,:);coords(3,:)];
        end
    end

    clear i cornersz cornersy cornersx

end


%% Projection function

function [X_proj,Y_proj,K] = project_to_plane(X_r,Y_r,Z_r,X,Y,Z,vector)
    % Define basis from reference points
    i_hat = [X_r(1)-X_r(4); Y_r(1)-Y_r(4); Z_r(1)-Z_r(4)];
    j_hat = [X_r(1)-X_r(2); Y_r(1)-Y_r(2); Z_r(1)-Z_r(2)];

    % Gram–Schmidt to ensure orthonormal basis
    i_hat = i_hat / norm(i_hat);
    j_hat = j_hat - dot(j_hat,i_hat)*i_hat;
    j_hat = j_hat / norm(j_hat);
    k_hat = cross(i_hat,j_hat);

    % Translate point into local frame
    origin = [X_r(1); Y_r(1); Z_r(1)];
    rel_P  = [X;Y;Z] - origin;

    % Local coordinates of point
    X_f = dot(i_hat,rel_P);
    Y_f = dot(j_hat,rel_P);
    Z_f = dot(k_hat,rel_P);

    % Local coordinates of projection vector
    V_i = dot(i_hat,vector);
    V_j = dot(j_hat,vector);
    V_k = dot(k_hat,vector);

    if abs(V_k) < 1e-12
        error('Vector is parallel to plane, no intersection.');
    end

    % Intersection with plane Z_f = 0
    K = Z_f / V_k;
    proj_local = [X_f, Y_f, Z_f] - K * [V_i, V_j, V_k];

    X_proj = proj_local(1);
    Y_proj = proj_local(2);
end
