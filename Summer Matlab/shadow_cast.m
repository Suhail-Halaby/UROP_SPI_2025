function [combined_EF,area_fraction,sun_component,X,Y,Z] = shadow_cast(plane,elv,azimuth,mode)
    elv = - elv;
    [X,Y,Z] = sph2cart(azimuth,elv,1);
    Sun_Vec = -[X;Y;Z];
     
    area_fraction = zeros(1,size(plane.solar,1));
    sun_component = zeros(1,size(plane.solar,1));


    for i = 1 :size(plane.solar,1)
        Vrts = plane.solar{i,6};
        X_r = Vrts(1,:);
        Y_r = Vrts(2,:);
        Z_r = Vrts(3,:);

        % find the normal vector of each cell -- assuming the cells are
        % planar quadrilaterals
        

        Tab = [X_r(1),Y_r(1),Z_r(1)] - [X_r(2),Y_r(2),Z_r(2)];
        Tac = [X_r(1),Y_r(1),Z_r(1)] - [X_r(4),Y_r(4),Z_r(4)];

        normal_cell_vector = cross(Tab,Tac);
        normal_cell_vector = normal_cell_vector ./ norm(normal_cell_vector);
        
    
        scell = zeros(2,4);
           
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
        sun_component(i) = max(0,dot(normal_cell_vector',-Sun_Vec));
        area_fraction(i) = area(trg_panel)/area(solar_panel);
        clear  G H X_r Y_r Z_r j v_cell h_cell wing_cell_temp  scell Vrts Vrts_w X_w Y_w Z_w p K_list K cut_out
    end
    
    if mode == 1
        combined_EF = sun_component.*area_fraction;
    else 
        combined_EF = area_fraction;
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























end