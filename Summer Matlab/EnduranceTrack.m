clear
clc
close all
%% Input 
    % Number of radial spokes
    num_spokes = 10;
    % Bounding polygon, specified counter-clockwise
    poly_bound_x = 10*[0,50,100,80,65,50] ;
    poly_bound_y = 10*[0,-50,-20,200,150,80];
%% Initial Bounding Box Processing

    box_geom = polyshape(poly_bound_x,poly_bound_y);
    [hub_x,hub_y] = centroid(box_geom);
    disp('Bounding Box Centroid:')
    disp([hub_x,hub_y])
    
    figure
    hold on
    plot(box_geom)
    plot(hub_x,hub_y,'xk',MarkerSize=10)
    title('Uncentered Bounding Polygon with Centroid')
    xlabel('x [m]')
    ylabel('y [m]')
    
    
    poly_bound_x = poly_bound_x - hub_x;
    poly_bound_y = poly_bound_y - hub_y;
    box_geom = polyshape(poly_bound_x,poly_bound_y);

    % find max spoke extension 
    radials = linspace(0,2*pi,num_spokes);
    [vertex_t,vertex_r] = cart2pol(poly_bound_x,poly_bound_y);

    for i = 1:num_spokes
        for j = 1:length(vertex_t)-1
            if (radials(i) > vertex_t(j)) && (radials(i) < vertex_t(j+1))
                rad_vec = [cos( radials(i) ); sin( radials(i) )];
            end
        end
    end
    
    figure
    hold on
    plot(box_geom)
    plot(0,0,'xk',MarkerSize=10)
    title('Centered Bounding Polygon with Centroid')
    xlabel('x [m]')
    ylabel('y [m]')
    
    

