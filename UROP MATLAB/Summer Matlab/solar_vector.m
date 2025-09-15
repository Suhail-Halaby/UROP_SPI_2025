function [SI] = solar_vector(time_info,lat,long,alt)
    
    lat = deg2rad(lat);
    long = deg2rad(long);


    % Year:
    YYYY = year(time_info);
    % Month
    MM = month(time_info);
    % Day:
    DD = day(time_info);
    % Hour
    HH = hour(time_info);
    % Minute 
    Min = minute(time_info);
    % Second
    Sec = second(time_info);

    % Computing time of year
    ref_Y = datetime(YYYY-1,12,21,0,0,0);
    ref_D = datetime(YYYY,MM,DD,0,0,0);
    ToY = years(time_info - ref_Y);

    angle_in_orbit = -(2*pi)*(ToY);
   
    % Earth radius and added altitude (m)
    R_e = 6.371e6 + alt;
    
    % Solar vector in orbital frame
    S_vec = [cos(angle_in_orbit),sin(angle_in_orbit),0];
    S_vec = S_vec / norm(S_vec);
    

    % position point on earth's surface
    P = R_e*[cos(lat)*cos(long),cos(lat)*sin(long),sin(lat)];
    P_quat = quaternion([0,P]);

    % local NEU frame
    north = [-sin(lat)*cos(long),-sin(lat)*sin(long),cos(lat)]; 
    east  = [-sin(long),cos(long),0];
    up    = [cos(lat)*cos(long),cos(lat)*sin(long),sin(lat)];
    frame = [north;east;up];

    % Earth axis of rotation
    tilt = deg2rad(-23.5);
    a = [sin(tilt),0,cos(tilt)];
    north_pole = a';

    q1 = quaternion([cos(tilt/2),0,sin(tilt/2),0]);
 
    % time of day rotational position:

    theta = 2*pi*days(time_info-ref_D);
    theta = theta + angle_in_orbit + pi;

    P_rot = zeros(3,1);
    F_rot = zeros(3,3,1);
    SI = zeros(3,1);

    % Canonical reference vectors
    P0 = [0;0;1]; % North pole
    U0 = [0;0;1]; % also used for Up calculation (radial from origin)

    q2 = quaternion([cos(theta/2), sin(theta/2)*a]);
    q_total = q2 * q1;  % apply q1 then q2

    % Position in ECEF
    r = rotatepoint(q_total, P)';     % 3x1 position vector
    P_rot(:,1) = r;

    % Up = normalized radial
    up = r / norm(r);

    % North = projection of rotated pole onto tangent plane
    north = north_pole - (dot(north_pole, up) * up);
    north = north / norm(north);

    % East = cross(north, up)
    east = cross(north, up);

    % Store NEU as rows
    F_rot(:,:,1) = [north'; east'; up'];

    % Sun dot prod
    SI(:,1) = [dot(S_vec,north) / norm(north) ;dot(S_vec,east) / norm(east);dot(S_vec,up) / norm(up)];
    
end