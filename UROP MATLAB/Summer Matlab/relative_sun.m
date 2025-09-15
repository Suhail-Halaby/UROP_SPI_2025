function [relative_elv, relative_azimuth,vec] = relative_sun(sun_vector,Euler)


% Sun vector points INTO the aircraft
% It points down when the sun is high
% North at noon in the north hemisphere
% West at sunrise

% Taking the negative (pointing out)

sun_vector = -sun_vector;

% Reformatting into global NED

sun_vector(3) = -sun_vector(3);


% Reframing the sun vector into the body fixed system


pitch = Euler(1);
roll = Euler(2);
heading = Euler(3);


R1E = [ cos(heading),     sin(heading),  0; ...   
        -sin(heading),    cos(heading),  0; ...
        0,                0,             1];

R21 = [ cos(pitch),     0,     -sin(pitch); ...
        0,              1,      0;          ...
        sin(pitch),     0,      cos(pitch)];

RB2 = [ 1,               0,               0; ...
        0,       cos(roll),       sin(roll); ...
        0,      -sin(roll),       cos(roll)];


sun_body_fixed = (RB2*R21*R1E)*sun_vector;
vec = sun_body_fixed;

sx = sun_body_fixed(1);
sy = sun_body_fixed(2);
sz = sun_body_fixed(3);

relative_azimuth = atan2(sy, sx);         % azimuth in xy-plane
relative_elv     = atan2(-sz, hypot(sx,sy)); % use -sz because body z is down




end