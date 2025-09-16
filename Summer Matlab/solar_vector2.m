function SI  = solar_vector2(time_info,lat,long,alt)

    % Convert to radians
    lat = deg2rad(lat);
    long = deg2rad(long);

    % Extract time
    YYYY = time_info(1);
    MM   = time_info(2);
    DD   = time_info(3);
    HH   = time_info(4);
    Min  = time_info(5);
    Sec  = time_info(6);

    %% --- Compute "time of year" fraction since winter solstice ---
    % Days from Jan 1 to current date
    current_days = month2days(MM,YYYY) + DD - 1 ...
                   + (HH + (Min + Sec/60)/60) / 24;
    % Days from Jan 1 to winter solstice in previous year (Dec 21)
    solstice_prev = month2days(12,YYYY-1) + 21 - 1;

    % Total days in previous year
    days_prev_year = 365 + leapcomp(YYYY-1);

    % Days since last solstice
    days_since_solstice = current_days + (days_prev_year - solstice_prev);

    % Total days in current year
    days_this_year = 365 + leapcomp(YYYY);

    % Fraction of year since solstice
    ToY = days_since_solstice / days_this_year;

    % Orbital angle
    angle_in_orbit = -(2*pi)*(ToY);

    %% --- Geometry ---
    % Earth radius with altitude
    R_e = 6.371e6 + alt;

    % Solar vector in orbital frame
    S_vec = [cos(angle_in_orbit), sin(angle_in_orbit), 0];
    S_vec = S_vec / norm(S_vec);

    % Position on Earth surface
    P = R_e * [cos(lat)*cos(long), cos(lat)*sin(long), sin(lat)];

    % Earth axis tilt
    tilt = deg2rad(-23.5);
    a = [sin(tilt), 0, cos(tilt)];
    north_pole = a';

    q1 = quaternion([cos(tilt/2), 0, sin(tilt/2), 0]);

    %% --- Time of day rotation ---
    % Fraction of current day
    frac_day = (HH + Min/60 + Sec/3600) / 24;

    % This matches: theta = 2*pi*days(time_info-ref_D)
    theta = 2*pi*frac_day + angle_in_orbit + pi;

    % Rotations
    q2 = quaternion([cos(theta/2), sin(theta/2)*a]);
    q_total = q2 * q1;

    % Position in ECEF
    r = rotatepoint(q_total, P)';

    % Up vector
    up = r / norm(r);

    % North vector
    north = north_pole - (dot(north_pole, up) * up);
    north = north / norm(north);

    % East vector
    east = cross(north, up);

    % Solar vector components in NEU
    SI = [dot(S_vec, north) / norm(north);
          dot(S_vec, east)  / norm(east);
          dot(S_vec, up)    / norm(up)];

    SI = -SI;
end


%% Helper: cumulative days until given month
function days = month2days(months,year)
    days_in_month = [31,28,31,30,31,30,31,31,30,31,30,31];
    if leapcomp(year)
        days_in_month(2) = 29;
    end
    days = sum(days_in_month(1:months-1));
end

%% Helper: leap year check
function comp = leapcomp(year)
    comp = (mod(year,4) == 0 && (mod(year,100) ~= 0 || mod(year,400) == 0));
end
