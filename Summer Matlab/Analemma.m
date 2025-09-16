function Analemma()
holdT = clock;
yr = holdT(1);
n = 1;
for m=1:12
    for d=1:eomday(yr,m)
        dt(n) = datetime(yr,m,d,12,0,0);
        n = n+1;
    end
end
% solar position for each day
dt.TimeZone = 'Etc/GMT';
[declin,~,~,eqTime] = EarthEphemeris(dt);
plot(eqTime,declin,'LineWidth',1)
xlabel('equation of time, min')
ylabel('declination, deg')
% plot month names along analemma
hold on;
for m=1:12
    thisDate = datetime([yr,m,15,12,0,0]);
    thisDate.TimeZone = dt.TimeZone;
    k = dt==thisDate;
    ds = datestr(thisDate,'mmm');
    text(eqTime(k),declin(k),ds)
end
% explain
fprintf('declination is the latitude at which the solar zenith angle is zero\n')
fprintf('equation of time is solar time minus clock time (in January, the sun seems to linger near sunset)\n')
end