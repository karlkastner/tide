% 2017-11-02 11:22:14.144945274 +0100 test_utm2latlon.m

lat = -90  + 180*rand();
lon = -180 + 360*rand();
%^zone = '49M';
zone = []

[x y zone]          = latlon2utm(lat,lon,zone)
[lat(:,2) lon(:,2)] = utm2latlon(x,y,zone)

int32(lat)
int32(lon)

