function [ distance ] = ppdistanceNZ(lat1,lon1,dep1,lat2,lon2,dep2 )
% This function takes coordinates and depths of two points and returns the
% distance between them

% Haversine formula for the angle between two coordinates
% http://www.movable-type.co.uk/scripts/latlong.html
a=(sin(pi*(lat2-lat1)/360))^2+cos(pi*lat1/180)*cos(pi*lat2/180)*(sin(pi*(lon2-lon1)/360))^2;
haver=2*atan2(sqrt(a),sqrt(1-a));

R=6371;  % Earth's mean radius in km

%Law of cosines
distance=sqrt((R-dep2)^2+(R-dep1)^2-(2*(R-dep2)*(R-dep1)*cos(haver)));

end

