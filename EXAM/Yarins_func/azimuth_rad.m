function az = azimuth_rad(x1, y1, x2, y2)
%AZIMUTH Summary of this function goes here
%   Detailed explanation goes here
delta_x = x2 - x1;
delta_y = y2 - y1;

az = atan2(delta_y, delta_x);

%if az<0
%    az=az+2*pi;
%end

end

