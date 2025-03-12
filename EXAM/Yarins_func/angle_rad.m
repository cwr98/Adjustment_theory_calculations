function alpha = angle_rad(xi,yi,xk,yk, xl,yl)
%ANGLE Summary of this function goes here
% i - middle, k - right, l - left
%   Detailed explanation goes here
alpha = azimuth_rad(xi,yi,xk,yk)-azimuth_rad(xi,yi,xl,yl);
end

