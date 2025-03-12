%--------------------------------------------------------------------------
%   
%   Function to calculate the direction from point ''i'' to ''k''
% 
%   Author         : Sven Weisbrich
%   Version        : April 13, 2012
%   Last changes   : April 13, 2012
%
%--------------------------------------------------------------------------

function r=direction_rad(x1, y1, x2, y2, w_rad)
  az=azimuth_rad(x1,y1,x2,y2);

  %if az<0
  %    az=az+2*pi;
  %end
  
  r=az-w_rad;
  
  %if r<0
  %  r=r+2*pi;
  %end
  

  
end