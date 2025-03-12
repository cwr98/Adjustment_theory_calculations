%--------------------------------------------------------------------------
%   
%   Function to calculate the direction from point ''i'' to ''k''
% 
%   Author         : Sven Weisbrich
%   Version        : April 13, 2012
%   Last changes   : April 13, 2012
%
%--------------------------------------------------------------------------

function r=direction(yi, xi, yk, xk, w)
  r=atan2((yk-yi),(xk-xi));
  
  if r<0
      r=r+2*pi;
  end
  r=r-w;
  
end