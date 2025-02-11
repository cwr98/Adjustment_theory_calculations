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
  
  r=r-w;

  
end