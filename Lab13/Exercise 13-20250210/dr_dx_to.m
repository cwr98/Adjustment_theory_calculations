%--------------------------------------------------------------------------
%   
%   Function to calculate the partial derivative for the measured direction
%   to point ''k''
%   Author         : Sven Weisbrich
%   Version        : April 13, 2012
%   Last changes   : April 13, 2012
%
%--------------------------------------------------------------------------

function result=dr_dx_to(yi,xi,yk,xk)
s=distance(yi,xi,yk,xk);
result=-(yk-yi)/s^2;
end