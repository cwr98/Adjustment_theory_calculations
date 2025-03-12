function derivative = distance_der(x1, y1, x2, y2, param)
%DISTANCE_DERIVATIVE Calculate the derivative of the distance between two points
%   with respect to either x1, y1, x2, or y2.
%
%   Inputs:
%   - x1, y1: Coordinates of the first point.
%   - x2, y2: Coordinates of the second point.
%   - param: Parameter with respect to which the derivative is calculated.
%            It can be 'x1', 'y1', 'x2', or 'y2'.
%
%   Output:
%   - derivative: The derivative of the distance.

% Calculate the distance between the two points
distance = sqrt((x2 - x1)^2 + (y2 - y1)^2);

% Calculate the derivative with respect to the specified parameter
switch param
    case 'x1'
        derivative = -(x2 - x1) / distance;
    case 'y1'
        derivative = -(y2 - y1) / distance;
    case 'x2'
        derivative = (x2 - x1) / distance;
    case 'y2'
        derivative = (y2 - y1) / distance;
    otherwise
        error('Invalid parameter. Please specify either ''x1'', ''y1'', ''x2'', or ''y2''.')
end
end


