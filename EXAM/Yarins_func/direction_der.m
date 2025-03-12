function derivative = direction_der(x1, y1, x2, y2, w, input_var)
%DIRECTION_DER Summary of this function goes here
%   Detailed explanation goes here
% Calculate the derivative based on the input_var
switch lower(input_var)
    case 'x1'
        derivative = azimuth_der(x1,y1,x2,y2, 'x1');
    case 'y1'
        derivative = azimuth_der(x1,y1,x2,y2, 'y1');
    case 'x2'
        derivative = azimuth_der(x1,y1,x2,y2, 'x2');
    case 'y2'
        derivative = azimuth_der(x1,y1,x2,y2, 'y2');
    case 'omega'
        derivative = -1;

    otherwise
        error('Invalid input_var. Please use ''x1'', ''y1'', ''x2'', or ''y2''.');
end
end
