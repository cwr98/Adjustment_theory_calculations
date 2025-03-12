function derivative = azimuth_der(x1, y1, x2, y2, input_var)
%AZIMUTH_DERIVATIVE Calculate the derivative of azimuth with respect to a specified input.
%   d_azimuth_dx = AZIMUTH_DERIVATIVE(x1, y1, x2, y2, input_var) calculates
%   the derivative of the azimuth with respect to the specified input_var,
%   where input_var can be 'x1', 'y1', 'x2', or 'y2'.

% Calculate azimuth
delta_x = x2 - x1;
delta_y = y2 - y1;

% Calculate the derivative based on the input_var
switch lower(input_var)
    case 'x1'
        derivative = (y2 - y1) / (delta_x^2 + delta_y^2);
    case 'y1'
        derivative = (x1 - x2) / (delta_x^2 + delta_y^2);
    case 'x2'
        derivative = (y1 - y2) / (delta_x^2 + delta_y^2);
    case 'y2'
        derivative = (x2 - x1) / (delta_x^2 + delta_y^2);
    otherwise
        error('Invalid input_var. Please use ''x1'', ''y1'', ''x2'', or ''y2''.');
end
end


