function angle_derivative = angle_der(xi, yi, xk, yk, xl, yl, input_var)
%ANGLE_DER Compute the derivative of angle with respect to a specified input.
%   angle_derivative = ANGLE_DER(xi, yi, xk, yk, xl, yl, input_var) calculates
%   the derivative of the angle with respect to the specified input_var,
%   where input_var can be 'xi', 'yi', 'xk', 'yk', 'xl', or 'yl'.


% Compute derivative of angle using chain rule
switch lower(input_var)
    case 'xi'
        angle_derivative = azimuth_der(xi, yi, xk, yk, 'x1') - azimuth_der(xi, yi, xl, yl, 'x1');
    case 'yi'
        angle_derivative = azimuth_der(xi, yi, xk, yk, 'y1') - azimuth_der(xi, yi, xl, yl, 'y1');
    case 'xk'
        angle_derivative = azimuth_der(xi, yi, xk, yk, 'x2');
    case 'yk'
        angle_derivative = azimuth_der(xi, yi, xk, yk, 'y2');
    case 'xl'
        angle_derivative = -azimuth_der(xi, yi, xl, yl, 'x2');
    case 'yl'
        angle_derivative = -azimuth_der(xi, yi, xl, yl, 'y2');
    otherwise
        error('Invalid input_var. Please use ''xi'', ''yi'', ''xk'', ''yk'', ''xl'', or ''yl''.');
end
end


