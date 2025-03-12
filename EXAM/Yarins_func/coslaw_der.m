function derivative = cosine_law_derivative(a, b, alpha, param)
%COSINE_LAW_DERIVATIVE Calculate the derivative of the cosine law equation
%   with respect to 'a', 'b', or 'alpha'.
%
%   Inputs:
%   - a: Length of side 'a'.
%   - b: Length of side 'b'.
%   - alpha: Angle between sides 'a' and 'b' in radians.
%   - param: Parameter with respect to which the derivative is calculated.
%            It can be 'a', 'b', or 'alpha'.
%
%   Output:
%   - derivative: The derivative of the cosine law equation.

% Calculate cosine term
cos_alpha = cos(alpha);

% Calculate c
c = sqrt(a^2 + b^2 - 2 * a * b * cos_alpha);

% Calculate derivatives
switch param
    case 'a'
        derivative = (a - b * cos_alpha) / c;
    case 'b'
        derivative = (b - a * cos_alpha) / c;
    case 'alpha'
        derivative = (a * b * sin(alpha)) / c;
    otherwise
        error('Invalid parameter. Please specify either ''a'', ''b'', or ''alpha''.')
end
end


