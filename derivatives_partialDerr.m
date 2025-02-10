clc;
clear all;
close all;

f = @(x) x^3 + 2*x^2 - 5*x + 4; % Define the function
x0 = 2; % Point at which to evaluate the derivative
h = 1e-5; % Small step size
df_numeric = (f(x0 + h) - f(x0 - h)) / (2*h); % Central difference method
disp('Numerical derivative at x = 2:');
disp(df_numeric);
