clc;
clear all;
close all;

format longG

%stochastic model
s_1 = 4/100; %m 
s_2 = 2/100; %m
s_3 = 1/100; %m
%s_4 = 2/100; %m

p1 = 1/s_1^2
p2 = 1/s_2^2
p3 = 1/s_3^2
%p4 = 1/s_4^2


% Define the symbolic variables
syms x y k

% Define the equations
eq1 = p1*2*(x+y-3)+4*p2*(2*x-y-1.5)+2*(x-y-0.2)*p3+0.9*k==0;
eq2 = p1*2*(x+y-3)-2*p2*(2*x-y-1.5)-2*p3*(x-y-0.2)-k==0;
eq3 = 0.9*x-y==0;

% Solve the system of equations
solution = solve([eq1, eq2, eq3], [x, y, k]);

disp(solution)

