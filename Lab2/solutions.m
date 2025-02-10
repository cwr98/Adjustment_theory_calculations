%--------------------------------------------------------------------------
%   
%          ADJUSTMENT THEORY I
%     Exercise 2: Fundamentals of statistics  
% 
%   Author         : Anastasia Pasioti
%   Version        : October 04, 2018
%   Last changes   : October 31, 2022
%
%--------------------------------------------------------------------------

clc;
clear all;
close all;

%--------------------------------------------------------------------------
%   Task 2
%--------------------------------------------------------------------------

% 1. find coefficient C
syms x a b
pdf = a*x + b; % Define the PDF formula

% Define the limits for the integral
lower_limit = 2;
upper_limit = 6;

% Compute the definite integral of the PDF between 2 and 6
integral_value = int(pdf, x, lower_limit, upper_limit);

solution = solve(integral_value == 2, a);

% Display the result
disp('The value of a is:');
disp(solution);

% 2. describle the graph of proboblity density function in a analytical way

% 3. what is distribution function F(x) for random var x

% 4. find expectation E(x)

% 5. find STD for random var x


%--------------------------------------------------------------------------
%   Task 3
%--------------------------------------------------------------------------
