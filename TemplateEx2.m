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
%   Task 3
%--------------------------------------------------------------------------
%Load the data
data = load('distances.txt'); 

%Number of measurements
n = length(data);

%Absolute frequency function/polygon
% figure(1)
% % [x y] 
% % bar
% hold on
% % plot
% grid 
% % title
% % xlabel
% % ylabel
% hold off

%Relative frequency function/polygon
% figure (2)
% [x y] 
% bar
% hold on
% plot
% grid
% title
% xlabel
% ylabel
% hold off
% 
% Cumulative frequency function/polygon
% figure (3)
% [x_cum y_cum] = 
% bar
% hold on
% plot
% grid
% title
% xlabel
% ylabel
% hold off

% Calculate the mean value, variance and standard deviation 
% of a single observation as well as for the arithmetic mean
element = data(1);
%Mean value
mean_value = mean(element);

%Variance of a single observation
var_value = var(element);

%Standard deviation of a single observation
std_value = sqrt(var_value);

%Variance of an arithmetic mean
var_of_mean = var_value / n;

%Standard deviation of an arithmetic mean
std_of_mean = std_value / sqrt(n);

%Question









