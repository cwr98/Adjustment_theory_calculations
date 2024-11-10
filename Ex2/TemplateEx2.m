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

absFreq = absFrequency(data);
relFreq = relFrequency(data);
cumFreq = cumFrequency(data);

%Absolute frequency function/polygon
figure(1);
%[x y] 
x = 1:length(absFreq);  % Assuming absFreq has one frequency per bin
y = absFreq;
bar(absFreq);
hold on
%plot
bar(x, y);  % Bar plot of absolute frequencies
hold on;
plot(x, y, '-o');  % Line plot on top of bars

title('Absolute Frequency');
xlabel('Data Values');
ylabel('Absolute Frequency');
grid on;

%Relative frequency function/polygon
figure (2);
% [x y] 
x = 1:length(relFreq);
y = relFreq;

bar(relFreq);
hold on
%plot
grid on;
title('Relative frequency')
xlabel('Data Values');
ylabel('Relative frequency function');
hold off


% Cumulative frequency function/polygon
figure (3);
bar(cumFreq);
hold on
% plot
grid on;
title('Cumulative frequency')
xlabel('Data Values') 
ylabel('Cumulative frequency function')
hold off

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









